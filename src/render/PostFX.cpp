#include "stellar/render/PostFX.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace stellar::render {

static const char* kFullscreenVS = R"GLSL(
#version 330 core
layout(location=0) in vec2 aPos;
layout(location=1) in vec2 aUv;
out vec2 vUv;
void main() {
  vUv = aUv;
  gl_Position = vec4(aPos, 0.0, 1.0);
}
)GLSL";

static const char* kBrightFS = R"GLSL(
#version 330 core
in vec2 vUv;
out vec4 oColor;

uniform sampler2D uScene;
uniform float uThreshold;
uniform float uKnee;
uniform float uBoost;

void main() {
  vec3 c = texture(uScene, vUv).rgb;

  // Simple highlight extraction with a soft knee.
  float br = max(max(c.r, c.g), c.b);
  float knee = max(1e-4, uKnee);
  float t = smoothstep(uThreshold, uThreshold + knee, br);

  oColor = vec4(c * t * uBoost, 1.0);
}
)GLSL";

static const char* kBlurFS = R"GLSL(
#version 330 core
in vec2 vUv;
out vec4 oColor;

uniform sampler2D uTex;
uniform vec2 uDir; // texel direction (e.g. vec2(1/w,0) or vec2(0,1/h))

// 9-tap Gaussian (separable, "optimized" offsets)
void main() {
  vec3 sum = vec3(0.0);

  sum += texture(uTex, vUv).rgb * 0.227027;

  sum += texture(uTex, vUv + uDir * 1.384615).rgb * 0.316216;
  sum += texture(uTex, vUv - uDir * 1.384615).rgb * 0.316216;

  sum += texture(uTex, vUv + uDir * 3.230769).rgb * 0.070270;
  sum += texture(uTex, vUv - uDir * 3.230769).rgb * 0.070270;

  oColor = vec4(sum, 1.0);
}
)GLSL";

static const char* kCompositeFS = R"GLSL(
#version 330 core
in vec2 vUv;
out vec4 oColor;

uniform sampler2D uScene;
uniform sampler2D uBloom;

uniform int   uBloomEnabled;
uniform float uBloomIntensity;

uniform float uExposure;
uniform float uGamma;

uniform float uVignette;
uniform float uGrain;
uniform float uAberration;
uniform float uWarp;
uniform float uTime;

// Hyperspace / jump tunnel (screen-space)
uniform float uHyper;
uniform float uHyperTwist;
uniform float uHyperDensity;
uniform float uHyperNoise;
uniform float uHyperIntensity;

// Cheap hash noise
float rand1(vec2 p) {
  // (not a great RNG, but good enough for grain/sparkles)
  return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453123);
}

// ACES-ish filmic curve (common approximation)
vec3 acesFilm(vec3 x) {
  const float a = 2.51;
  const float b = 0.03;
  const float c = 2.43;
  const float d = 0.59;
  const float e = 0.14;
  return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

void main() {
  vec2 uv0 = vUv;
  vec2 uv = uv0;

  vec2 center = vec2(0.5, 0.5);
  vec2 d = uv - center;
  float dist = length(d);
  vec2 dir = (dist > 1e-6) ? (d / dist) : vec2(0.0, 0.0);

  // Preserve stable radial values for vignette + hyperspace patterns (avoid "swimming" vignette).
  float dist0 = dist;
  vec2 dir0 = dir;

  float hyper = clamp(uHyper, 0.0, 1.0);
  float hyperDensity = clamp(uHyperDensity, 0.0, 1.0);
  float hyperNoise = clamp(uHyperNoise, 0.0, 1.0);
  float hyperTwist = clamp(uHyperTwist, 0.0, 1.0);

  // Hyperspace UV distortion: a cheap swirl in polar space, blended back toward the original UV
  // so the edges don't tear apart.
  if (hyper > 0.0) {
    // Max distance to a corner is sqrt(0.5^2 + 0.5^2) ~= 0.7071.
    float rN = clamp(dist / 0.70710678, 0.0, 1.0);
    float ang = atan(d.y, d.x);

    float density = mix(8.0, 42.0, hyperDensity);
    float t = uTime;

    float twist = hyperTwist * hyper;

    // Twist harder near the center; keep edges mostly stable.
    float ang2 = ang + twist * (1.0 - rN) * (0.55 + 0.11 * density) + t * (0.30 + 0.20 * twist);
    vec2 rot = vec2(cos(ang2), sin(ang2)) * dist;

    // Small radial wobble so it doesn't look perfectly mathematical.
    float wobble = sin(t * 1.6 + rN * density) * (0.002 + 0.010 * hyperNoise) * hyper;

    uv = center + rot + dir * wobble;
    uv = mix(uv0, uv, hyper * 0.85);

    // Update radial values for sampling.
    d = uv - center;
    dist = length(d);
    dir = (dist > 1e-6) ? (d / dist) : vec2(0.0, 0.0);
  }

  // Chromatic aberration: offset R/B slightly along radial direction.
  vec3 scene;
  float ca = uAberration;
  if (hyper > 0.0) {
    // Tiny extra CA during hyperspace.
    ca += hyper * 0.0020 * (0.25 + 0.75 * hyperNoise);
  }

  if (ca > 0.0) {
    vec2 o = dir * ca;
    scene.r = texture(uScene, clamp(uv + o, 0.0, 1.0)).r;
    scene.g = texture(uScene, uv).g;
    scene.b = texture(uScene, clamp(uv - o, 0.0, 1.0)).b;
  } else {
    scene = texture(uScene, uv).rgb;
  }

  // Warp streaks: sample a few taps "backwards" along radial direction.
  float warp = uWarp;
  if (hyper > 0.0) {
    // If hyperspace is enabled (manual or auto), ensure at least a small warp streak.
    warp = max(warp, hyper * 0.055);
  }

  if (warp > 0.0) {
    vec3 streak = vec3(0.0);
    const int N = 6;
    for (int i = 1; i <= N; ++i) {
      float t = float(i) / float(N);
      vec2 u = uv - dir * warp * t;
      streak += texture(uScene, clamp(u, 0.0, 1.0)).rgb * (1.0 - t);
    }
    streak /= float(N);
    scene = scene + streak * 0.65;
  }

  vec3 bloom = (uBloomEnabled != 0) ? texture(uBloom, uv).rgb : vec3(0.0);
  vec3 color = scene + bloom * uBloomIntensity;

  // Hyperspace overlay (additive, pre-tonemap).
  if (hyper > 0.0) {
    float rN = clamp(dist0 / 0.70710678, 0.0, 1.0);
    float ang = atan(dir0.y, dir0.x);

    float density = mix(8.0, 42.0, hyperDensity);
    float t = uTime;

    // Moving ring/stripe pattern (strongest near center).
    float phase = (1.0 - rN) * density * 1.25 + t * (7.0 + 4.0 * hyperDensity);
    float rings = 0.5 + 0.5 * sin(phase);
    rings = pow(max(rings, 0.0), 3.0);

    float mask = pow(1.0 - rN, 0.55);

    // Sparkles: occasional bright specks (also strongest near center).
    float n = rand1(uv0 * vec2(1200.0, 900.0) + vec2(t * 0.41, -t * 0.73));
    float spark = smoothstep(0.995 - 0.012 * hyperNoise, 1.0, n);
    spark *= mask * (0.35 + 0.65 * hyperNoise);

    float e = (rings + spark * 1.8) * mask * hyper;

    // Two-color ramp (blue -> magenta) that shifts with time and angle.
    vec3 cA = vec3(0.08, 0.60, 1.35);
    vec3 cB = vec3(1.25, 0.22, 0.95);
    float hue = 0.5 + 0.5 * sin(t * 0.55 + ang * 2.2);
    vec3 hc = mix(cA, cB, hue);

    color += hc * e * max(0.0, uHyperIntensity);
  }

  // Exposure -> tonemap -> vignette -> grain -> gamma
  color *= max(0.0001, uExposure);
  color = acesFilm(color);

  // Vignette
  float vig = 1.0 - uVignette * smoothstep(0.20, 0.95, dist0);
  color *= clamp(vig, 0.0, 1.0);

  // Grain (add in post-tonemap space; subtle)
  float g = uGrain;
  if (g > 0.0) {
    float gn = rand1(uv0 * vec2(1403.1, 911.7) + uTime);
    color += (gn - 0.5) * g;
  }

  // Gamma correction
  float invGamma = 1.0 / max(0.0001, uGamma);
  color = pow(max(color, vec3(0.0)), vec3(invGamma));

  oColor = vec4(color, 1.0);
}
)GLSL";

static unsigned int makeTex2D(int w, int h, int internalFormat, unsigned int format, unsigned int type) {
  unsigned int tex = 0;
  gl::GenTextures(1, &tex);
  gl::BindTexture(GL_TEXTURE_2D, tex);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  gl::TexImage2D(GL_TEXTURE_2D, 0, internalFormat, w, h, 0, format, type, nullptr);
  gl::BindTexture(GL_TEXTURE_2D, 0);
  return tex;
}

PostFX::~PostFX() {
  destroy();
}

bool PostFX::init(std::string* outError) {
  std::string err;

  if (!bright_.build(kFullscreenVS, kBrightFS, &err)) {
    if (outError) *outError = err;
    return false;
  }
  if (!blur_.build(kFullscreenVS, kBlurFS, &err)) {
    if (outError) *outError = err;
    return false;
  }
  if (!composite_.build(kFullscreenVS, kCompositeFS, &err)) {
    if (outError) *outError = err;
    return false;
  }

  // Fullscreen triangle
  const float tri[] = {
    // pos      uv
    -1.0f, -1.0f,  0.0f, 0.0f,
     3.0f, -1.0f,  2.0f, 0.0f,
    -1.0f,  3.0f,  0.0f, 2.0f,
  };

  gl::GenVertexArrays(1, &vao_);
  gl::GenBuffers(1, &vbo_);

  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferData(GL_ARRAY_BUFFER, sizeof(tri), tri, GL_STATIC_DRAW);

  gl::EnableVertexAttribArray(0);
  gl::VertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
  gl::EnableVertexAttribArray(1);
  gl::VertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));

  gl::BindVertexArray(0);
  gl::BindBuffer(GL_ARRAY_BUFFER, 0);

  return true;
}

void PostFX::destroy() {
  if (vbo_) {
    gl::DeleteBuffers(1, &vbo_);
    vbo_ = 0;
  }
  if (vao_) {
    gl::DeleteVertexArrays(1, &vao_);
    vao_ = 0;
  }

  if (sceneTex_) {
    gl::DeleteTextures(1, &sceneTex_);
    sceneTex_ = 0;
  }
  if (brightTex_) {
    gl::DeleteTextures(1, &brightTex_);
    brightTex_ = 0;
  }
  if (pingTex_[0]) gl::DeleteTextures(1, &pingTex_[0]);
  if (pingTex_[1]) gl::DeleteTextures(1, &pingTex_[1]);
  pingTex_[0] = pingTex_[1] = 0;

  if (depthRbo_) {
    gl::DeleteRenderbuffers(1, &depthRbo_);
    depthRbo_ = 0;
  }

  if (sceneFbo_) {
    gl::DeleteFramebuffers(1, &sceneFbo_);
    sceneFbo_ = 0;
  }
  if (brightFbo_) {
    gl::DeleteFramebuffers(1, &brightFbo_);
    brightFbo_ = 0;
  }
  if (pingFbo_[0]) gl::DeleteFramebuffers(1, &pingFbo_[0]);
  if (pingFbo_[1]) gl::DeleteFramebuffers(1, &pingFbo_[1]);
  pingFbo_[0] = pingFbo_[1] = 0;

  w_ = h_ = bloomW_ = bloomH_ = 0;
}

void PostFX::ensureSize(int w, int h) {
  if (w <= 0 || h <= 0) return;
  if (w == w_ && h == h_) return;
  createOrResize(w, h);
}

void PostFX::createOrResize(int w, int h) {
  // Delete only the size-dependent pieces.
  if (sceneTex_) { gl::DeleteTextures(1, &sceneTex_); sceneTex_ = 0; }
  if (brightTex_) { gl::DeleteTextures(1, &brightTex_); brightTex_ = 0; }
  if (pingTex_[0]) gl::DeleteTextures(1, &pingTex_[0]);
  if (pingTex_[1]) gl::DeleteTextures(1, &pingTex_[1]);
  pingTex_[0] = pingTex_[1] = 0;

  if (depthRbo_) { gl::DeleteRenderbuffers(1, &depthRbo_); depthRbo_ = 0; }

  if (!sceneFbo_) gl::GenFramebuffers(1, &sceneFbo_);
  if (!brightFbo_) gl::GenFramebuffers(1, &brightFbo_);
  if (!pingFbo_[0]) gl::GenFramebuffers(1, &pingFbo_[0]);
  if (!pingFbo_[1]) gl::GenFramebuffers(1, &pingFbo_[1]);

  w_ = w;
  h_ = h;
  bloomW_ = std::max(1, w / 2);
  bloomH_ = std::max(1, h / 2);

  // HDR scene texture + depth
  sceneTex_ = makeTex2D(w_, h_, GL_RGBA16F, GL_RGBA, GL_FLOAT);

  gl::GenRenderbuffers(1, &depthRbo_);
  gl::BindRenderbuffer(GL_RENDERBUFFER, depthRbo_);
  gl::RenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w_, h_);
  gl::BindRenderbuffer(GL_RENDERBUFFER, 0);

  gl::BindFramebuffer(GL_FRAMEBUFFER, sceneFbo_);
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, sceneTex_, 0);
  gl::FramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthRbo_);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);

  if (gl::CheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
    // If incomplete, keep going but the scene will render black; better than crash.
  }

  // Bright-pass texture
  brightTex_ = makeTex2D(bloomW_, bloomH_, GL_RGBA16F, GL_RGBA, GL_FLOAT);

  gl::BindFramebuffer(GL_FRAMEBUFFER, brightFbo_);
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, brightTex_, 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);

  // Ping-pong textures
  pingTex_[0] = makeTex2D(bloomW_, bloomH_, GL_RGBA16F, GL_RGBA, GL_FLOAT);
  pingTex_[1] = makeTex2D(bloomW_, bloomH_, GL_RGBA16F, GL_RGBA, GL_FLOAT);

  gl::BindFramebuffer(GL_FRAMEBUFFER, pingFbo_[0]);
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pingTex_[0], 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);

  gl::BindFramebuffer(GL_FRAMEBUFFER, pingFbo_[1]);
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, pingTex_[1], 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);

  gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
}

void PostFX::beginScene(int w, int h) const {
  (void)w;
  (void)h;
  gl::BindFramebuffer(GL_FRAMEBUFFER, sceneFbo_);
}

void PostFX::drawFullscreen() const {
  gl::BindVertexArray(vao_);
  glDrawArrays(GL_TRIANGLES, 0, 3);
  gl::BindVertexArray(0);
}

void PostFX::present(int w, int h, const PostFXSettings& s, float timeSeconds) {
  if (!sceneFbo_ || !sceneTex_) return;

  ensureSize(w, h);

  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);
  glDisable(GL_BLEND);

  // 1) Bright pass
  if (s.bloomEnabled) {
    gl::BindFramebuffer(GL_FRAMEBUFFER, brightFbo_);
    glViewport(0, 0, bloomW_, bloomH_);
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);

    bright_.bind();
    gl::ActiveTexture(GL_TEXTURE0);
    gl::BindTexture(GL_TEXTURE_2D, sceneTex_);
    bright_.setUniform1i("uScene", 0);
    bright_.setUniform1f("uThreshold", s.bloomThreshold);
    bright_.setUniform1f("uKnee", s.bloomKnee);
    bright_.setUniform1f("uBoost", s.bloomBoost);
    drawFullscreen();
  }

  // 2) Blur ping-pong
  unsigned int bloomTex = brightTex_;
  if (s.bloomEnabled && s.bloomPasses > 0) {
    blur_.bind();
    blur_.setUniform1i("uTex", 0);

    bool horizontal = true;
    unsigned int input = brightTex_;

    for (int i = 0; i < s.bloomPasses; ++i) {
      const int idx = horizontal ? 0 : 1;
      gl::BindFramebuffer(GL_FRAMEBUFFER, pingFbo_[idx]);
      glViewport(0, 0, bloomW_, bloomH_);
      glClearColor(0, 0, 0, 1);
      glClear(GL_COLOR_BUFFER_BIT);

      gl::ActiveTexture(GL_TEXTURE0);
      gl::BindTexture(GL_TEXTURE_2D, input);

      if (horizontal) {
        blur_.setUniform2f("uDir", 1.0f / (float)bloomW_, 0.0f);
      } else {
        blur_.setUniform2f("uDir", 0.0f, 1.0f / (float)bloomH_);
      }

      drawFullscreen();

      input = pingTex_[idx];
      bloomTex = input;
      horizontal = !horizontal;
    }
  }

  // 3) Composite
  gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
  glViewport(0, 0, w, h);

  composite_.bind();
  gl::ActiveTexture(GL_TEXTURE0);
  gl::BindTexture(GL_TEXTURE_2D, sceneTex_);
  composite_.setUniform1i("uScene", 0);

  gl::ActiveTexture(GL_TEXTURE1);
  gl::BindTexture(GL_TEXTURE_2D, bloomTex);
  composite_.setUniform1i("uBloom", 1);

  composite_.setUniform1i("uBloomEnabled", s.bloomEnabled ? 1 : 0);
  composite_.setUniform1f("uBloomIntensity", s.bloomIntensity);

  composite_.setUniform1f("uExposure", s.exposure);
  composite_.setUniform1f("uGamma", s.gamma);

  composite_.setUniform1f("uVignette", s.vignette);
  composite_.setUniform1f("uGrain", s.grain);
  composite_.setUniform1f("uAberration", s.chromaticAberration);
  composite_.setUniform1f("uWarp", s.warp);
  composite_.setUniform1f("uTime", timeSeconds);

  composite_.setUniform1f("uHyper", s.hyperspace);
  composite_.setUniform1f("uHyperTwist", s.hyperspaceTwist);
  composite_.setUniform1f("uHyperDensity", s.hyperspaceDensity);
  composite_.setUniform1f("uHyperNoise", s.hyperspaceNoise);
  composite_.setUniform1f("uHyperIntensity", s.hyperspaceIntensity);

  drawFullscreen();

  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);
}

} // namespace stellar::render
