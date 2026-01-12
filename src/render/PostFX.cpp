#include "stellar/render/PostFX.h"

#include "stellar/render/AutoExposure.h"
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

// Luminance extraction: output RG where
//  R = log(luminance) * weight
//  G = weight
//
// This allows weighted log-average metering via downsampling:
//    avgLogLum = avg(R) / avg(G)
// and avgLum = exp(avgLogLum).
static const char* kLuminanceFS = R"GLSL(
#version 330 core
in vec2 vUv;
out vec4 oColor;

uniform sampler2D uScene;
uniform float uAspect;
uniform float uCenterWeight; // 0..1

void main() {
  vec3 c = texture(uScene, vUv).rgb;
  float lum = dot(c, vec3(0.2126, 0.7152, 0.0722));
  lum = max(lum, 1e-4);

  float logLum = log(lum);

  // Optional center-weighted metering (reduces bright edge influence).
  float w = 1.0;
  float cw = clamp(uCenterWeight, 0.0, 1.0);
  if (cw > 0.0) {
    vec2 p = (vUv - 0.5) * vec2(max(0.001, uAspect), 1.0);
    float d2 = dot(p, p);
    // Gaussian-ish falloff (sigma tuned for "spot/center" metering feel).
    float g = exp(-d2 * 6.5);
    w = mix(1.0, g, cw);
  }

  oColor = vec4(logLum * w, w, 0.0, 1.0);
}
)GLSL";

// 2x2 downsample by averaging, implemented via texelFetch (no filtering).
// Works for any input size; edge texels are clamped.
static const char* kReduce2x2FS = R"GLSL(
#version 330 core
in vec2 vUv;
out vec4 oColor;

uniform sampler2D uIn;
uniform vec2 uInSize;

ivec2 inSizeI() {
  return max(ivec2(1, 1), ivec2(uInSize + 0.5));
}

ivec2 clampCoord(ivec2 c) {
  ivec2 s = inSizeI();
  return ivec2(clamp(c.x, 0, s.x - 1), clamp(c.y, 0, s.y - 1));
}

void main() {
  ivec2 outPx = ivec2(gl_FragCoord.xy);
  ivec2 base = outPx * 2;

  vec4 a = texelFetch(uIn, clampCoord(base + ivec2(0, 0)), 0);
  vec4 b = texelFetch(uIn, clampCoord(base + ivec2(1, 0)), 0);
  vec4 c = texelFetch(uIn, clampCoord(base + ivec2(0, 1)), 0);
  vec4 d = texelFetch(uIn, clampCoord(base + ivec2(1, 1)), 0);

  oColor = (a + b + c + d) * 0.25;
}
)GLSL";

static const char* kCompositeFS = R"GLSL(
#version 330 core

// -----------------------------------------------------------------------------
// Procedural compositor shader variants
// -----------------------------------------------------------------------------
//
// PostFX can inject additional `#define` lines after the #version directive at
// runtime (see PostFX::buildCompositeVariant). These defines allow the game to
// build deterministic shader permutations from a seed.

#ifndef STELLAR_TONEMAP_MODE
#define STELLAR_TONEMAP_MODE 0
#endif

#ifndef STELLAR_GRADE_ENABLED
#define STELLAR_GRADE_ENABLED 0
#endif

#ifndef STELLAR_GRADE_SAT
#define STELLAR_GRADE_SAT 1.0
#endif
#ifndef STELLAR_GRADE_CONTRAST
#define STELLAR_GRADE_CONTRAST 1.0
#endif
#ifndef STELLAR_GRADE_HUE
#define STELLAR_GRADE_HUE 0.0
#endif

#ifndef STELLAR_GRADE_LIFT_R
#define STELLAR_GRADE_LIFT_R 0.0
#endif
#ifndef STELLAR_GRADE_LIFT_G
#define STELLAR_GRADE_LIFT_G 0.0
#endif
#ifndef STELLAR_GRADE_LIFT_B
#define STELLAR_GRADE_LIFT_B 0.0
#endif

#ifndef STELLAR_GRADE_GAMMA_R
#define STELLAR_GRADE_GAMMA_R 1.0
#endif
#ifndef STELLAR_GRADE_GAMMA_G
#define STELLAR_GRADE_GAMMA_G 1.0
#endif
#ifndef STELLAR_GRADE_GAMMA_B
#define STELLAR_GRADE_GAMMA_B 1.0
#endif

#ifndef STELLAR_GRADE_GAIN_R
#define STELLAR_GRADE_GAIN_R 1.0
#endif
#ifndef STELLAR_GRADE_GAIN_G
#define STELLAR_GRADE_GAIN_G 1.0
#endif
#ifndef STELLAR_GRADE_GAIN_B
#define STELLAR_GRADE_GAIN_B 1.0
#endif

#ifndef STELLAR_SPLIT_TONE_ENABLED
#define STELLAR_SPLIT_TONE_ENABLED 0
#endif
#ifndef STELLAR_SPLIT_TONE_STRENGTH
#define STELLAR_SPLIT_TONE_STRENGTH 0.0
#endif
#ifndef STELLAR_SPLIT_SHADOW_R
#define STELLAR_SPLIT_SHADOW_R 1.0
#endif
#ifndef STELLAR_SPLIT_SHADOW_G
#define STELLAR_SPLIT_SHADOW_G 1.0
#endif
#ifndef STELLAR_SPLIT_SHADOW_B
#define STELLAR_SPLIT_SHADOW_B 1.0
#endif
#ifndef STELLAR_SPLIT_HI_R
#define STELLAR_SPLIT_HI_R 1.0
#endif
#ifndef STELLAR_SPLIT_HI_G
#define STELLAR_SPLIT_HI_G 1.0
#endif
#ifndef STELLAR_SPLIT_HI_B
#define STELLAR_SPLIT_HI_B 1.0
#endif

#ifndef STELLAR_VIGNETTE_MUL
#define STELLAR_VIGNETTE_MUL 1.0
#endif
#ifndef STELLAR_GRAIN_MUL
#define STELLAR_GRAIN_MUL 1.0
#endif

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

// Viewport size in pixels.
uniform vec2 uResolution;

// Retro / CRT compositor mode.
uniform int   uRetroEnabled;
uniform int   uRetroPixelSize;
uniform int   uRetroSteps;
uniform float uRetroDither;
uniform float uRetroScanlines;
uniform float uRetroCurvature;
uniform float uRetroJitter;

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

// 8x8 Bayer matrix (0..63). Used for ordered dithering.
// Reference pattern: https://en.wikipedia.org/wiki/Ordered_dithering (Bayer matrix)
const int kBayer8[64] = int[64](
  0, 48, 12, 60,  3, 51, 15, 63,
 32, 16, 44, 28, 35, 19, 47, 31,
  8, 56,  4, 52, 11, 59,  7, 55,
 40, 24, 36, 20, 43, 27, 39, 23,
  2, 50, 14, 62,  1, 49, 13, 61,
 34, 18, 46, 30, 33, 17, 45, 29,
 10, 58,  6, 54,  9, 57,  5, 53,
 42, 26, 38, 22, 41, 25, 37, 21
);

float bayer8x8(ivec2 p) {
  int x = p.x & 7;
  int y = p.y & 7;
  int idx = y * 8 + x;
  return (float(kBayer8[idx]) + 0.5) / 64.0;
}

// Cheap CRT barrel distortion.
vec2 crtWarp(vec2 uv, float k) {
  // Map to -1..1.
  vec2 p = uv * 2.0 - 1.0;
  float r2 = dot(p, p);
  // Push corners outward.
  p += p * r2 * (0.18 * k);
  return p * 0.5 + 0.5;
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

// Reinhard tonemap (simple, robust).
vec3 reinhardTone(vec3 x) {
  return x / (1.0 + x);
}

// Uncharted 2 filmic tonemap (Hable curve).
// (Implemented as a common approximation; white point fixed at W=11.2.)
vec3 uncharted2Tonemap(vec3 x) {
  const float A = 0.15;
  const float B = 0.50;
  const float C = 0.10;
  const float D = 0.20;
  const float E = 0.02;
  const float F = 0.30;
  return ((x * (A * x + C * B) + D * E) / (x * (A * x + B) + D * F)) - E / F;
}

vec3 uncharted2Filmic(vec3 x) {
  const float W = 11.2;
  vec3 curr = uncharted2Tonemap(x);
  vec3 whiteScale = 1.0 / uncharted2Tonemap(vec3(W));
  return clamp(curr * whiteScale, 0.0, 1.0);
}

// Compile-time tonemap selection (driven by injected #defines).
vec3 tonemap(vec3 x) {
#if STELLAR_TONEMAP_MODE == 0
  return acesFilm(x);
#elif STELLAR_TONEMAP_MODE == 1
  return reinhardTone(x);
#elif STELLAR_TONEMAP_MODE == 2
  return uncharted2Filmic(x);
#else
  return clamp(x, 0.0, 1.0);
#endif
}

// Approximate hue rotation by rotating around the neutral axis in RGB space.
// This is not a strict HSV hue rotation, but it is cheap and "looks right" for
// stylized palettes.
vec3 hueRotate(vec3 c, float angle) {
  vec3 axis = normalize(vec3(1.0, 1.0, 1.0));
  float ca = cos(angle);
  float sa = sin(angle);
  // Rodrigues' rotation formula.
  return c * ca + cross(axis, c) * sa + axis * dot(axis, c) * (1.0 - ca);
}

vec3 applyGrade(vec3 c) {
#if STELLAR_GRADE_ENABLED == 0
  return c;
#else
  const vec3 lumaW = vec3(0.2126, 0.7152, 0.0722);

  vec3 lift = vec3(STELLAR_GRADE_LIFT_R, STELLAR_GRADE_LIFT_G, STELLAR_GRADE_LIFT_B);
  vec3 gamma = vec3(STELLAR_GRADE_GAMMA_R, STELLAR_GRADE_GAMMA_G, STELLAR_GRADE_GAMMA_B);
  vec3 gain = vec3(STELLAR_GRADE_GAIN_R, STELLAR_GRADE_GAIN_G, STELLAR_GRADE_GAIN_B);

  // Lift / Gamma / Gain (applied in linear-ish space after tonemap).
  c = c * gain + lift;
  c = max(c, vec3(0.0));
  c = pow(c, vec3(1.0) / max(gamma, vec3(0.0001)));

  // Saturation (lerp to luma).
  float y = dot(c, lumaW);
  c = mix(vec3(y), c, STELLAR_GRADE_SAT);

  // Contrast around mid-grey.
  c = (c - 0.5) * STELLAR_GRADE_CONTRAST + 0.5;

  // Hue shift (optional).
  float h = STELLAR_GRADE_HUE;
  if (abs(h) > 1e-6) {
    c = hueRotate(c, h);
  }

  // Split tone (optional).
#if STELLAR_SPLIT_TONE_ENABLED == 1
  float l = dot(c, lumaW);
  vec3 sh = c * vec3(STELLAR_SPLIT_SHADOW_R, STELLAR_SPLIT_SHADOW_G, STELLAR_SPLIT_SHADOW_B);
  vec3 hi = c * vec3(STELLAR_SPLIT_HI_R, STELLAR_SPLIT_HI_G, STELLAR_SPLIT_HI_B);
  float t = smoothstep(0.15, 0.85, l);
  vec3 st = mix(sh, hi, t);
  c = mix(c, st, clamp(STELLAR_SPLIT_TONE_STRENGTH, 0.0, 1.0));
#endif

  return clamp(c, 0.0, 1.0);
#endif
}

void main() {
  vec2 uv0 = vUv;

  // ---------------------------------------------------------------------------
  // Retro / CRT compositor path (optional)
  // ---------------------------------------------------------------------------
  if (uRetroEnabled != 0) {
    vec2 uv = uv0;
    const vec2 center = vec2(0.5, 0.5);

    // CRT curvature (barrel distortion)
    float curv = clamp(uRetroCurvature, 0.0, 1.0);
    if (curv > 0.0) {
      uv = crtWarp(uv, curv);
    }

    // Horizontal per-scanline jitter (very subtle; VHS-ish)
    float jit = clamp(uRetroJitter, 0.0, 1.0);
    if (jit > 0.0) {
      float line = floor(gl_FragCoord.y);
      float t = floor(uTime * 60.0);
      float j = (rand1(vec2(line, t)) - 0.5) * 2.0;
      uv.x += j * 0.008 * jit;
    }

    uv = clamp(uv, 0.0, 1.0);

    // Pixelation grid and exact sampling (texelFetch bypasses bilinear filtering).
    vec2 resF = max(uResolution, vec2(1.0));
    ivec2 res = ivec2(resF);

    int px = max(1, uRetroPixelSize);
    ivec2 grid = ivec2(max(1, res.x / px), max(1, res.y / px));
    vec2 g = vec2(grid);
    ivec2 ip = ivec2(clamp(floor(uv * g), vec2(0.0), g - 1.0));

    ivec2 tex = ip * px + ivec2(px / 2);
    tex = clamp(tex, ivec2(0), res - 1);

    vec2 uvS = (vec2(tex) + 0.5) / vec2(res);

    vec3 scene = texelFetch(uScene, tex, 0).rgb;
    vec3 bloom = (uBloomEnabled != 0) ? texture(uBloom, uvS).rgb : vec3(0.0);
    vec3 color = scene + bloom * uBloomIntensity;

    // Exposure -> tonemap
    color *= max(0.0001, uExposure);
    color = tonemap(color);
    color = applyGrade(color);

    // Optional vignette (based on warped UV).
    float dist = length(uv - center);
    float vig = 1.0 - (uVignette * STELLAR_VIGNETTE_MUL) * smoothstep(0.20, 0.95, dist);
    color *= clamp(vig, 0.0, 1.0);

    // Grain (pre-quantization)
    float gAmt = uGrain * STELLAR_GRAIN_MUL;
    if (gAmt > 0.0) {
      float gn = rand1(uvS * vec2(1403.1, 911.7) + uTime);
      color += (gn - 0.5) * gAmt;
    }

    // Ordered dithering + quantization
    int stepsI = max(2, uRetroSteps);
    float steps = float(stepsI);
    color = clamp(color, 0.0, 1.0);

    float d = 0.0;
    float dith = clamp(uRetroDither, 0.0, 1.0);
    if (dith > 0.0) {
      float th = bayer8x8(ip);
      d = (th - 0.5) * dith;
    }
    color = floor(color * steps + 0.5 + d) / steps;

    // Scanlines (darken every other line)
    float sc = clamp(uRetroScanlines, 0.0, 1.0);
    if (sc > 0.0) {
      float m = mod(gl_FragCoord.y, 2.0);
      float dark = 1.0 - 0.18 * sc;
      color *= mix(dark, 1.0, m);
    }

    // Gamma correction
    float invGamma = 1.0 / max(0.0001, uGamma);
    color = pow(max(color, vec3(0.0)), vec3(invGamma));

    oColor = vec4(color, 1.0);
    return;
  }

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
  color = tonemap(color);
  color = applyGrade(color);

  // Vignette
  float vig = 1.0 - (uVignette * STELLAR_VIGNETTE_MUL) * smoothstep(0.20, 0.95, dist0);
  color *= clamp(vig, 0.0, 1.0);

  // Grain (add in post-tonemap space; subtle)
  float g = uGrain * STELLAR_GRAIN_MUL;
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
  if (!luminance_.build(kFullscreenVS, kLuminanceFS, &err)) {
    if (outError) *outError = err;
    return false;
  }
  if (!reduce_.build(kFullscreenVS, kReduce2x2FS, &err)) {
    if (outError) *outError = err;
    return false;
  }
  if (!composite_.build(kFullscreenVS, kCompositeFS, &err)) {
    if (outError) *outError = err;
    return false;
  }

  hasCustomComposite_ = false;
  compositeFsBuilt_ = kCompositeFS;

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

  // 1x1 black texture (HDR) used as a stable fallback input.
  if (!blackTex_) {
    blackTex_ = makeTex2D(1, 1, GL_RGBA16F, GL_RGBA, GL_FLOAT);
    const float z[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    gl::BindTexture(GL_TEXTURE_2D, blackTex_);
    gl::TexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 1, 1, GL_RGBA, GL_FLOAT, z);
    gl::BindTexture(GL_TEXTURE_2D, 0);
  }

  return true;
}

namespace {

static std::string injectDefines(std::string_view baseFs, std::string_view defines) {
  // Insert defines right after the first newline (after '#version ...').
  std::string fs(baseFs);
  const std::size_t nl = fs.find('\n');
  const std::size_t insertPos = (nl == std::string::npos) ? fs.size() : (nl + 1);

  std::string block;
  block.reserve(defines.size() + 64);
  block += "\n// --- injected variant defines ---\n";
  block.append(defines.data(), defines.size());
  if (!block.empty() && block.back() != '\n') block.push_back('\n');
  block += "// --- end injected defines ---\n\n";

  fs.insert(insertPos, block);
  return fs;
}

} // namespace

bool PostFX::buildCompositeVariant(std::string_view defines, std::string* outError) {
  std::string err;
  ShaderProgram tmp;
  const std::string fs = injectDefines(kCompositeFS, defines);

  if (!tmp.build(kFullscreenVS, fs, &err)) {
    if (outError) *outError = err;
    return false;
  }

  composite_ = std::move(tmp);
  hasCustomComposite_ = true;
  compositeFsBuilt_ = fs;
  return true;
}

bool PostFX::resetCompositeShader(std::string* outError) {
  std::string err;
  ShaderProgram tmp;
  if (!tmp.build(kFullscreenVS, kCompositeFS, &err)) {
    if (outError) *outError = err;
    return false;
  }
  composite_ = std::move(tmp);
  hasCustomComposite_ = false;
  compositeFsBuilt_ = kCompositeFS;
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

  if (blackTex_) {
    gl::DeleteTextures(1, &blackTex_);
    blackTex_ = 0;
  }

  if (depthRbo_) {
    gl::DeleteRenderbuffers(1, &depthRbo_);
    depthRbo_ = 0;
  }

  if (sceneFbo_) {
    gl::DeleteFramebuffers(1, &sceneFbo_);
    sceneFbo_ = 0;
  }

  w_ = h_ = 0;
}

void PostFX::ensureSize(int w, int h) {
  if (w <= 0 || h <= 0) return;
  if (w == w_ && h == h_) return;
  createOrResize(w, h);
}

void PostFX::createOrResize(int w, int h) {
  // Delete only the size-dependent pieces.
  if (sceneTex_) { gl::DeleteTextures(1, &sceneTex_); sceneTex_ = 0; }

  if (depthRbo_) { gl::DeleteRenderbuffers(1, &depthRbo_); depthRbo_ = 0; }

  if (!sceneFbo_) gl::GenFramebuffers(1, &sceneFbo_);
  w_ = w;
  h_ = h;

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

  // Track a dt for temporal effects (auto exposure adaptation).
  float dt = 0.0f;
  if (haveLastPresentTime_) {
    dt = timeSeconds - lastPresentTime_;
  }
  lastPresentTime_ = timeSeconds;
  haveLastPresentTime_ = true;

  // ---------------------------------------------------------------------------
  // FrameGraph-driven post-processing
  // ---------------------------------------------------------------------------
  frameGraph_.reset();
  frameGraph_.setViewport(w, h);

  // Manual exposure can be combined with the previous frame's auto-exposure.
  // The current frame's luminance is measured and applied on the NEXT frame.
  const float effectiveExposure = s.exposure * ((s.autoExposureEnabled) ? std::max(0.0001f, autoExposure_.exposure) : 1.0f);

  // External inputs.
  const auto scene = frameGraph_.importTexture("SceneHDR", sceneTex_, w_, h_);
  const auto black = frameGraph_.importTexture("Black", blackTex_, 1, 1);

  // --- Auto exposure: measure scene luminance (log-average) -----------------
  FrameGraph::TextureHandle lumFinal{};
  if (s.autoExposureEnabled) {
    // Start from a capped quarter-res buffer for performance.
    int lumW = std::max(1, w_ / 4);
    int lumH = std::max(1, h_ / 4);

    const int cap = std::max(16, s.autoExposureMaxSize);
    if (lumW > cap || lumH > cap) {
      const float fx = (lumW > 0) ? (float)cap / (float)lumW : 1.0f;
      const float fy = (lumH > 0) ? (float)cap / (float)lumH : 1.0f;
      const float f = std::min(fx, fy);
      lumW = std::max(1, (int)std::floor((float)lumW * f));
      lumH = std::max(1, (int)std::floor((float)lumH * f));
    }

    FrameGraph::TextureDesc lumD{};
    lumD.width = lumW;
    lumD.height = lumH;
    lumD.internalFormat = GL_RG16F;
    lumD.format = GL_RG;
    lumD.type = GL_FLOAT;
    lumD.linearFilter = false;
    lumD.clampToEdge = true;

    const auto lum0 = frameGraph_.createTexture("AE_Lum0", lumD);

    frameGraph_.addPass(
        "AE_Extract",
        {scene},
        lum0,
        [&](const FrameGraph::PassContext& ctx) {
          (void)ctx;
          luminance_.bind();
          gl::ActiveTexture(GL_TEXTURE0);
          gl::BindTexture(GL_TEXTURE_2D, ctx.texture(scene));
          luminance_.setUniform1i("uScene", 0);
          const float asp = (h_ > 0) ? ((float)w_ / (float)h_) : 1.0f;
          luminance_.setUniform1f("uAspect", asp);
          luminance_.setUniform1f("uCenterWeight", s.autoExposureCenterWeight);
          drawFullscreen();
        });

    FrameGraph::TextureHandle prev = lum0;
    int curW = lumW;
    int curH = lumH;

    // Downsample chain until 1x1.
    int level = 0;
    while (curW > 1 || curH > 1) {
      const int nxtW = std::max(1, (curW + 1) / 2);
      const int nxtH = std::max(1, (curH + 1) / 2);

      FrameGraph::TextureDesc d = lumD;
      d.width = nxtW;
      d.height = nxtH;

      const auto out = frameGraph_.createTexture("AE_Reduce_" + std::to_string(level), d);
      frameGraph_.addPass(
          "AE_Reduce_" + std::to_string(level),
          {prev},
          out,
          [&, prev, out](const FrameGraph::PassContext& ctx) {
            (void)out;
            reduce_.bind();
            gl::ActiveTexture(GL_TEXTURE0);
            gl::BindTexture(GL_TEXTURE_2D, ctx.texture(prev));
            reduce_.setUniform1i("uIn", 0);
            reduce_.setUniform2f("uInSize", (float)ctx.textureWidth(prev), (float)ctx.textureHeight(prev));
            drawFullscreen();
          });

      prev = out;
      curW = nxtW;
      curH = nxtH;
      ++level;
      if (level > 16) break; // safety
    }

    lumFinal = prev;
  }

  FrameGraph::TextureHandle bloomTex = black;

  // Half-res HDR transient textures.
  FrameGraph::TextureDesc half{};
  half.scale = 0.5f;
  half.internalFormat = GL_RGBA16F;
  half.format = GL_RGBA;
  half.type = GL_FLOAT;
  half.linearFilter = true;
  half.clampToEdge = true;

  if (s.bloomEnabled) {
    const auto brightOut = frameGraph_.createTexture("Bright", half);

    frameGraph_.addPass(
      "BrightPass",
      {scene},
      brightOut,
      [&](const FrameGraph::PassContext& ctx) {
        (void)ctx;
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT);

        bright_.bind();
        gl::ActiveTexture(GL_TEXTURE0);
        gl::BindTexture(GL_TEXTURE_2D, ctx.texture(scene));
        bright_.setUniform1i("uScene", 0);
        bright_.setUniform1f("uThreshold", s.bloomThreshold);
        bright_.setUniform1f("uKnee", s.bloomKnee);
        bright_.setUniform1f("uBoost", s.bloomBoost);
        drawFullscreen();
      });

    FrameGraph::TextureHandle prev = brightOut;
    bool horizontal = true;
    const int passes = std::max(0, s.bloomPasses);

    for (int i = 0; i < passes; ++i) {
      const bool hdir = horizontal;
      const std::string pname = hdir ? ("BlurH_" + std::to_string(i)) : ("BlurV_" + std::to_string(i));
      const auto out = frameGraph_.createTexture(pname, half);

      frameGraph_.addPass(
        pname,
        {prev},
        out,
        [&, hdir, prev, out](const FrameGraph::PassContext& ctx) {
          (void)out;
          glClearColor(0, 0, 0, 1);
          glClear(GL_COLOR_BUFFER_BIT);

          blur_.bind();
          blur_.setUniform1i("uTex", 0);

          gl::ActiveTexture(GL_TEXTURE0);
          gl::BindTexture(GL_TEXTURE_2D, ctx.texture(prev));

          const float invW = (ctx.outW > 0) ? (1.0f / (float)ctx.outW) : 0.0f;
          const float invH = (ctx.outH > 0) ? (1.0f / (float)ctx.outH) : 0.0f;
          blur_.setUniform2f("uDir", hdir ? invW : 0.0f, hdir ? 0.0f : invH);
          drawFullscreen();
        });

      prev = out;
      horizontal = !horizontal;
    }

    bloomTex = prev;
  }

  frameGraph_.addPass(
    "Composite",
    {scene, bloomTex},
    FrameGraph::Backbuffer(),
    [&](const FrameGraph::PassContext& ctx) {
      (void)ctx;
      composite_.bind();

      gl::ActiveTexture(GL_TEXTURE0);
      gl::BindTexture(GL_TEXTURE_2D, ctx.texture(scene));
      composite_.setUniform1i("uScene", 0);

      gl::ActiveTexture(GL_TEXTURE1);
      gl::BindTexture(GL_TEXTURE_2D, ctx.texture(bloomTex));
      composite_.setUniform1i("uBloom", 1);

      composite_.setUniform1i("uBloomEnabled", s.bloomEnabled ? 1 : 0);
      composite_.setUniform1f("uBloomIntensity", s.bloomIntensity);

      composite_.setUniform1f("uExposure", effectiveExposure);
      composite_.setUniform1f("uGamma", s.gamma);

      composite_.setUniform1f("uVignette", s.vignette);
      composite_.setUniform1f("uGrain", s.grain);
      composite_.setUniform1f("uAberration", s.chromaticAberration);
      composite_.setUniform1f("uWarp", s.warp);
      composite_.setUniform1f("uTime", timeSeconds);

      composite_.setUniform2f("uResolution", (float)w, (float)h);

      // Retro compositor
      composite_.setUniform1i("uRetroEnabled", s.retroEnabled ? 1 : 0);
      composite_.setUniform1i("uRetroPixelSize", s.retroPixelSize);
      composite_.setUniform1i("uRetroSteps", s.retroColorSteps);
      composite_.setUniform1f("uRetroDither", s.retroDitherStrength);
      composite_.setUniform1f("uRetroScanlines", s.retroScanlines);
      composite_.setUniform1f("uRetroCurvature", s.retroCurvature);
      composite_.setUniform1f("uRetroJitter", s.retroJitter);

      composite_.setUniform1f("uHyper", s.hyperspace);
      composite_.setUniform1f("uHyperTwist", s.hyperspaceTwist);
      composite_.setUniform1f("uHyperDensity", s.hyperspaceDensity);
      composite_.setUniform1f("uHyperNoise", s.hyperspaceNoise);
      composite_.setUniform1f("uHyperIntensity", s.hyperspaceIntensity);

      drawFullscreen();
    });

  std::string gErr;
  if (!frameGraph_.compile(&gErr)) {
    // Fallback: if the graph fails to compile, at least draw the raw scene.
    gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, w, h);

    composite_.bind();
    gl::ActiveTexture(GL_TEXTURE0);
    gl::BindTexture(GL_TEXTURE_2D, sceneTex_);
    composite_.setUniform1i("uScene", 0);

    gl::ActiveTexture(GL_TEXTURE1);
    gl::BindTexture(GL_TEXTURE_2D, blackTex_);
    composite_.setUniform1i("uBloom", 1);

    composite_.setUniform1i("uBloomEnabled", 0);
    composite_.setUniform1f("uBloomIntensity", 0.0f);

    composite_.setUniform1f("uExposure", effectiveExposure);
    composite_.setUniform1f("uGamma", s.gamma);
    composite_.setUniform1f("uVignette", s.vignette);
    composite_.setUniform1f("uGrain", s.grain);
    composite_.setUniform1f("uAberration", s.chromaticAberration);
    composite_.setUniform1f("uWarp", s.warp);
    composite_.setUniform1f("uTime", timeSeconds);
    composite_.setUniform2f("uResolution", (float)w, (float)h);
    composite_.setUniform1i("uRetroEnabled", s.retroEnabled ? 1 : 0);
    composite_.setUniform1i("uRetroPixelSize", s.retroPixelSize);
    composite_.setUniform1i("uRetroSteps", s.retroColorSteps);
    composite_.setUniform1f("uRetroDither", s.retroDitherStrength);
    composite_.setUniform1f("uRetroScanlines", s.retroScanlines);
    composite_.setUniform1f("uRetroCurvature", s.retroCurvature);
    composite_.setUniform1f("uRetroJitter", s.retroJitter);
    composite_.setUniform1f("uHyper", s.hyperspace);
    composite_.setUniform1f("uHyperTwist", s.hyperspaceTwist);
    composite_.setUniform1f("uHyperDensity", s.hyperspaceDensity);
    composite_.setUniform1f("uHyperNoise", s.hyperspaceNoise);
    composite_.setUniform1f("uHyperIntensity", s.hyperspaceIntensity);
    drawFullscreen();
  } else {
    frameGraph_.execute();
  }

  // --- Auto exposure update (read back 1x1 log-luminance) -------------------
  if (s.autoExposureEnabled && lumFinal.valid() && !lumFinal.isBackbuffer()) {
    const unsigned int lt = frameGraph_.glTexture(lumFinal);
    if (lt) {
      float rg[2] = {0.0f, 1.0f};
      gl::BindTexture(GL_TEXTURE_2D, lt);
      glGetTexImage(GL_TEXTURE_2D, 0, GL_RG, GL_FLOAT, rg);
      gl::BindTexture(GL_TEXTURE_2D, 0);

      const float wAvg = rg[1];
      const float avgLogLum = (wAvg > 1e-6f) ? (rg[0] / wAvg) : rg[0];
      if (std::isfinite(avgLogLum)) {
        const float avgLum = std::exp(avgLogLum);
        AutoExposureConfig cfg{};
        cfg.key = s.autoExposureKey;
        cfg.minExposure = s.autoExposureMin;
        cfg.maxExposure = s.autoExposureMax;
        cfg.speedUp = s.autoExposureSpeedUp;
        cfg.speedDown = s.autoExposureSpeedDown;
        stepAutoExposure(autoExposure_, avgLum, avgLogLum, dt, cfg);
      }
    }
  }

  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);
}

} // namespace stellar::render
