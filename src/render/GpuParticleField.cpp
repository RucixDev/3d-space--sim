#include "stellar/render/GpuParticleField.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>

namespace stellar::render {

namespace {

static const char* kSimVS = R"GLSL(
#version 330 core
out vec2 vUv;

// Full-screen triangle (no VBO). We still output UV for debugging, but
// simulation uses gl_FragCoord for exact texel addressing.
void main() {
  vec2 verts[3] = vec2[3](vec2(-1.0, -1.0), vec2(3.0, -1.0), vec2(-1.0, 3.0));
  vec2 p = verts[gl_VertexID];
  gl_Position = vec4(p, 0.0, 1.0);
  vUv = 0.5 * p + vec2(0.5);
}
)GLSL";

// Init shader outputs either position or velocity depending on uMode.
// uMode: 0=pos, 1=vel
static const char* kInitFS = R"GLSL(
#version 330 core
in vec2 vUv;
out vec4 FragColor;

uniform int   uSeed;
uniform int   uMode;
uniform float uBounds;
uniform float uInitSpeed;

// PCG-style hash (cheap, decent quality)
uint pcg_hash(uint v) {
  uint state = v * 747796405u + 2891336453u;
  uint word  = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
  return (word >> 22u) ^ word;
}

float u01(uint v) {
  return float(pcg_hash(v)) / 4294967295.0;
}

vec3 rand3(uvec2 p, uint salt) {
  uint b = uint(uSeed);
  uint h0 = pcg_hash(p.x + pcg_hash(p.y + pcg_hash(b + salt)));
  uint h1 = pcg_hash(h0 + 0x9E3779B9u);
  uint h2 = pcg_hash(h1 + 0x9E3779B9u);
  return vec3(float(h0), float(h1), float(h2)) / 4294967295.0;
}

// Uniformly sample inside a sphere via rejection sampling (few iterations).
vec3 randInSphere(uvec2 p) {
  vec3 r = vec3(0.0);
  for (int i = 0; i < 6; ++i) {
    vec3 x = rand3(p, uint(i) * 1013u) * 2.0 - 1.0;
    if (dot(x,x) <= 1.0) {
      r = x;
      break;
    }
  }
  return r;
}

void main() {
  uvec2 ip = uvec2(ivec2(gl_FragCoord.xy));

  if (uMode == 0) {
    vec3 p = randInSphere(ip) * uBounds;
    FragColor = vec4(p, 1.0);
  } else {
    vec3 v = (rand3(ip, 777u) * 2.0 - 1.0) * uInitSpeed;
    FragColor = vec4(v, 1.0);
  }
}
)GLSL";

// Velocity update: v += dt * (noise - attract*pos) and apply drag.
static const char* kVelFS = R"GLSL(
#version 330 core
out vec4 FragColor;

uniform sampler2D uPos;
uniform sampler2D uVel;

uniform float uDt;
uniform float uTime;
uniform vec3  uShift;

uniform int   uSeed;
uniform float uBounds;
uniform float uDrag;
uniform float uNoiseFreq;
uniform float uNoiseStrength;
uniform float uAttract;
uniform float uMaxSpeed;

uint pcg_hash(uint v) {
  uint state = v * 747796405u + 2891336453u;
  uint word  = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
  return (word >> 22u) ^ word;
}

float fade(float t) {
  return t*t*t*(t*(t*6.0 - 15.0) + 10.0);
}

float hash01(uvec3 p, uint seed) {
  uint h = pcg_hash(p.x + pcg_hash(p.y + pcg_hash(p.z + seed)));
  return float(h) / 4294967295.0;
}

float noise3(vec3 x, uint seed) {
  vec3 i = floor(x);
  vec3 f = fract(x);
  vec3 u = vec3(fade(f.x), fade(f.y), fade(f.z));

  uvec3 i0 = uvec3(i);
  float n000 = hash01(i0 + uvec3(0,0,0), seed);
  float n100 = hash01(i0 + uvec3(1,0,0), seed);
  float n010 = hash01(i0 + uvec3(0,1,0), seed);
  float n110 = hash01(i0 + uvec3(1,1,0), seed);
  float n001 = hash01(i0 + uvec3(0,0,1), seed);
  float n101 = hash01(i0 + uvec3(1,0,1), seed);
  float n011 = hash01(i0 + uvec3(0,1,1), seed);
  float n111 = hash01(i0 + uvec3(1,1,1), seed);

  float x00 = mix(n000, n100, u.x);
  float x10 = mix(n010, n110, u.x);
  float x01 = mix(n001, n101, u.x);
  float x11 = mix(n011, n111, u.x);

  float y0 = mix(x00, x10, u.y);
  float y1 = mix(x01, x11, u.y);
  return mix(y0, y1, u.z);
}

float fbm3(vec3 p, uint seed) {
  float amp = 0.5;
  float freq = 1.0;
  float sum = 0.0;
  for (int i = 0; i < 4; ++i) {
    sum += amp * noise3(p * freq, seed + uint(i) * 1013u);
    freq *= 2.0;
    amp *= 0.5;
  }
  return sum;
}

vec3 noiseVec(vec3 p, uint seed) {
  float nx = fbm3(p + vec3(0.0, 0.0, 0.0), seed);
  float ny = fbm3(p + vec3(19.1, 7.7, 3.3), seed + 17u);
  float nz = fbm3(p + vec3(5.2, 13.9, 23.7), seed + 31u);
  return (vec3(nx, ny, nz) * 2.0 - 1.0);
}

void main() {
  ivec2 ip = ivec2(gl_FragCoord.xy);
  vec3 pos = texelFetch(uPos, ip, 0).xyz - uShift;
  vec3 vel = texelFetch(uVel, ip, 0).xyz;

  uint seed = uint(uSeed);

  // Coherent noise in a moving field.
  vec3 p = pos * uNoiseFreq + vec3(0.07, 0.03, 0.05) * uTime;
  vec3 n = noiseVec(p, seed);

  // Mild "eddy" bias so it doesn't look like random jitter.
  vec3 swirl = vec3(n.y - n.z, n.z - n.x, n.x - n.y);
  vec3 acc = swirl * uNoiseStrength;

  // Optional pull toward local origin.
  acc += (-pos) * uAttract;

  // Semi-implicit drag.
  float drag = max(0.0, uDrag);
  float damp = exp(-drag * max(0.0, uDt));
  vel = vel * damp + acc * uDt;

  // Cap speed.
  float sp = length(vel);
  if (sp > uMaxSpeed && sp > 1e-6) {
    vel = vel * (uMaxSpeed / sp);
  }

  FragColor = vec4(vel, 1.0);
}
)GLSL";

// Position update: apply shift, integrate with new vel, and wrap into bounds.
static const char* kPosFS = R"GLSL(
#version 330 core
out vec4 FragColor;

uniform sampler2D uPos;
uniform sampler2D uVel;
uniform float uDt;
uniform vec3  uShift;
uniform float uBounds;

vec3 wrapBox(vec3 p, float b) {
  // Wrap each component to [-b, b] (toroidal).
  vec3 d = vec3(2.0 * b);
  p = mod(p + vec3(b), d) - vec3(b);
  return p;
}

void main() {
  ivec2 ip = ivec2(gl_FragCoord.xy);
  vec3 pos = texelFetch(uPos, ip, 0).xyz - uShift;
  vec3 vel = texelFetch(uVel, ip, 0).xyz;
  pos += vel * uDt;
  pos = wrapBox(pos, max(0.001, uBounds));
  FragColor = vec4(pos, 1.0);
}
)GLSL";

// Draw: sample particle positions with gl_VertexID.
static const char* kDrawVS = R"GLSL(
#version 330 core

uniform mat4 uView;
uniform mat4 uProj;

uniform sampler2D uPos;
uniform sampler2D uVel;

uniform int   uDim;
uniform int   uSeed;
uniform float uBounds;
uniform float uPointSize;
uniform float uAlpha;
uniform float uIntensity;

out vec3  vColor;
out float vAlpha;

uint pcg_hash(uint v) {
  uint state = v * 747796405u + 2891336453u;
  uint word  = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
  return (word >> 22u) ^ word;
}

float u01(uint v) {
  return float(pcg_hash(v)) / 4294967295.0;
}

vec3 hsv2rgb(vec3 c) {
  vec4 K = vec4(1.0, 2.0/3.0, 1.0/3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main() {
  int id = gl_VertexID;
  int x = id % uDim;
  int y = id / uDim;
  ivec2 ip = ivec2(x, y);

  vec3 pos = texelFetch(uPos, ip, 0).xyz;
  vec3 vel = texelFetch(uVel, ip, 0).xyz;

  // Fade particles near the box boundary to avoid visible tiling.
  float r = length(pos);
  float edge = smoothstep(uBounds * 0.85, uBounds, r);
  float fade = 1.0 - edge;

  // Color variation by stable hash.
  uint h = pcg_hash(uint(id) ^ uint(uSeed));
  float hue = u01(h);
  float sat = 0.25 + 0.35 * u01(h + 0xA511E9B3u);
  float val = 0.75 + 0.35 * u01(h + 0x63D83595u);
  vec3 col = hsv2rgb(vec3(hue, sat, val));

  // Speed -> brighter (helps sell motion).
  float sp = length(vel);
  float spBoost = 0.65 + 0.35 * clamp(sp / max(0.001, 8.0), 0.0, 1.0);

  vColor = col * (uIntensity * spBoost);
  vAlpha = uAlpha * fade;

  gl_Position = uProj * uView * vec4(pos, 1.0);

  float sizeRnd = 0.65 + 0.7 * u01(h + 0xB5297A4Du);
  gl_PointSize = uPointSize * sizeRnd;
}
)GLSL";

static const char* kDrawFS = R"GLSL(
#version 330 core
in vec3  vColor;
in float vAlpha;
out vec4 FragColor;

uniform int uUseSprite;
uniform sampler2D uSprite;

void main() {
  // Soft circular mask.
  vec2 p = gl_PointCoord * 2.0 - 1.0;
  float d = dot(p,p);
  float circle = smoothstep(1.0, 0.7, d);

  if (uUseSprite == 0) {
    FragColor = vec4(vColor, vAlpha * circle);
  } else {
    vec4 t = texture(uSprite, gl_PointCoord);
    float a = t.a * vAlpha * circle;
    FragColor = vec4(t.rgb * vColor, a);
  }
}
)GLSL";

static void setIdentity(float m[16]) {
  for (int i = 0; i < 16; ++i) m[i] = (i % 5 == 0) ? 1.0f : 0.0f;
}

static void appendError(std::string* outError, const std::string& msg) {
  if (!outError) return;
  if (!outError->empty()) *outError += "\n";
  *outError += msg;
}

} // namespace

GpuParticleField::~GpuParticleField() { destroyGL(); }

void GpuParticleField::destroyGL() {
  if (posTex_[0]) gl::DeleteTextures(2, posTex_);
  if (velTex_[0]) gl::DeleteTextures(2, velTex_);
  posTex_[0] = posTex_[1] = 0;
  velTex_[0] = velTex_[1] = 0;

  if (fbo_) {
    gl::DeleteFramebuffers(1, &fbo_);
    fbo_ = 0;
  }
  if (vao_) {
    gl::DeleteVertexArrays(1, &vao_);
    vao_ = 0;
  }

  dim_ = 0;
  ping_ = 0;
  inited_ = false;
  needsInit_ = true;
}

bool GpuParticleField::init(std::string* outError) {
  destroyGL();

  std::string err;
  if (!initShader_.build(kSimVS, kInitFS, &err)) {
    appendError(outError, "GpuParticleField: init shader failed: " + err);
    return false;
  }
  if (!velShader_.build(kSimVS, kVelFS, &err)) {
    appendError(outError, "GpuParticleField: vel shader failed: " + err);
    return false;
  }
  if (!posShader_.build(kSimVS, kPosFS, &err)) {
    appendError(outError, "GpuParticleField: pos shader failed: " + err);
    return false;
  }
  if (!drawShader_.build(kDrawVS, kDrawFS, &err)) {
    appendError(outError, "GpuParticleField: draw shader failed: " + err);
    return false;
  }

  gl::GenVertexArrays(1, &vao_);
  gl::GenFramebuffers(1, &fbo_);

  setIdentity(view_);
  setIdentity(proj_);

  // Default sampler bindings.
  drawShader_.bind();
  drawShader_.setUniform1i("uPos", 0);
  drawShader_.setUniform1i("uVel", 1);
  drawShader_.setUniform1i("uSprite", 2);

  velShader_.bind();
  velShader_.setUniform1i("uPos", 0);
  velShader_.setUniform1i("uVel", 1);

  posShader_.bind();
  posShader_.setUniform1i("uPos", 0);
  posShader_.setUniform1i("uVel", 1);

  inited_ = true;
  return true;
}

void GpuParticleField::setViewProj(const float* view, const float* proj) {
  if (view) std::memcpy(view_, view, sizeof(float) * 16);
  if (proj) std::memcpy(proj_, proj, sizeof(float) * 16);
}

void GpuParticleField::reset() {
  needsInit_ = true;
}

static void allocFloatTex2D(unsigned int tex, int w, int h) {
  gl::BindTexture(GL_TEXTURE_2D, tex);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  gl::TexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  gl::TexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
}

bool GpuParticleField::ensureSize(int dim, std::string* outError) {
  dim = std::clamp(dim, 8, 2048);
  if (dim_ == dim && posTex_[0] && velTex_[0]) return true;

  // Reallocate.
  if (posTex_[0]) gl::DeleteTextures(2, posTex_);
  if (velTex_[0]) gl::DeleteTextures(2, velTex_);
  posTex_[0] = posTex_[1] = 0;
  velTex_[0] = velTex_[1] = 0;

  gl::GenTextures(2, posTex_);
  gl::GenTextures(2, velTex_);
  for (int i = 0; i < 2; ++i) {
    allocFloatTex2D(posTex_[i], dim, dim);
    allocFloatTex2D(velTex_[i], dim, dim);
  }

  // Verify FBO completeness with one attachment.
  gl::BindFramebuffer(GL_FRAMEBUFFER, fbo_);
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, posTex_[0], 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);
  const GLenum status = gl::CheckFramebufferStatus(GL_FRAMEBUFFER);
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
  gl::BindFramebuffer(GL_FRAMEBUFFER, 0);
  if (status != GL_FRAMEBUFFER_COMPLETE) {
    std::ostringstream ss;
    ss << "GpuParticleField: framebuffer incomplete (status=0x" << std::hex << (unsigned int)status << ")";
    appendError(outError, ss.str());
    return false;
  }

  dim_ = dim;
  ping_ = 0;
  needsInit_ = true;
  return true;
}

bool GpuParticleField::initState(const Settings& s, std::string* outError) {
  if (!ensureSize(s.dim, outError)) return false;

  // Save state we might clobber.
  GLint prevFbo = 0;
  GLint prevViewport[4] = {0,0,0,0};
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
  glGetIntegerv(GL_VIEWPORT, prevViewport);

  gl::BindFramebuffer(GL_FRAMEBUFFER, fbo_);
  glViewport(0, 0, dim_, dim_);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);

  gl::BindVertexArray(vao_);
  initShader_.bind();
  initShader_.setUniform1i("uSeed", s.seed);
  initShader_.setUniform1f("uBounds", std::max(0.001f, s.boundsU));
  initShader_.setUniform1f("uInitSpeed", std::max(0.0f, std::min(s.maxSpeed, 25.0f)) * 0.35f);

  // Position
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, posTex_[0], 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);
  initShader_.setUniform1i("uMode", 0);
  glDrawArrays(GL_TRIANGLES, 0, 3);

  // Velocity
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, velTex_[0], 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);
  initShader_.setUniform1i("uMode", 1);
  glDrawArrays(GL_TRIANGLES, 0, 3);

  // Detach.
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
  gl::BindVertexArray(0);

  // Restore state.
  gl::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
  glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);

  lastSeed_ = s.seed;
  needsInit_ = false;
  ping_ = 0;
  return true;
}

void GpuParticleField::update(double dtSec,
                              const float shiftDeltaU[3],
                              float timeSec,
                              const Settings& s) {
  if (!inited_) return;
  if (!s.enabled) return;

  std::string err;
  if (s.dim != dim_ || !posTex_[0] || !velTex_[0]) {
    if (!ensureSize(s.dim, &err)) {
      // If allocation fails, disable by marking not inited.
      inited_ = false;
      return;
    }
  }

  if (needsInit_ || s.seed != lastSeed_) {
    (void)initState(s, &err);
  }

  const float dt = (float)std::max(0.0, dtSec);
  if (dt <= 0.0f && (!shiftDeltaU || (shiftDeltaU[0] == 0.0f && shiftDeltaU[1] == 0.0f && shiftDeltaU[2] == 0.0f))) {
    return;
  }

  // Save state.
  GLint prevFbo = 0;
  GLint prevViewport[4] = {0,0,0,0};
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
  glGetIntegerv(GL_VIEWPORT, prevViewport);

  gl::BindFramebuffer(GL_FRAMEBUFFER, fbo_);
  glViewport(0, 0, dim_, dim_);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);

  const int r = readIndex();
  const int w = writeIndex();

  const float shift[3] = {
      shiftDeltaU ? shiftDeltaU[0] : 0.0f,
      shiftDeltaU ? shiftDeltaU[1] : 0.0f,
      shiftDeltaU ? shiftDeltaU[2] : 0.0f,
  };

  gl::BindVertexArray(vao_);

  // --- Velocity pass (writes vel[w]) ---
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, velTex_[w], 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);

  velShader_.bind();
  velShader_.setUniform1f("uDt", dt);
  velShader_.setUniform1f("uTime", timeSec);
  velShader_.setUniform3f("uShift", shift[0], shift[1], shift[2]);
  velShader_.setUniform1i("uSeed", s.seed);
  velShader_.setUniform1f("uBounds", std::max(0.001f, s.boundsU));
  velShader_.setUniform1f("uDrag", std::max(0.0f, s.drag));
  velShader_.setUniform1f("uNoiseFreq", std::max(0.00001f, s.noiseFrequency));
  velShader_.setUniform1f("uNoiseStrength", std::max(0.0f, s.noiseStrength));
  velShader_.setUniform1f("uAttract", std::max(0.0f, s.attractStrength));
  velShader_.setUniform1f("uMaxSpeed", std::max(0.01f, s.maxSpeed));

  gl::ActiveTexture(GL_TEXTURE0);
  gl::BindTexture(GL_TEXTURE_2D, posTex_[r]);
  gl::ActiveTexture(GL_TEXTURE1);
  gl::BindTexture(GL_TEXTURE_2D, velTex_[r]);
  glDrawArrays(GL_TRIANGLES, 0, 3);

  // --- Position pass (writes pos[w]) ---
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, posTex_[w], 0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);

  posShader_.bind();
  posShader_.setUniform1f("uDt", dt);
  posShader_.setUniform3f("uShift", shift[0], shift[1], shift[2]);
  posShader_.setUniform1f("uBounds", std::max(0.001f, s.boundsU));

  gl::ActiveTexture(GL_TEXTURE0);
  gl::BindTexture(GL_TEXTURE_2D, posTex_[r]);
  gl::ActiveTexture(GL_TEXTURE1);
  gl::BindTexture(GL_TEXTURE_2D, velTex_[w]);
  glDrawArrays(GL_TRIANGLES, 0, 3);

  // Detach.
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);
  gl::BindVertexArray(0);

  // Restore state.
  gl::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
  glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);

  ping_ ^= 1;
}

void GpuParticleField::draw(unsigned int spriteTexHandle, const Settings& s) {
  if (!inited_ || !s.enabled) return;
  if (!posTex_[readIndex()] || !velTex_[readIndex()]) return;

  const int r = readIndex();

  glEnable(GL_PROGRAM_POINT_SIZE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE); // additive for dust

  drawShader_.bind();
  drawShader_.setUniformMat4("uView", view_);
  drawShader_.setUniformMat4("uProj", proj_);
  drawShader_.setUniform1i("uDim", dim_);
  drawShader_.setUniform1i("uSeed", s.seed);
  drawShader_.setUniform1f("uBounds", std::max(0.001f, s.boundsU));
  drawShader_.setUniform1f("uPointSize", std::max(0.25f, s.pointSizePx));
  drawShader_.setUniform1f("uAlpha", std::clamp(s.alpha, 0.0f, 1.0f));
  drawShader_.setUniform1f("uIntensity", std::max(0.0f, s.intensity));

  const int useSprite = (s.textured && spriteTexHandle != 0) ? 1 : 0;
  drawShader_.setUniform1i("uUseSprite", useSprite);

  gl::ActiveTexture(GL_TEXTURE0);
  gl::BindTexture(GL_TEXTURE_2D, posTex_[r]);
  gl::ActiveTexture(GL_TEXTURE1);
  gl::BindTexture(GL_TEXTURE_2D, velTex_[r]);

  if (useSprite) {
    gl::ActiveTexture(GL_TEXTURE2);
    gl::BindTexture(GL_TEXTURE_2D, spriteTexHandle);
  }

  gl::BindVertexArray(vao_);
  const GLsizei count = (GLsizei)std::min<std::int64_t>((std::int64_t)dim_ * (std::int64_t)dim_, (std::int64_t)0x7fffffff);
  glDrawArrays(GL_POINTS, 0, count);
  gl::BindVertexArray(0);
}

} // namespace stellar::render
