#include "stellar/render/ShaderToy.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cstddef>
#include <string>

namespace stellar::render {

namespace {

static const char* kFullscreenTriVS = R"GLSL(
#version 330 core
out vec2 vUv;

const vec2 kPos[3] = vec2[](
  vec2(-1.0, -1.0),
  vec2( 3.0, -1.0),
  vec2(-1.0,  3.0)
);

void main() {
  gl_Position = vec4(kPos[gl_VertexID], 0.0, 1.0);
  vUv = 0.5 * (gl_Position.xy + 1.0);
}
)GLSL";

// Clamp user snippet length to a sane maximum to avoid pathological allocations
// if someone pastes megabytes into the in-game editor.
static std::string_view clampUserCode(std::string_view s) {
  constexpr std::size_t kMax = 256u * 1024u;
  if (s.size() <= kMax) return s;
  return s.substr(0, kMax);
}

} // namespace

std::string buildShaderToyFragmentSource(std::string_view userCode, std::string_view extraHeader) {
  userCode = clampUserCode(userCode);

  std::string out;
  out.reserve(userCode.size() + 32 * 1024);

  out += "#version 330 core\n";
  out += "in vec2 vUv;\n";
  out += "out vec4 oColor;\n\n";

  // --- ShaderToy-ish uniforms ---
  out += "uniform vec2  iResolution;\n"; // (w,h)
  out += "uniform float iTime;\n";
  out += "uniform float iTimeDelta;\n";
  out += "uniform int   iFrame;\n";
  out += "uniform vec4  iMouse;\n";
  out += "uniform int   iSeed;\n";
  out += "uniform int   iPass;\n\n";

  // Optional input channels for multi-pass sketches.
  out += "uniform sampler2D iChannel0;\n";
  out += "uniform sampler2D iChannel1;\n";
  out += "uniform sampler2D iChannel2;\n";
  out += "uniform sampler2D iChannel3;\n";
  out += "uniform vec3  iChannelResolution[4];\n\n";

  // Camera uniforms for raymarch sketches.
  out += "uniform vec3  iCamPos;\n";
  out += "uniform vec3  iCamRight;\n";
  out += "uniform vec3  iCamUp;\n";
  out += "uniform vec3  iCamForward;\n";
  out += "uniform float iTanHalfFovY;\n";
  out += "uniform float iAspect;\n\n";

  // --- Tiny GLSL utility library (procedural building blocks) ---
  out += R"GLSL(

// ---- Math helpers ----
#define PI 3.14159265359

float saturate(float x) { return clamp(x, 0.0, 1.0); }
vec2  saturate(vec2 x)  { return clamp(x, 0.0, 1.0); }
vec3  saturate(vec3 x)  { return clamp(x, 0.0, 1.0); }

mat2 rot2(float a) {
  float s = sin(a), c = cos(a);
  return mat2(c, -s, s, c);
}

// Cosine palette.
vec3 palette(float t, vec3 a, vec3 b, vec3 c, vec3 d) {
  return a + b * cos(2.0 * PI * (c * t + d));
}

// ---- Integer hash (stable, no sin()) ----
uint hash1(uint x) {
  x ^= x >> 16u;
  x *= 0x7feb352du;
  x ^= x >> 15u;
  x *= 0x846ca68bu;
  x ^= x >> 16u;
  return x;
}

uint hash2(uvec2 v) {
  return hash1(v.x ^ hash1(v.y));
}

uint hash3(uvec3 v) {
  return hash1(v.x ^ hash1(v.y ^ hash1(v.z)));
}

float u01(uint x) {
  return float(x) * (1.0 / 4294967295.0);
}

float hash21(ivec2 p) {
  uvec2 u = uvec2(uint(p.x), uint(p.y));
  uint s = uint(iSeed);
  uint h = hash2(u + uvec2(s * 1664525u, s * 1013904223u));
  return u01(h);
}

vec2 hash22(ivec2 p) {
  float a = hash21(p);
  float b = hash21(p + ivec2(17, 113));
  return vec2(a, b);
}

float hash31(ivec3 p) {
  uvec3 u = uvec3(uint(p.x), uint(p.y), uint(p.z));
  uint s = uint(iSeed);
  uint h = hash3(u + uvec3(s * 1103515245u));
  return u01(h);
}

vec3 hash33(ivec3 p) {
  float a = hash31(p);
  float b = hash31(p + ivec3(19, 73, 7));
  float c = hash31(p + ivec3(41, 29, 97));
  return vec3(a, b, c);
}

// Quintic curve for lattice interpolation (C2 continuous).
float fade(float t) {
  return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

// ---- Value noise ----
float noise2(vec2 x) {
  ivec2 i = ivec2(floor(x));
  vec2  f = fract(x);
  vec2  u = vec2(fade(f.x), fade(f.y));

  float a = hash21(i);
  float b = hash21(i + ivec2(1, 0));
  float c = hash21(i + ivec2(0, 1));
  float d = hash21(i + ivec2(1, 1));

  return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

float noise3(vec3 x) {
  ivec3 i = ivec3(floor(x));
  vec3  f = fract(x);
  vec3  u = vec3(fade(f.x), fade(f.y), fade(f.z));

  float n000 = hash31(i + ivec3(0,0,0));
  float n100 = hash31(i + ivec3(1,0,0));
  float n010 = hash31(i + ivec3(0,1,0));
  float n110 = hash31(i + ivec3(1,1,0));
  float n001 = hash31(i + ivec3(0,0,1));
  float n101 = hash31(i + ivec3(1,0,1));
  float n011 = hash31(i + ivec3(0,1,1));
  float n111 = hash31(i + ivec3(1,1,1));

  float nx00 = mix(n000, n100, u.x);
  float nx10 = mix(n010, n110, u.x);
  float nx01 = mix(n001, n101, u.x);
  float nx11 = mix(n011, n111, u.x);
  float nxy0 = mix(nx00, nx10, u.y);
  float nxy1 = mix(nx01, nx11, u.y);
  return mix(nxy0, nxy1, u.z);
}

float fbm2(vec2 p) {
  float v = 0.0;
  float a = 0.5;
  mat2 m = mat2(1.6, 1.2, -1.2, 1.6);
  for (int i = 0; i < 6; ++i) {
    v += a * noise2(p);
    p = m * p;
    a *= 0.5;
  }
  return v;
}

float fbm3(vec3 p) {
  float v = 0.0;
  float a = 0.5;
  mat3 m = mat3( 1.6,  1.2,  0.0,
                -1.2,  1.6,  0.0,
                 0.0,  0.0,  1.35);
  for (int i = 0; i < 6; ++i) {
    v += a * noise3(p);
    p = m * p;
    a *= 0.5;
  }
  return v;
}

// ---- Dynamic LOD / band-limiting helpers ----
//
// These helpers use screen-space derivatives (dFdx/dFdy/fwidth) to estimate
// how quickly a procedural coordinate varies across a pixel. This allows you
// to *band-limit* procedural detail (especially FBM octaves) automatically
// when it becomes sub-pixel during zoom-out, reducing shimmer/aliasing.
//
// featurePx(): approximate feature size in pixels (bigger => smoother).
float featurePx(vec2 p) {
  vec2 dx = dFdx(p);
  vec2 dy = dFdy(p);
  float fw = max(length(dx), length(dy));
  return 1.0 / max(1.0e-6, fw);
}

float featurePx(vec3 p) {
  vec3 dx = dFdx(p);
  vec3 dy = dFdy(p);
  float fw = max(length(dx), length(dy));
  return 1.0 / max(1.0e-6, fw);
}

// Convert a feature size (px) into a [0..1] LOD weight.
// - pxLo: mostly filtered out below this size
// - pxHi: fully visible above this size
float lodWeightFromPx(float featPx, float pxLo, float pxHi) {
  float lo = min(pxLo, pxHi);
  float hi = max(pxLo, pxHi);
  hi = max(hi, lo + 1.0e-3);
  return smoothstep(lo, hi, featPx);
}

// Filtered FBM: automatically fades out high-frequency octaves as they become
// sub-pixel.
//
// energyComp in [0..1]:
//   0 => allow energy to drop (smoother / less contrast at distance)
//   1 => renormalize by the remaining octave weights (crisper / more constant)
float fbm2_lod(vec2 p, float pxLo, float pxHi, float energyComp) {
  float v = 0.0;
  float a = 0.5;
  float wsum = 0.0;
  mat2 m = mat2(1.6, 1.2, -1.2, 1.6);

  for (int i = 0; i < 6; ++i) {
    float w = lodWeightFromPx(featurePx(p), pxLo, pxHi);
    v += a * noise2(p) * w;
    wsum += a * w;
    p = m * p;
    a *= 0.5;
  }

  float vNorm = (wsum > 1.0e-6) ? (v / wsum) : v;
  return mix(v, vNorm, saturate(energyComp));
}

float fbm3_lod(vec3 p, float pxLo, float pxHi, float energyComp) {
  float v = 0.0;
  float a = 0.5;
  float wsum = 0.0;
  mat3 m = mat3( 1.6,  1.2,  0.0,
                -1.2,  1.6,  0.0,
                 0.0,  0.0,  1.35);

  for (int i = 0; i < 6; ++i) {
    float w = lodWeightFromPx(featurePx(p), pxLo, pxHi);
    v += a * noise3(p) * w;
    wsum += a * w;
    p = m * p;
    a *= 0.5;
  }

  float vNorm = (wsum > 1.0e-6) ? (v / wsum) : v;
  return mix(v, vNorm, saturate(energyComp));
}

// Worley / Voronoi distance (distance to nearest feature point).
float worley2(vec2 p) {
  ivec2 ip = ivec2(floor(p));
  float d2 = 1e9;
  for (int y = -1; y <= 1; ++y) {
    for (int x = -1; x <= 1; ++x) {
      ivec2 cell = ip + ivec2(x, y);
      vec2  rnd = hash22(cell);
      vec2  fp = vec2(cell) + rnd;
      vec2  r = fp - p;
      d2 = min(d2, dot(r, r));
    }
  }
  return sqrt(d2);
}

// Domain warp (2D) using two decorrelated fbm fields.
vec2 warp2(vec2 p, float strength) {
  vec2 q;
  q.x = fbm2(p + vec2(17.1, 3.2));
  q.y = fbm2(p + vec2( 5.7, 9.2));
  return p + strength * (q * 2.0 - 1.0);
}

// ---- SDF helpers (for 3D sketches) ----
float sdSphere(vec3 p, float r) { return length(p) - r; }

float sdBox(vec3 p, vec3 b) {
  vec3 d = abs(p) - b;
  return min(max(d.x, max(d.y, d.z)), 0.0) + length(max(d, 0.0));
}

// Camera ray reconstruction (uv in [0,1]).
vec3 rayDirFromUv(vec2 uv) {
  vec2 ndc = uv * 2.0 - 1.0;
  float x = ndc.x * iAspect * iTanHalfFovY;
  float y = ndc.y * iTanHalfFovY;
  return normalize(iCamForward + x * iCamRight + y * iCamUp);
}

// Cheap tonemap for HDR-ish sketches.
vec3 tonemapSimple(vec3 x) {
  return x / (1.0 + x);
}

)GLSL";

  // --- Optional user header (extra uniforms, macros, etc.) ---
  // Insert this BEFORE `#line 1` so snippet error logs keep correct line numbers.
  if (!extraHeader.empty()) {
    out.append(extraHeader.data(), extraHeader.size());
    if (!out.empty() && out.back() != '\n') out += "\n";
    out += "\n";
  }

  // --- User snippet ---
  out += "#line 1\n";
  out.append(userCode.data(), userCode.size());
  out += "\n\n";

  // --- Entry point ---
  out += R"GLSL(

// Expected snippet signature:
//   vec4 shaderMain(vec2 uv);
// where uv is in [0,1].

void main() {
  oColor = shaderMain(vUv);
}
)GLSL";

  return out;
}

ShaderToy::~ShaderToy() {
  if (vao_ != 0) {
    gl::DeleteVertexArrays(1, &vao_);
    vao_ = 0;
  }
}

bool ShaderToy::init(std::string* outError) {
  if (vao_ == 0) {
    gl::GenVertexArrays(1, &vao_);
  }

  // Default texture used for unconnected iChannel*s.
  if (blackTex_.handle() == 0) {
    const unsigned char px[4] = {0, 0, 0, 255};
    // nearestFilter=true avoids sampling garbage for feedback buffers; clamp prevents edge bleed.
    blackTex_.createRGBA(1, 1, px, /*generateMips=*/false, /*nearestFilter=*/true, /*clampToEdge=*/true);
  }

  // Build a tiny default shader so draw() can be called before build().
  const char* kDefaultSnippet = R"GLSL(
// Default snippet
vec4 shaderMain(vec2 uv) {
  vec2 p = uv * 2.0 - 1.0;
  p.x *= iAspect;
  float g = 0.5 + 0.5 * sin(iTime + 6.0 * length(p));
  vec3 col = palette(g, vec3(0.35), vec3(0.55), vec3(1.0), vec3(0.0, 0.33, 0.67));
  return vec4(col, 1.0);
}
)GLSL";

  lastVertex_ = kFullscreenTriVS;
  lastFragment_ = buildShaderToyFragmentSource(kDefaultSnippet, {});
  return shader_.build(lastVertex_, lastFragment_, outError);
}

bool ShaderToy::build(std::string_view userCode, std::string* outError) {
  return buildWithHeader(userCode, {}, outError);
}

bool ShaderToy::buildWithHeader(std::string_view userCode, std::string_view extraHeader, std::string* outError) {
  if (vao_ == 0) {
    // init() not called yet; do a minimal init.
    if (!init(outError)) return false;
  }

  lastVertex_ = kFullscreenTriVS;
  lastFragment_ = buildShaderToyFragmentSource(userCode, extraHeader);
  return shader_.build(lastVertex_, lastFragment_, outError);
}

void ShaderToy::draw(const ShaderToyInputs& in) const {
  if (shader_.handle() == 0 || vao_ == 0) return;

  const int w = std::max(1, in.width);
  const int h = std::max(1, in.height);

  ShaderToyInputs u = in;
  u.aspect = (h > 0) ? ((float)w / (float)h) : 1.0f;

  shader_.bind();
  shader_.setUniform2f("iResolution", (float)w, (float)h);
  shader_.setUniform1f("iTime", u.timeSec);
  shader_.setUniform1f("iTimeDelta", u.timeDeltaSec);
  shader_.setUniform1i("iFrame", u.frame);
  shader_.setUniform4f("iMouse", u.mouse[0], u.mouse[1], u.mouse[2], u.mouse[3]);
  shader_.setUniform1i("iSeed", u.seed);
  shader_.setUniform1i("iPass", u.passIndex);

  // Input channels (iChannel0..3). If a channel is not provided, bind a 1x1 black texture.
  shader_.setUniform1i("iChannel0", 0);
  shader_.setUniform1i("iChannel1", 1);
  shader_.setUniform1i("iChannel2", 2);
  shader_.setUniform1i("iChannel3", 3);

  const Texture2D* ch0 = u.channels[0] ? u.channels[0] : &blackTex_;
  const Texture2D* ch1 = u.channels[1] ? u.channels[1] : &blackTex_;
  const Texture2D* ch2 = u.channels[2] ? u.channels[2] : &blackTex_;
  const Texture2D* ch3 = u.channels[3] ? u.channels[3] : &blackTex_;

  shader_.setUniform3f("iChannelResolution[0]", (float)ch0->width(), (float)ch0->height(), 1.0f);
  shader_.setUniform3f("iChannelResolution[1]", (float)ch1->width(), (float)ch1->height(), 1.0f);
  shader_.setUniform3f("iChannelResolution[2]", (float)ch2->width(), (float)ch2->height(), 1.0f);
  shader_.setUniform3f("iChannelResolution[3]", (float)ch3->width(), (float)ch3->height(), 1.0f);

  ch0->bind(0);
  ch1->bind(1);
  ch2->bind(2);
  ch3->bind(3);

  shader_.setUniform3f("iCamPos", u.camPos[0], u.camPos[1], u.camPos[2]);
  shader_.setUniform3f("iCamRight", u.camRight[0], u.camRight[1], u.camRight[2]);
  shader_.setUniform3f("iCamUp", u.camUp[0], u.camUp[1], u.camUp[2]);
  shader_.setUniform3f("iCamForward", u.camForward[0], u.camForward[1], u.camForward[2]);
  shader_.setUniform1f("iTanHalfFovY", u.tanHalfFovY);
  shader_.setUniform1f("iAspect", u.aspect);

  // User parameters (if any).
  if (u.userParams) {
    u.userParams->applyToShader(shader_);
  }

  gl::BindVertexArray(vao_);
  glDrawArrays(GL_TRIANGLES, 0, 3);
}

} // namespace stellar::render
