#include "stellar/render/StarCoronaRenderer.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace stellar::render {

StarCoronaRenderer::~StarCoronaRenderer() {
  if (instanceVbo_) {
    gl::DeleteBuffers(1, &instanceVbo_);
    instanceVbo_ = 0;
  }
}

static const char* kVS = R"GLSL(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;

layout(location=3) in vec3 iPos;
layout(location=4) in vec3 iScale;
layout(location=5) in vec4 iQuat; // (x,y,z,w)
layout(location=6) in vec3 iColor;

uniform mat4 uView;
uniform mat4 uProj;

out vec3 vColor;
out vec3 vNrm;
out vec3 vWorldPos;
out vec3 vInstPos;

vec3 quatRotate(vec4 q, vec3 v) {
  // q assumed normalized
  vec3 t = 2.0 * cross(q.xyz, v);
  return v + q.w * t + cross(q.xyz, t);
}

void main() {
  vec3 local = aPos * iScale;
  vec3 pos = quatRotate(iQuat, local) + iPos;
  gl_Position = uProj * uView * vec4(pos, 1.0);
  vColor = iColor;
  vNrm = quatRotate(iQuat, aNrm);
  vWorldPos = pos;
  vInstPos = iPos;
}
)GLSL";

static const char* kFS = R"GLSL(
#version 330 core
in vec3 vColor;
in vec3 vNrm;
in vec3 vWorldPos;
in vec3 vInstPos;

uniform vec3 uCamPos;
uniform float uTime;

uniform int   uSeed;
uniform float uIntensity;
uniform float uRimPower;

uniform float uNoiseFreq;
uniform float uNoiseStrength;
uniform float uFlowSpeed;

uniform float uStreamerStrength;
uniform float uStreamerSharpness;

uniform float uPromStrength;
uniform int   uPromCount;
uniform float uPromSharpness;

uniform float uPulseStrength;

out vec4 FragColor;

float clamp01(float x) { return clamp(x, 0.0, 1.0); }

// Quintic fade
float fade(float t) {
  return t*t*t*(t*(t*6.0 - 15.0) + 10.0);
}

uint pcg_hash(uint v) {
  uint state = v * 747796405u + 2891336453u;
  uint word  = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
  return (word >> 22u) ^ word;
}

// Seeded float in [0,1].
float rand01(uint s) {
  return float(pcg_hash(s)) / 4294967295.0;
}

// Deterministic random unit vector from a seed.
vec3 randUnit(uint s) {
  float z = rand01(s + 11u) * 2.0 - 1.0;
  float a = rand01(s + 29u) * 6.28318530718;
  float r = sqrt(max(0.0, 1.0 - z*z));
  return vec3(r * cos(a), z, r * sin(a));
}

float hash01(uvec3 p, uint seed) {
  uint h = pcg_hash(p.x + pcg_hash(p.y + pcg_hash(p.z + seed)));
  return float(h) / 4294967295.0;
}

float smoothNoise3D(vec3 x, uint seed) {
  vec3 p0 = floor(x);
  vec3 f  = fract(x);

  uvec3 i0 = uvec3(p0);

  float n000 = hash01(i0 + uvec3(0,0,0), seed);
  float n100 = hash01(i0 + uvec3(1,0,0), seed);
  float n010 = hash01(i0 + uvec3(0,1,0), seed);
  float n110 = hash01(i0 + uvec3(1,1,0), seed);
  float n001 = hash01(i0 + uvec3(0,0,1), seed);
  float n101 = hash01(i0 + uvec3(1,0,1), seed);
  float n011 = hash01(i0 + uvec3(0,1,1), seed);
  float n111 = hash01(i0 + uvec3(1,1,1), seed);

  float u = fade(f.x);
  float v = fade(f.y);
  float w = fade(f.z);

  float x00 = mix(n000, n100, u);
  float x10 = mix(n010, n110, u);
  float x01 = mix(n001, n101, u);
  float x11 = mix(n011, n111, u);
  float y0v = mix(x00, x10, v);
  float y1v = mix(x01, x11, v);
  return mix(y0v, y1v, w);
}

float fbm3D(vec3 p, uint seed) {
  float amp = 0.5;
  float freq = 1.0;
  float sum = 0.0;
  for (int i = 0; i < 6; ++i) {
    sum += amp * smoothNoise3D(p * freq, seed + uint(i) * 1013u);
    freq *= 2.07;
    amp *= 0.52;
  }
  return sum;
}

// Helper: build an orthonormal basis given an axis.
void basisFromAxis(vec3 axis, out vec3 bx, out vec3 by) {
  vec3 up = (abs(axis.y) < 0.95) ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
  bx = normalize(cross(up, axis));
  by = normalize(cross(axis, bx));
}

void main() {
  vec3 n = normalize(vNrm);
  vec3 v = normalize(uCamPos - vWorldPos);
  float ndotv = clamp(dot(n, v), 0.0, 1.0);

  // Rim/limb factor. Corona should mostly appear at grazing angles.
  float rim = pow(1.0 - ndotv, max(0.25, uRimPower));

  // Direction from star center (seam-free domain).
  vec3 dir = normalize(vWorldPos - vInstPos);

  uint s = uint(max(uSeed, 1));
  vec3 axis = randUnit(s ^ 0xA57Eu);
  vec3 bx, by;
  basisFromAxis(axis, bx, by);
  float lat = dot(dir, axis);
  float az  = atan(dot(dir, by), dot(dir, bx));

  float phase = rand01(s ^ 0xB16Du) * 6.28318530718;

  // Animated flow: advect the noise domain along a deterministic direction.
  vec3 flowDir = randUnit(s ^ 0xC0DEu);
  float t = uTime * max(uFlowSpeed, 0.0);

  float freq = max(0.001, uNoiseFreq);
  float n0 = fbm3D(dir * freq + flowDir * t, s ^ 0x51A7u);
  float n1 = fbm3D(dir * (freq * 1.9) + flowDir * (t * 1.35), s ^ 0x70ADu);
  float n = clamp01(0.65 * n0 + 0.35 * n1);

  // Edge fray/filaments.
  float edge = 1.0 + (n - 0.5) * 2.0 * clamp(uNoiseStrength, 0.0, 3.0);
  edge = max(edge, 0.0);

  // Streamers: azimuthal rays modulated by an equatorial band.
  float eq = pow(clamp01(1.0 - abs(lat)), 1.25);
  float rays = 0.5 + 0.5 * cos(float(max(uPromCount, 1)) * az + phase + (n - 0.5) * 2.2);
  rays = pow(clamp01(rays), max(0.5, uStreamerSharpness));
  float streamer = rays * eq;

  // Prominence lobes: fewer, fatter arcs.
  int pc = clamp(uPromCount, 1, 24);
  float lobes = 0.5 + 0.5 * sin(float(pc) * 0.65 * az + phase * 1.3 + (n1 - 0.5) * 3.0);
  lobes = pow(clamp01(lobes), max(0.5, uPromSharpness));
  float prom = lobes * eq;

  // Subtle pulsation so the corona feels alive.
  float pulse = 1.0 + clamp(uPulseStrength, 0.0, 0.5) * (0.5 + 0.5 * sin(uTime * 0.65 + phase));

  float base = rim * edge;
  float streamAdd = rim * clamp(uStreamerStrength, 0.0, 4.0) * streamer;
  float promAdd   = rim * clamp(uPromStrength, 0.0, 4.0) * prom;

  float a = (base + streamAdd + promAdd) * pulse;
  vec3 col = vColor * (uIntensity * a);

  // Additive blend uses alpha as a multiplier; keep alpha=1 for simplicity.
  FragColor = vec4(col, 1.0);
}
)GLSL";

bool StarCoronaRenderer::init(std::string* outError) {
  if (!shader_.build(kVS, kFS, outError)) return false;

  gl::GenBuffers(1, &instanceVbo_);

  // Default identity matrices.
  for (int i = 0; i < 16; ++i) {
    view_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
    proj_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
  }

  return true;
}

void StarCoronaRenderer::setViewProj(const float* view, const float* proj) {
  std::memcpy(view_, view, sizeof(float) * 16);
  std::memcpy(proj_, proj, sizeof(float) * 16);
}

void StarCoronaRenderer::drawInstances(const std::vector<InstanceData>& instances) {
  if (!mesh_ || instances.empty()) return;

  // Additive halo.
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE);

  const int seed = std::max(settings_.seed, 1);

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);
  shader_.setUniform3f("uCamPos", camPos_[0], camPos_[1], camPos_[2]);
  shader_.setUniform1f("uTime", timeSeconds_);

  shader_.setUniform1i("uSeed", seed);
  shader_.setUniform1f("uIntensity", settings_.intensity);
  shader_.setUniform1f("uRimPower", settings_.rimPower);
  shader_.setUniform1f("uNoiseFreq", settings_.noiseFrequency);
  shader_.setUniform1f("uNoiseStrength", settings_.noiseStrength);
  shader_.setUniform1f("uFlowSpeed", settings_.flowSpeed);
  shader_.setUniform1f("uStreamerStrength", settings_.streamerStrength);
  shader_.setUniform1f("uStreamerSharpness", settings_.streamerSharpness);
  shader_.setUniform1f("uPromStrength", settings_.prominenceStrength);
  shader_.setUniform1i("uPromCount", settings_.prominenceCount);
  shader_.setUniform1f("uPromSharpness", settings_.prominenceSharpness);
  shader_.setUniform1f("uPulseStrength", settings_.pulseStrength);

  mesh_->bind();

  gl::BindBuffer(GL_ARRAY_BUFFER, instanceVbo_);
  gl::BufferData(GL_ARRAY_BUFFER,
                 static_cast<GLsizeiptr>(instances.size() * sizeof(InstanceData)),
                 instances.data(),
                 GL_DYNAMIC_DRAW);

  // location 3: vec3 position
  gl::EnableVertexAttribArray(3);
  gl::VertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, px));
  gl::VertexAttribDivisor(3, 1);

  // location 4: vec3 scale
  gl::EnableVertexAttribArray(4);
  gl::VertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, sx));
  gl::VertexAttribDivisor(4, 1);

  // location 5: vec4 quaternion
  gl::EnableVertexAttribArray(5);
  gl::VertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, qx));
  gl::VertexAttribDivisor(5, 1);

  // location 6: vec3 color
  gl::EnableVertexAttribArray(6);
  gl::VertexAttribPointer(6, 3, GL_FLOAT, GL_FALSE, sizeof(InstanceData), (void*)offsetof(InstanceData, cr));
  gl::VertexAttribDivisor(6, 1);

  mesh_->drawInstanced(static_cast<std::uint32_t>(instances.size()));
}

} // namespace stellar::render
