#include "stellar/render/MeshRenderer.h"

#include "stellar/math/Quat.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>

namespace stellar::render {

// Some Windows OpenGL headers omit these tokens even though glMapBufferRange is
// available via extension loading.
#ifndef GL_MAP_WRITE_BIT
#define GL_MAP_WRITE_BIT 0x0002
#endif
#ifndef GL_MAP_INVALIDATE_RANGE_BIT
#define GL_MAP_INVALIDATE_RANGE_BIT 0x0004
#endif


MeshRenderer::~MeshRenderer() {
  if (instanceVbo_) {
    gl::DeleteBuffers(1, &instanceVbo_);
    instanceVbo_ = 0;
  }
}

static const char* kVS = R"GLSL(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNrm;
layout(location=2) in vec2 aUv;

layout(location=3) in vec3 iPos;
layout(location=4) in vec3 iScale;
layout(location=5) in vec4 iQuat; // (x,y,z,w)
layout(location=6) in vec3 iColor;

uniform mat4 uView;
uniform mat4 uProj;

out vec2 vUv;
out vec3 vColor;
out vec3 vNrm;
out vec3 vWorldPos;

vec3 quatRotate(vec4 q, vec3 v) {
  // q assumed normalized
  vec3 t = 2.0 * cross(q.xyz, v);
  return v + q.w * t + cross(q.xyz, t);
}

void main() {
  vec3 local = aPos * iScale;
  vec3 pos = quatRotate(iQuat, local) + iPos;
  gl_Position = uProj * uView * vec4(pos, 1.0);
  vUv = aUv;
  vColor = iColor;
  vNrm = quatRotate(iQuat, aNrm);
  vWorldPos = pos;
}
)GLSL";

static const char* kFS = R"GLSL(
#version 330 core
in vec2 vUv;
in vec3 vColor;
in vec3 vNrm;
in vec3 vWorldPos;

uniform sampler2D uTex;
uniform sampler2D uNormalTex;
uniform float uUseNormalMap;
uniform float uNormalStrength;

uniform sampler2D uEmissiveTex;
uniform float uUseEmissive;
uniform float uEmissiveStrength;
uniform float uEmissiveNightFadeStart;
uniform float uEmissiveNightFadeEnd;

uniform float uUnlit;
uniform vec3 uLightPos;
uniform vec3 uCameraPos;

uniform float uSpecularStrength;
uniform float uShininess;

uniform float uAlphaFromTex;
uniform float uAlphaMul;

out vec4 FragColor;

mat3 cotangentFrame(vec3 N, vec3 p, vec2 uv) {
  // Cotangent frame / derivative-based tangent space.
  // This avoids requiring per-vertex tangents for normal mapping.
  vec3 dp1 = dFdx(p);
  vec3 dp2 = dFdy(p);
  vec2 duv1 = dFdx(uv);
  vec2 duv2 = dFdy(uv);

  vec3 dp2perp = cross(dp2, N);
  vec3 dp1perp = cross(N, dp1);
  vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
  vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;

  float invMax = inversesqrt(max(dot(T, T), dot(B, B)));
  return mat3(T * invMax, B * invMax, N);
}

void main() {
  vec3 n = normalize(vNrm);

  if (uUseNormalMap > 0.5) {
    vec3 tn = texture(uNormalTex, vUv).xyz * 2.0 - 1.0;
    tn.xy *= max(uNormalStrength, 0.0);
    n = normalize(cotangentFrame(n, vWorldPos, vUv) * tn);
  }

  // Point-light position in world space (the star sits at the origin in the prototype).
  vec3 l = normalize(uLightPos - vWorldPos);
  float ndotl = dot(n, l);
  float diff = max(ndotl, 0.0);

  vec4 t = texture(uTex, vUv);
  vec3 tex = t.rgb;

  float lit = (0.35 + 0.65 * diff);
  float shade = mix(lit, 1.0, clamp(uUnlit, 0.0, 1.0));
  vec3 col = tex * vColor * shade;

  // Simple Blinn-Phong specular (disabled for unlit meshes).
  if (uUnlit < 0.5 && uSpecularStrength > 0.0) {
    vec3 v = normalize(uCameraPos - vWorldPos);
    vec3 h = normalize(l + v);
    float s = pow(max(dot(n, h), 0.0), max(uShininess, 1.0));
    col += vec3(s * uSpecularStrength);
  }

  // Optional emissive pass (e.g., night-side city lights).
  if (uUseEmissive > 0.5 && uEmissiveStrength > 0.0) {
    vec3 e = texture(uEmissiveTex, vUv).rgb;
    // Reveal emissive primarily on the night side; fade across the terminator.
    float night = smoothstep(uEmissiveNightFadeStart, uEmissiveNightFadeEnd, -ndotl);
    col += e * (uEmissiveStrength * night);
  }


  float a = 1.0;
  if (uAlphaFromTex > 0.5) {
    a = clamp(t.a * uAlphaMul, 0.0, 1.0);
  }
  FragColor = vec4(col, a);
}
)GLSL";

bool MeshRenderer::init(std::string* outError) {
  if (!shader_.build(kVS, kFS, outError)) return false;

  gl::GenBuffers(1, &instanceVbo_);

  // Streamed instance buffer (grown on demand).
  gl::BindBuffer(GL_ARRAY_BUFFER, instanceVbo_);
  gl::BufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STREAM_DRAW);
  instanceVboCapacityBytes_ = 0;

  // default identity matrices
  for (int i = 0; i < 16; ++i) {
    view_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
    proj_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
  }

  return true;
}

void MeshRenderer::ensureInstanceAttribLayoutBound_() {
  if (!mesh_) return;
  if (instanceLayoutBoundForMesh_ == mesh_) return;

  // We assume the mesh VAO is already bound via mesh_->bind(). Instance
  // attribute pointers are VAO state, so we only need to set them once per
  // mesh.
  gl::BindBuffer(GL_ARRAY_BUFFER, instanceVbo_);

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

  instanceLayoutBoundForMesh_ = mesh_;
}

void MeshRenderer::setViewProj(const float* view, const float* proj) {
  std::memcpy(view_, view, sizeof(float) * 16);
  std::memcpy(proj_, proj, sizeof(float) * 16);

  // Cache the frustum once per view/projection update (used by drawInstancesCulled).
  stellar::math::Mat4d v{};
  stellar::math::Mat4d p{};
  for (int i = 0; i < 16; ++i) {
    v.m[i] = static_cast<double>(view_[i]);
    p.m[i] = static_cast<double>(proj_[i]);
  }
  frustum_ = stellar::math::frustumFromViewProjection(v, p);
  frustumValid_ = frustum_.isFinite();
}

void MeshRenderer::drawInstances(const std::vector<InstanceData>& instances) {
  if (!mesh_ || instances.empty()) return;

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);
  shader_.setUniform1i("uTex", 0);
  shader_.setUniform1i("uNormalTex", 1);
  shader_.setUniform1f("uUseNormalMap", normalTex_ ? 1.0f : 0.0f);
  shader_.setUniform1f("uNormalStrength", normalStrength_);
  shader_.setUniform1i("uEmissiveTex", 2);
  shader_.setUniform1f("uUseEmissive", emissiveTex_ ? 1.0f : 0.0f);
  shader_.setUniform1f("uEmissiveStrength", emissiveStrength_);
  shader_.setUniform1f("uEmissiveNightFadeStart", emissiveNightFadeStart_);
  shader_.setUniform1f("uEmissiveNightFadeEnd", emissiveNightFadeEnd_);
  shader_.setUniform1f("uUnlit", unlit_ ? 1.0f : 0.0f);
  shader_.setUniform3f("uLightPos", lightPos_[0], lightPos_[1], lightPos_[2]);
  shader_.setUniform3f("uCameraPos", cameraPos_[0], cameraPos_[1], cameraPos_[2]);
  shader_.setUniform1f("uSpecularStrength", specularStrength_);
  shader_.setUniform1f("uShininess", shininess_);
  shader_.setUniform1f("uAlphaFromTex", alphaFromTexture_ ? 1.0f : 0.0f);
  shader_.setUniform1f("uAlphaMul", alphaMul_);

  if (tex_) tex_->bind(0);
  if (normalTex_) normalTex_->bind(1);
  if (emissiveTex_) emissiveTex_->bind(2);
  else if (tex_) tex_->bind(2);

  // Bind mesh VAO and ensure its instance attribute layout is bound once.
  mesh_->bind();
  ensureInstanceAttribLayoutBound_();

  gl::BindBuffer(GL_ARRAY_BUFFER, instanceVbo_);

  const std::size_t bytes = instances.size() * sizeof(InstanceData);
  if (bytes > instanceVboCapacityBytes_) {
    // Grow once, then stream into the existing store to avoid per-draw reallocations.
    gl::BufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(bytes), nullptr, GL_STREAM_DRAW);
    instanceVboCapacityBytes_ = bytes;
  }

  bool wrote = false;
  if (gl::MapBufferRange) {
    void* dst = gl::MapBufferRange(GL_ARRAY_BUFFER,
                                  0,
                                  static_cast<GLsizeiptr>(bytes),
                                  GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_RANGE_BIT);
    if (dst) {
      std::memcpy(dst, instances.data(), bytes);
      gl::UnmapBuffer(GL_ARRAY_BUFFER);
      wrote = true;
    }
  }
  if (!wrote) {
    gl::BufferSubData(GL_ARRAY_BUFFER, 0, static_cast<GLsizeiptr>(bytes), instances.data());
  }

  mesh_->drawInstanced(static_cast<std::uint32_t>(instances.size()));
}

void MeshRenderer::drawInstancesCulled(const std::vector<InstanceData>& instances) {
  if (!mesh_ || instances.empty()) return;

  // Small batches aren't worth the CPU cost.
  constexpr std::size_t kMinCullCount = 64;
  if (instances.size() < kMinCullCount) {
    drawInstances(instances);
    return;
  }

  if (!frustumValid_ || !mesh_->hasLocalBounds()) {
    drawInstances(instances);
    return;
  }

  const auto& localBounds = mesh_->localBounds();
  const stellar::math::Vec3d localCenter = localBounds.center();
  const double localRadius = mesh_->localBoundingRadius();

  culledScratch_.clear();
  culledScratch_.reserve(instances.size());

  for (const auto& inst : instances) {
    const double sx = static_cast<double>(inst.sx);
    const double sy = static_cast<double>(inst.sy);
    const double sz = static_cast<double>(inst.sz);

    const double maxScale = std::max({std::abs(sx), std::abs(sy), std::abs(sz)});
    const double radius = localRadius * maxScale;

    // World center = T + R * (S * localCenter)
    const stellar::math::Vec3d scaledCenter{localCenter.x * sx, localCenter.y * sy, localCenter.z * sz};

    stellar::math::Quatd q{static_cast<double>(inst.qw),
                           static_cast<double>(inst.qx),
                           static_cast<double>(inst.qy),
                           static_cast<double>(inst.qz)};

    if (!std::isfinite(q.w) || !std::isfinite(q.x) || !std::isfinite(q.y) || !std::isfinite(q.z)) {
      q = stellar::math::Quatd::identity();
    }
    q = q.normalized();

    const stellar::math::Vec3d worldCenter = q.rotate(scaledCenter) +
                                             stellar::math::Vec3d{static_cast<double>(inst.px),
                                                                  static_cast<double>(inst.py),
                                                                  static_cast<double>(inst.pz)};

    if (frustum_.intersectsSphere(worldCenter, radius)) {
      culledScratch_.push_back(inst);
    }
  }

  if (!culledScratch_.empty()) {
    drawInstances(culledScratch_);
  }
}

} // namespace stellar::render
