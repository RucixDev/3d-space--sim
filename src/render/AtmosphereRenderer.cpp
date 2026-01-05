#include "stellar/render/AtmosphereRenderer.h"

#include "stellar/render/Gl.h"

#include <cstring>

namespace stellar::render {

AtmosphereRenderer::~AtmosphereRenderer() {
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
}
)GLSL";

static const char* kFS = R"GLSL(
#version 330 core
in vec3 vColor;
in vec3 vNrm;
in vec3 vWorldPos;

uniform vec3 uCamPos;
uniform vec3 uSunPos;
uniform float uIntensity;
uniform float uPower;
uniform float uSunLitBoost;
uniform float uForwardScatter;

out vec4 FragColor;

void main() {
  vec3 n = normalize(vNrm);
  vec3 v = normalize(uCamPos - vWorldPos);

  // Sun position is provided by the caller.
  vec3 l = normalize(uSunPos - vWorldPos);

  float ndotv = clamp(dot(n, v), 0.0, 1.0);
  float ndotl = clamp(dot(n, l), 0.0, 1.0);

  // Fresnel-style rim factor.
  float rim = pow(1.0 - ndotv, max(0.25, uPower));

  // Make the rim brighter on the day-side limb.
  float sun = mix(1.0, 0.35 + 0.65 * ndotl, clamp(uSunLitBoost, 0.0, 1.0));

  // Forward-scatter highlight when looking roughly toward the sun.
  float forward = pow(clamp(dot(v, -l), 0.0, 1.0), 6.0);

  float a = rim * sun + forward * uForwardScatter * rim;
  a = clamp(a, 0.0, 1.0);

  // Additive blend: contribution = rgb * alpha.
  FragColor = vec4(vColor * uIntensity, a);
}
)GLSL";

bool AtmosphereRenderer::init(std::string* outError) {
  if (!shader_.build(kVS, kFS, outError)) return false;

  gl::GenBuffers(1, &instanceVbo_);

  // Default identity matrices.
  for (int i = 0; i < 16; ++i) {
    view_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
    proj_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
  }

  return true;
}

void AtmosphereRenderer::setViewProj(const float* view, const float* proj) {
  std::memcpy(view_, view, sizeof(float) * 16);
  std::memcpy(proj_, proj, sizeof(float) * 16);
}

void AtmosphereRenderer::drawInstances(const std::vector<InstanceData>& instances) {
  if (!mesh_ || instances.empty()) return;

  // Additive limb glow
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE);

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);
  shader_.setUniform3f("uCamPos", camPos_[0], camPos_[1], camPos_[2]);
  shader_.setUniform3f("uSunPos", sunPos_[0], sunPos_[1], sunPos_[2]);
  shader_.setUniform1f("uIntensity", intensity_);
  shader_.setUniform1f("uPower", power_);
  shader_.setUniform1f("uSunLitBoost", sunLitBoost_);
  shader_.setUniform1f("uForwardScatter", forwardScatter_);

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
