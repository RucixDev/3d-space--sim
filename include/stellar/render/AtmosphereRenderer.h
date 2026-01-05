#pragma once

#include "stellar/render/MeshRenderer.h" // InstanceData
#include "stellar/render/Shader.h"

#include <string>
#include <vector>

namespace stellar::render {

// A tiny additive atmosphere / limb glow renderer.
//
// It draws instanced UV-spheres (or any PN mesh) with a Fresnel-ish rim term and
// an optional sun-lit boost (caller supplies the sun position; the game sets it to the main star).
//
// Usage pattern (caller controls when to draw):
//   - depth test enabled
//   - depth writes disabled
//   - additive blend (SRC_ALPHA, ONE)
//
// The per-instance color (InstanceData::cr/cg/cb) is treated as the atmosphere tint.
class AtmosphereRenderer {
public:
  AtmosphereRenderer() = default;
  ~AtmosphereRenderer();

  AtmosphereRenderer(const AtmosphereRenderer&) = delete;
  AtmosphereRenderer& operator=(const AtmosphereRenderer&) = delete;

  bool init(std::string* outError = nullptr);

  void setMesh(const Mesh* mesh) { mesh_ = mesh; }

  void setViewProj(const float* view, const float* proj);
  void setCameraPos(float x, float y, float z) {
    camPos_[0] = x;
    camPos_[1] = y;
    camPos_[2] = z;
  }

  // Sun/light position in world/render units.
  // (Default is the origin for backwards compatibility.)
  void setSunPos(float x, float y, float z) {
    sunPos_[0] = x;
    sunPos_[1] = y;
    sunPos_[2] = z;
  }

  // Scales the added atmosphere color (HDR-friendly).
  void setIntensity(float v) { intensity_ = v; }
  // Fresnel power exponent; higher = thinner rim.
  void setPower(float v) { power_ = v; }
  // How strongly sun lighting modulates the rim (0..1).
  void setSunLitBoost(float v) { sunLitBoost_ = v; }
  // Forward-scatter highlight when looking roughly toward the sun (0..1).
  void setForwardScatter(float v) { forwardScatter_ = v; }

  void drawInstances(const std::vector<InstanceData>& instances);

private:
  const Mesh* mesh_{nullptr};
  ShaderProgram shader_{};
  unsigned int instanceVbo_{0};

  float view_[16]{};
  float proj_[16]{};
  float camPos_[3]{0.0f, 0.0f, 0.0f};

  float sunPos_[3]{0.0f, 0.0f, 0.0f};

  float intensity_{1.0f};
  float power_{5.0f};
  float sunLitBoost_{0.85f};
  float forwardScatter_{0.25f};
};

} // namespace stellar::render
