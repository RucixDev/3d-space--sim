#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/Mesh.h"
#include "stellar/render/Shader.h"
#include "stellar/render/Texture.h"

#include <string>
#include <vector>

namespace stellar::render {

// Volumetric single-scattering atmosphere renderer.
//
// This is an optional, more physically-motivated alternative to the
// AtmosphereRenderer limb-glow shader. It raymarches through a spherical shell
// and integrates a cheap single-scattering approximation using an exponential
// density profile.
//
// Design goals:
//  - deterministic & self-contained (no external LUTs required)
//  - fast enough for a prototype (low step counts)
//  - supports the project's spectral Mie phase LUT when available
//
// Usage pattern (caller controlled):
//  - depth test enabled
//  - depth writes disabled
//  - additive blend (SRC_ALPHA, ONE)
//
// Notes:
//  - units are "render units" (same as MeshRenderer / AtmosphereRenderer)
//  - the mesh is expected to be a unit sphere scaled per-instance
struct VolumetricAtmosphereInstance {
  // Center position (render units)
  float px, py, pz;

  // Outer radius scale (render units). For a unit sphere mesh, sx==sy==sz==Rout.
  float sx, sy, sz;

  // Orientation quaternion (x,y,z,w). Unused for perfect spheres but kept for
  // parity with other instanced renderers.
  float qx, qy, qz, qw;

  // Scattering tint (roughly an albedo / color bias)
  float cr, cg, cb;

  // Inner radius (planet radius) in render units.
  float innerRadius;
  // Exponential scale height in render units.
  float scaleHeight;
  // Relative density multiplier (1 ~= Earth-ish). Used as a cheap thickness/
  // extinction scaler.
  float densityMul;
  // Relative Mie multiplier (1 default). Lets callers bias aerosol haze per body.
  float mieMul;
};

class VolumetricAtmosphereRenderer {
public:
  VolumetricAtmosphereRenderer() = default;
  ~VolumetricAtmosphereRenderer();

  VolumetricAtmosphereRenderer(const VolumetricAtmosphereRenderer&) = delete;
  VolumetricAtmosphereRenderer& operator=(const VolumetricAtmosphereRenderer&) = delete;

  bool init(std::string* outError = nullptr);

  void setMesh(const Mesh* mesh) { mesh_ = mesh; }

  void setViewProj(const float* view, const float* proj);
  void setCameraPos(float x, float y, float z) {
    camPos_[0] = x;
    camPos_[1] = y;
    camPos_[2] = z;
  }

  void setSunPos(float x, float y, float z) {
    sunPos_[0] = x;
    sunPos_[1] = y;
    sunPos_[2] = z;
  }

  void setIntensity(float v) { intensity_ = v; }
  void setRayleighStrength(float v) { rayleighStrength_ = v; }
  void setMieStrength(float v) { mieStrength_ = v; }
  void setScaleHeightMul(float v) { scaleHeightMul_ = v; }
  void setDensityMul(float v) { densityMul_ = v; }

  void setSteps(int viewSteps) { viewSteps_ = viewSteps; }
  void setSunSteps(int sunSteps) { sunSteps_ = sunSteps; }

  // Jittered sampling reduces banding. 0 disables.
  void setDitherStrength(float v) { ditherStrength_ = v; }

  // Fallback phase function when no LUT is bound.
  void setMieG(float g) { mieG_ = g; }

  // Optional: spectral Mie phase LUT (width=N samples over μ∈[-1,1], height=1).
  void setUseMiePhaseLut(bool v) { useMiePhaseLut_ = v; }
  void setMiePhaseLut(const Texture2D* lut) { miePhaseLut_ = lut; }
  // Blends between HG fallback (0) and the LUT (1).
  void setMiePhaseStrength(float v) { miePhaseStrength_ = v; }

  // Analytic multiple scattering approximation (cheap).
  //
  // This adds a second in-scattering term shaped by a broader phase function
  // and boosted by a geometric-series energy approximation.
  void setMultipleScatteringStrength(float v) { msStrength_ = v; }
  void setMultipleScatteringAlbedo(float v) { msAlbedo_ = v; }

  // Optional: separate phase LUT for the multiple-scattering lobe.
  // If not set, a HG fallback with g^2 is used.
  void setUseMultipleScatteringPhaseLut(bool v) { useMsPhaseLut_ = v; }
  void setMultipleScatteringPhaseLut(const Texture2D* lut) { msPhaseLut_ = lut; }
  // Blends between HG fallback (0) and the LUT (1).
  void setMultipleScatteringPhaseStrength(float v) { msPhaseStrength_ = v; }

  void drawInstances(const std::vector<VolumetricAtmosphereInstance>& instances);

private:
  const Mesh* mesh_{nullptr};
  ShaderProgram shader_{};
  unsigned int instanceVbo_{0};

  float view_[16]{};
  float proj_[16]{};
  float camPos_[3]{0.0f, 0.0f, 0.0f};
  float sunPos_[3]{0.0f, 0.0f, 0.0f};

  float intensity_{1.0f};
  float rayleighStrength_{1.0f};
  float mieStrength_{1.0f};
  float scaleHeightMul_{1.0f};
  float densityMul_{1.0f};

  int viewSteps_{10};
  int sunSteps_{6};
  float ditherStrength_{0.35f};
  float mieG_{0.80f};

  bool useMiePhaseLut_{false};
  const Texture2D* miePhaseLut_{nullptr};
  float miePhaseStrength_{1.0f};

  // Multiple scattering approximation.
  float msStrength_{0.0f};
  float msAlbedo_{0.92f};
  bool useMsPhaseLut_{false};
  const Texture2D* msPhaseLut_{nullptr};
  float msPhaseStrength_{1.0f};
};

} // namespace stellar::render
