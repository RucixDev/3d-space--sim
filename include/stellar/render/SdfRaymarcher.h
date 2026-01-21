#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/RenderTarget.h"
#include "stellar/render/SdfGraph.h"
#include "stellar/render/Shader.h"

#include <array>
#include <string>

namespace stellar::render {

// GPU SDF previewer: builds a GLSL fragment shader from an SdfGraph and
// renders it via sphere tracing into an offscreen RenderTarget2D.
//
// This is intentionally a *preview* path: it is meant to be fast to iterate,
// not a perfect match for CPU meshing. The mesher still uses evalSdfGraph(...)
// from SdfGraph.cpp.

struct SdfRaymarchCamera {
  float pos[3]{0.0f, 0.0f, 3.0f};
  float forward[3]{0.0f, 0.0f, -1.0f};
  float right[3]{1.0f, 0.0f, 0.0f};
  float up[3]{0.0f, 1.0f, 0.0f};
  float fovYRadians{1.0f};
};

enum class SdfRaymarchDebug : core::u8 {
  Shaded = 0,
  Normals = 1,
  Steps = 2,
  Distance = 3,
};

struct SdfRaymarchSettings {
  // Core sphere-tracing parameters.
  int maxSteps{128};
  float epsilon{0.0015f};
  float maxDistance{25.0f};

  // Clip rays to an AABB around the origin to avoid wasting steps.
  bool useBoundsAabb{true};

  // Lighting.
  float lightDir[3]{0.40f, 0.75f, 0.20f};
  float baseColor[3]{0.85f, 0.90f, 1.00f};
  float backgroundColor[3]{0.02f, 0.02f, 0.03f};
  float ambient{0.20f};
  float diffuse{0.85f};
  float specular{0.18f};
  float shininess{64.0f};

  // Soft shadows (IQ-style). This is approximate and tuned for previews.
  bool softShadows{true};
  int shadowSteps{32};
  float shadowMaxDistance{8.0f};
  float shadowK{16.0f};

  // Ambient occlusion (cheap sample-based). Also approximate.
  bool ambientOcclusion{true};
  int aoSteps{5};
  float aoStepSize{0.12f};
  float aoStrength{1.0f};

  SdfRaymarchDebug debug{SdfRaymarchDebug::Shaded};
};

// Optional surface material inputs for the preview shader.
// When albedoTex is provided, shaded mode samples it using spherical UVs.
struct SdfRaymarchMaterial {
  const Texture2D* albedoTex{nullptr};
  // UV transform applied to the sampled coordinates.
  float uvScale{1.0f};
  float uvRotateDeg{0.0f};
  float uvOffset[2]{0.0f, 0.0f};

  // Projection / filtering.
  //
  // The old preview path selected a single dominant axis which can create visible
  // seams where the normal transitions between axes. Enabling triplanar blending
  // samples all three planar projections and blends them based on abs(normal).
  bool triplanarBlend{true};

  // Exponent applied to the abs(normal) weights before normalization.
  // Higher values reduce "layering" artifacts and make transitions tighter.
  float triplanarSharpness{8.0f};

  // Optional micro-normal perturbation derived from a height channel in the albedo
  // texture.
  //
  // This is intended for tooling / previews. It is *not* physically correct normal
  // mapping, but it adds convincing high-frequency detail when prototyping.
  bool heightFromAlpha{true};
  float microNormalStrength{0.0f};
  float microNormalStepTexels{1.0f};
};

struct SdfRaymarchStats {
  bool shaderRebuilt{false};
  double shaderBuildMs{0.0};
  double drawMs{0.0};
};

class SdfRaymarcher {
public:
  SdfRaymarcher() = default;
  ~SdfRaymarcher();

  SdfRaymarcher(const SdfRaymarcher&) = delete;
  SdfRaymarcher& operator=(const SdfRaymarcher&) = delete;

  // Render graph to an offscreen texture.
  // `bounds` is the half-extent of the preview AABB ([-bounds,+bounds]).
  // `iso` is provided so the preview matches the mesher's level-set.
  bool render(const SdfGraph& g,
              int width,
              int height,
              float timeSec,
              float bounds,
              float iso,
              const SdfRaymarchCamera& cam,
              const SdfRaymarchSettings& s,
              std::string* outError = nullptr,
              const SdfRaymarchMaterial* mat = nullptr);

  const Texture2D& texture() const { return target_.color(); }
  bool isReady() const { return target_.isInited() && shader_.handle() != 0; }
  const std::string& lastFragmentSource() const { return lastFragSrc_; }
  const SdfRaymarchStats& stats() const { return stats_; }

private:
  static core::u64 structureKey(const SdfGraph& g);
  static std::string buildFragmentShader(const SdfGraph& g);
  bool ensureShader(const SdfGraph& g, std::string* outError);

  RenderTarget2D target_{};
  ShaderProgram shader_{};
  core::u64 shaderKey_{0};
  unsigned int vao_{0};

  std::array<float, kSdfGraphMaxNodes * 8> packedParams_{};
  std::array<int, kSdfGraphMaxNodes> packedSeeds_{};

  SdfRaymarchStats stats_{};
  std::string lastFragSrc_{};
};

} // namespace stellar::render
