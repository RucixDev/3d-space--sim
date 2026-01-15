#pragma once

#include "stellar/render/MeshRenderer.h" // InstanceData
#include "stellar/render/Mesh.h"
#include "stellar/render/Shader.h"

#include <string>
#include <vector>

namespace stellar::render {

// A lightweight procedural star corona renderer.
//
// Renders an additive emissive shell (typically a slightly scaled UV-sphere) with:
//  - Fresnel-style rim brightening (so it reads as a corona/halo rather than a foggy overlay)
//  - seeded, seam-free FBM noise on the unit sphere to break up the rim
//  - animated "streamers" and "prominence" lobes aligned to a deterministic rotation axis
//
// The intent is *not* physical accuracy; this is an art-directable, cheap cue that makes
// the system's star feel alive and HDR/bloom-friendly.
//
// Usage (caller controls render state):
//  - depth test enabled
//  - depth writes disabled
//  - additive blend (SRC_ALPHA, ONE)
//
// Per-instance color (InstanceData::cr/cg/cb) acts as the corona tint.
struct StarCoronaSettings {
  int seed{1};
  float intensity{1.0f};

  // Rim thickness: higher = thinner/Sharper limb.
  float rimPower{4.8f};

  // 3D noise sampled on the unit sphere (direction). Frequency controls feature size.
  float noiseFrequency{3.25f};
  float noiseStrength{0.75f};

  // Animated flow speed (0 disables animation).
  float flowSpeed{0.12f};

  // Streamer rays around the limb.
  float streamerStrength{0.55f};
  float streamerSharpness{5.0f};

  // Prominence lobes (equatorial-ish arcs). Count is clamped in shader.
  float prominenceStrength{0.65f};
  int prominenceCount{7};
  float prominenceSharpness{3.2f};

  // Global brightness pulsation.
  float pulseStrength{0.08f};
};

class StarCoronaRenderer {
public:
  StarCoronaRenderer() = default;
  ~StarCoronaRenderer();

  StarCoronaRenderer(const StarCoronaRenderer&) = delete;
  StarCoronaRenderer& operator=(const StarCoronaRenderer&) = delete;

  bool init(std::string* outError = nullptr);

  void setMesh(const Mesh* mesh) { mesh_ = mesh; }
  void setViewProj(const float* view, const float* proj);
  void setCameraPos(float x, float y, float z) {
    camPos_[0] = x;
    camPos_[1] = y;
    camPos_[2] = z;
  }

  void setTimeSeconds(float t) { timeSeconds_ = t; }

  void setSettings(const StarCoronaSettings& s) { settings_ = s; }
  const StarCoronaSettings& settings() const { return settings_; }

  void drawInstances(const std::vector<InstanceData>& instances);

private:
  const Mesh* mesh_{nullptr};
  ShaderProgram shader_{};
  unsigned int instanceVbo_{0};

  float view_[16]{};
  float proj_[16]{};
  float camPos_[3]{0.0f, 0.0f, 0.0f};
  float timeSeconds_{0.0f};

  StarCoronaSettings settings_{};
};

} // namespace stellar::render
