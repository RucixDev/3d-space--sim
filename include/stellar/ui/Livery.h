#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/ProceduralLivery.h"

#include <string>

namespace stellar::ui {

// User-facing livery configuration (saved to a small text file).
//
// This is intentionally separate from the procedural generator's descriptor so we can
// evolve the UI without breaking the render module API.
struct LiveryConfig {
  int version{1};

  render::LiveryPattern pattern{render::LiveryPattern::Stripes};
  core::u64 seed{1337};

  float base[3]{0.22f, 0.26f, 0.33f};
  float accent1[3]{0.92f, 0.52f, 0.18f};
  float accent2[3]{0.86f, 0.86f, 0.88f};

  float scale{1.0f};
  float angleDeg{35.0f};
  float detail{0.65f};
  float wear{0.15f};
  float contrast{0.85f};
  bool decal{true};

  // Texture resolution for the generated livery.
  int textureSizePx{512};

  // If false, the in-world ship falls back to the checker texture (useful for debugging).
  bool applyInWorld{true};
};

std::string defaultLiveryPath();
LiveryConfig makeDefaultLivery();

bool saveLiveryToFile(const LiveryConfig& cfg, const std::string& path);
bool loadLiveryFromFile(const std::string& path, LiveryConfig& out);

} // namespace stellar::ui
