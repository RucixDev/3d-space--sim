#pragma once

#include "stellar/core/Types.h"

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace stellar::render {

// Procedural ship/station livery textures.
//
// These are generic 2D RGBA8 textures intended for UV-mapped meshes.
// The prototype uses a cube mesh for ships/contacts, but the generator is
// written to be usable for any future mesh with meaningful UVs.

enum class LiveryPattern : core::u8 {
  Solid = 0,
  Stripes,
  Hazard,
  Camo,
  Hex,
  Digital,

  Count
};

const char* toString(LiveryPattern p);
std::optional<LiveryPattern> patternFromString(const std::string& s);

struct LiveryDesc {
  LiveryPattern pattern{LiveryPattern::Stripes};
  core::u64 seed{1};

  // Colors are in linear-ish 0..1 floats. (The renderer treats them as sRGB-ish
  // because we keep the shader extremely simple.)
  float base[3]{0.22f, 0.26f, 0.33f};
  float accent1[3]{0.92f, 0.52f, 0.18f};
  float accent2[3]{0.86f, 0.86f, 0.88f};

  // Pattern controls (interpretation depends on the selected pattern).
  float scale{1.0f};        // 0.25..4.0
  float angleDeg{35.0f};    // stripes / hazard
  float detail{0.65f};      // noise detail / cell size
  float wear{0.15f};        // 0..1 (scratches, dirt)
  float contrast{0.85f};    // 0..1

  // Optional "registration" decal.
  bool decal{true};
};

struct LiveryImage {
  int w{0};
  int h{0};
  std::vector<std::uint8_t> rgba; // size = w*h*4
};

// Generate a deterministic square livery texture (RGBA8).
// sizePx controls both width and height.
LiveryImage generateLiveryTexture(const LiveryDesc& desc, int sizePx);

} // namespace stellar::render
