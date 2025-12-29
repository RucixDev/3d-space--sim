#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/Texture.h"

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace stellar::render {

enum class SpriteKind : core::u8 {
  Commodity = 0,
  Faction,
  Mission,
  Ship,
  Station,
  Planet,
  Star,

  // In-world objects / HUD markers.
  Cargo,
  Asteroid,
  Signal,

  // HUD / flight symbology.
  HudReticle,
  HudLead,
  HudVelocity
};

struct SpriteImage {
  int w{0};
  int h{0};
  std::vector<std::uint8_t> rgba; // size = w*h*4
};

// Generate a deterministic sprite image for a given kind/seed.
// The output is RGBA8 with a transparent background.
SpriteImage generateSprite(SpriteKind kind, core::u64 seed, int size);

// Small runtime cache that uploads generated sprites as OpenGL textures.
// Intended for UI icons (Dear ImGui) and simple debug visualizations.
class SpriteCache {
public:
  SpriteCache() = default;

  void clear();

  void setMaxEntries(std::size_t max) { maxEntries_ = max; }
  std::size_t maxEntries() const { return maxEntries_; }
  std::size_t size() const { return cache_.size(); }

  const Texture2D& get(SpriteKind kind, core::u64 seed, int size = 32);

private:
  struct Entry {
    Texture2D tex{};
    core::u64 lastUseTick{0};
  };

  std::unordered_map<core::u64, Entry> cache_;
  core::u64 tick_{0};
  std::size_t maxEntries_{512};

  void evictIfNeeded();
};

// UV rectangle in normalized texture coordinates.
// (u0,v0) is the lower-left corner, (u1,v1) is upper-right.
struct SpriteUvRect {
  float u0{0.0f};
  float v0{0.0f};
  float u1{1.0f};
  float v1{1.0f};
};

// A texture atlas for procedural sprites.
//
// Why an atlas?
// - HUD overlays can draw dozens/hundreds of markers per frame.
// - Binding a unique GL texture per marker (even if cached) can cause overhead.
// - An atlas batches all icons into a single texture while preserving
//   deterministic procedural generation.
//
// The atlas stores each sprite in a fixed-size cell with padding to reduce
// filtering bleed when scaled.
class SpriteAtlas {
public:
  SpriteAtlas() = default;

  void clear();

  // (Re)initialize atlas parameters. This clears any existing entries.
  // atlasSizePx: texture width/height in pixels (square)
  // cellSizePx:  sprite size per cell (excluding padding)
  // paddingPx:   padding replicated from sprite edges on all sides
  // nearestFilter: if true, uses nearest sampling (crisper UI); otherwise linear.
  void init(int atlasSizePx = 1024, int cellSizePx = 64, int paddingPx = 2, bool nearestFilter = true);

  void setMaxEntries(std::size_t max) { maxEntries_ = max; }
  std::size_t maxEntries() const { return maxEntries_; }
  std::size_t size() const { return entries_.size(); }

  int atlasSizePx() const { return atlasSizePx_; }
  int cellSizePx() const { return cellSizePx_; }
  int paddingPx() const { return paddingPx_; }
  int capacity() const { return capacity_; }

  const Texture2D& texture() const { return tex_; }

  // Get the UV rect for a given (kind, seed). Generates/uploads into the atlas
  // if needed. May evict least-recently-used entries when full.
  SpriteUvRect get(SpriteKind kind, core::u64 seed);

private:
  struct Entry {
    SpriteUvRect uv{};
    int cellIndex{0};
    core::u64 lastUseTick{0};
  };

  Texture2D tex_{};
  bool inited_{false};

  int atlasSizePx_{1024};
  int cellSizePx_{64};
  int paddingPx_{2};
  int stridePx_{68};
  int cols_{0};
  int rows_{0};
  int capacity_{0};

  std::vector<int> freeCells_{};
  std::unordered_map<core::u64, Entry> entries_{};

  core::u64 tick_{0};
  std::size_t maxEntries_{256};

  core::u64 atlasKey(SpriteKind kind, core::u64 seed) const;
  int allocCell();
  void uploadCell(int cellIndex, SpriteKind kind, core::u64 seed);
  void evictIfNeeded();
};

} // namespace stellar::render
