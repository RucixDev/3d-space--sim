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
  Star
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

} // namespace stellar::render
