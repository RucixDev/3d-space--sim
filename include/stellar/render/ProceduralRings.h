#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/Texture.h"

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace stellar::render {

// Simple procedural planet ring texture generator.
//
// The texture is intended for an annulus mesh with UVs:
//   u = angle around the ring  [0..1]
//   v = radial coordinate      [0..1] (inner->outer)
//
// The generator is deterministic per (seed, widthPx, heightPx).
struct RingImage {
  int w{0};
  int h{0};
  std::vector<std::uint8_t> rgba; // RGBA8 row-major
};

RingImage generateRingTexture(core::u64 seed, int widthPx = 1024, int heightPx = 256);

// Small LRU-ish cache for GL textures generated from procedural ring images.
//
// Note: This is intentionally tiny and simple (linear eviction). If rings become
// a hot path, upgrade to an actual LRU list.
class RingTextureCache {
public:
  RingTextureCache() = default;

  void clear();

  std::size_t size() const { return cache_.size(); }
  std::size_t maxEntries() const { return maxEntries_; }
  void setMaxEntries(std::size_t m);

  const Texture2D& get(core::u64 seed, int widthPx = 1024, int heightPx = 256);

private:
  struct Entry {
    Texture2D tex;
    core::u64 lastUseTick{0};
  };

  core::u64 tick_{0};
  std::size_t maxEntries_{64};
  std::unordered_map<core::u64, Entry> cache_;

  static core::u64 makeKey(core::u64 seed, int w, int h);
  void evictIfNeeded();
};

} // namespace stellar::render
