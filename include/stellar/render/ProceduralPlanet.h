#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/Texture.h"

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace stellar::render {

// Procedural surface textures for spherical bodies.
//
// These textures are designed for a standard UV sphere with equirectangular UVs:
//   - U = longitude (0..1)
//   - V = latitude  (0..1)
//
// The generator uses 3D noise evaluated on the unit sphere so the texture is
// seam-free at U=0/1.

enum class SurfaceKind : core::u8 {
  Rocky = 0,
  Desert,
  Ocean,
  Ice,
  GasGiant,
  Star,
  // Cloud layer alpha mask (RGB is near-white, A encodes density).
  // Intended to be rendered as a slightly larger UV-sphere shell with alpha blending.
  Clouds
};

struct SurfaceImage {
  int w{0};
  int h{0};
  std::vector<std::uint8_t> rgba; // size = w*h*4
};

// Generate a deterministic equirectangular RGBA8 texture for a given surface kind.
// widthPx controls the output width; height is width/2 (2:1 equirectangular).
//
// For SurfaceKind::Clouds, RGB is near-white and alpha encodes cloud density.
// Render the resulting texture as a slightly larger shell around the planet.
SurfaceImage generateSurfaceTexture(SurfaceKind kind, core::u64 seed, int widthPx);

// Generate a deterministic tangent-space normal map (RGBA8) for a given surface kind.
// The RGB channels encode XYZ in tangent space (mapped from [-1,1] to [0,255]).
//
// Note: strength/scaling is applied in the shader (MeshRenderer) so this output stays cacheable.
SurfaceImage generateSurfaceNormalMap(SurfaceKind kind, core::u64 seed, int widthPx);

// Small runtime cache that uploads generated surface textures as OpenGL textures.
// Intended for stars/planets (low count) and UI previews.
class SurfaceTextureCache {
public:
  SurfaceTextureCache() = default;

  void clear();

  void setMaxEntries(std::size_t max) { maxEntries_ = max; }
  std::size_t maxEntries() const { return maxEntries_; }
  std::size_t size() const { return cache_.size(); }

  // Get a GPU texture for (kind, seed). Generates and uploads if needed.
  // widthPx controls the texture resolution (height is width/2).
  const Texture2D& get(SurfaceKind kind, core::u64 seed, int widthPx = 512);

private:
  struct Entry {
    Texture2D tex{};
    core::u64 lastUseTick{0};
  };

  std::unordered_map<core::u64, Entry> cache_;
  core::u64 tick_{0};
  std::size_t maxEntries_{64};

  core::u64 makeKey(SurfaceKind kind, core::u64 seed, int widthPx) const;
  void evictIfNeeded();
};

// Small runtime cache that uploads generated surface normal maps as OpenGL textures.
// Normal maps are encoded as RGBA8 tangent-space normals and are intended to be sampled
// with mipmapping and linear filtering.
//
// This cache is separate from SurfaceTextureCache so albedo and normals can be managed
// independently (and disabled at runtime).
class SurfaceNormalMapCache {
public:
  SurfaceNormalMapCache() = default;

  void clear();

  void setMaxEntries(std::size_t max) { maxEntries_ = max; }
  std::size_t maxEntries() const { return maxEntries_; }
  std::size_t size() const { return cache_.size(); }

  // Get a GPU texture for (kind, seed). Generates and uploads if needed.
  // widthPx controls the texture resolution (height is width/2).
  const Texture2D& get(SurfaceKind kind, core::u64 seed, int widthPx = 512);

private:
  struct Entry {
    Texture2D tex{};
    core::u64 lastUseTick{0};
  };

  std::unordered_map<core::u64, Entry> cache_;
  core::u64 tick_{0};
  std::size_t maxEntries_{64};

  core::u64 makeKey(SurfaceKind kind, core::u64 seed, int widthPx) const;
  void evictIfNeeded();
};

} // namespace stellar::render
