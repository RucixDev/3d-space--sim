#pragma once

#include "stellar/core/Types.h"

#include "stellar/render/ProceduralPlanet.h" // SurfaceKind
#include "stellar/render/Shader.h"
#include "stellar/render/Texture.h"

#include <cstddef>
#include <string>
#include <unordered_map>

namespace stellar::render {

// GPU-side procedural surface texture generator + cache.
//
// This is an alternative to the CPU SurfaceTextureCache/SurfaceNormalMapCache.
// Instead of generating RGBA pixel buffers on the CPU and uploading them, we
// render procedurally into an RGBA8 texture via an FBO + fragment shader.
//
// The core design goal is to make high-resolution procedural surfaces (and their
// matching normal maps) much cheaper to generate, while keeping the same public
// "get texture for (kind, seed, width)" API.
class GpuSurfaceCache {
public:
  GpuSurfaceCache() = default;
  ~GpuSurfaceCache();

  GpuSurfaceCache(const GpuSurfaceCache&) = delete;
  GpuSurfaceCache& operator=(const GpuSurfaceCache&) = delete;

  // Must be called after an OpenGL context is created.
  bool init(std::string* outError = nullptr);
  bool isInited() const { return inited_; }

  void clear();

  void setMaxEntries(std::size_t max) { maxEntries_ = max; }
  std::size_t maxEntries() const { return maxEntries_; }
  std::size_t albedoSize() const { return albedoCache_.size(); }
  std::size_t normalSize() const { return normalCache_.size(); }

  // Get an albedo texture (RGBA8) for (kind, seed). Generates on demand.
  const Texture2D& albedo(SurfaceKind kind, core::u64 seed, int widthPx = 512);

  // Get a tangent-space normal map (RGBA8) for (kind, seed). Generates on demand.
  const Texture2D& normal(SurfaceKind kind, core::u64 seed, int widthPx = 512);

private:
  struct Entry {
    Texture2D tex{};
    core::u64 lastUseTick{0};
  };

  void destroy();

  core::u64 makeKey(const char* tag, SurfaceKind kind, core::u64 seed, int widthPx) const;
  void evictIfNeeded(std::unordered_map<core::u64, Entry>& cache);

  // mode: 0=albedo, 1=normal
  bool render(Texture2D& out, SurfaceKind kind, core::u64 seed, int w, int h, int mode);

  bool inited_{false};
  unsigned int fbo_{0};
  unsigned int vao_{0};
  ShaderProgram shader_{};

  std::unordered_map<core::u64, Entry> albedoCache_;
  std::unordered_map<core::u64, Entry> normalCache_;
  core::u64 tick_{0};
  std::size_t maxEntries_{96};
};

} // namespace stellar::render
