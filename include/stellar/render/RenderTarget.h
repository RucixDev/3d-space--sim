#pragma once

#include "stellar/render/Texture.h"

#include <string>

namespace stellar::render {

// Minimal OpenGL render target that renders into an RGBA8 texture (+ optional depth buffer).
//
// Primary use-case: render small 3D previews (ship/hangar, material preview, etc.)
// into an offscreen texture that can be displayed in Dear ImGui with ImGui::Image.
struct RenderTarget2DDesc {
  // If true, allocate a depth/stencil renderbuffer (GL_DEPTH24_STENCIL8).
  // Set false for lightweight fullscreen passes (procedural graph baking, blits, etc.).
  bool hasDepth{true};

  // If true, allocate mipmaps for the color texture.
  // Note: the mip chain must be regenerated after rendering via RenderTarget2D::generateMips().
  bool generateMips{false};

  // If true, use nearest filtering (useful for pixel art / debug buffers).
  // Otherwise use linear filtering.
  bool nearestFilter{false};

  // If true, clamp UVs to edge. If false, repeat.
  bool clampToEdge{true};
};

class RenderTarget2D {
public:
  RenderTarget2D() = default;
  ~RenderTarget2D();

  RenderTarget2D(const RenderTarget2D&) = delete;
  RenderTarget2D& operator=(const RenderTarget2D&) = delete;

  // Backwards-compatible init (RGBA8 + depth).
  bool init(int w, int h, std::string* outError = nullptr);

  // Init with an explicit descriptor.
  bool init(int w, int h, const RenderTarget2DDesc& desc, std::string* outError = nullptr);

  // Resize preserving the current descriptor.
  void ensureSize(int w, int h);

  void begin() const;
  static void end();

  // If generateMips is enabled in the descriptor, call this after rendering to refresh the mip chain.
  void generateMips() const;

  const Texture2D& color() const { return color_; }
  int width() const { return w_; }
  int height() const { return h_; }

  const RenderTarget2DDesc& desc() const { return desc_; }

  bool isInited() const { return inited_; }

private:
  void destroy();
  bool createOrResize(int w, int h, const RenderTarget2DDesc& desc, std::string* outError = nullptr);

  unsigned int fbo_{0};
  unsigned int depthRbo_{0};
  Texture2D color_{};
  int w_{0};
  int h_{0};
  bool inited_{false};
  RenderTarget2DDesc desc_{};
};

} // namespace stellar::render
