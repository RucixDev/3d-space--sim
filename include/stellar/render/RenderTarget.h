#pragma once

#include "stellar/render/Texture.h"

#include <string>

namespace stellar::render {

// Minimal OpenGL render target that renders into an RGBA8 texture + depth buffer.
//
// Primary use-case: render small 3D previews (ship/hangar, material preview, etc.)
// into an offscreen texture that can be displayed in Dear ImGui with ImGui::Image.
class RenderTarget2D {
public:
  RenderTarget2D() = default;
  ~RenderTarget2D();

  RenderTarget2D(const RenderTarget2D&) = delete;
  RenderTarget2D& operator=(const RenderTarget2D&) = delete;

  bool init(int w, int h, std::string* outError = nullptr);
  void ensureSize(int w, int h);

  void begin() const;
  static void end();

  const Texture2D& color() const { return color_; }
  int width() const { return w_; }
  int height() const { return h_; }

private:
  void destroy();
  bool createOrResize(int w, int h, std::string* outError = nullptr);

  unsigned int fbo_{0};
  unsigned int depthRbo_{0};
  Texture2D color_{};
  int w_{0};
  int h_{0};
  bool inited_{false};
};

} // namespace stellar::render
