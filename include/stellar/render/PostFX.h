#pragma once

#include "stellar/render/Shader.h"

#include <string>

namespace stellar::render {

struct PostFXSettings {
  bool enabled{true};

  // Bloom
  bool bloomEnabled{true};
  float bloomThreshold{1.05f};   // HDR threshold
  float bloomKnee{0.65f};        // soft threshold range
  float bloomBoost{1.10f};       // scales extracted highlights
  float bloomIntensity{0.85f};   // scales blurred bloom when composited
  int bloomPasses{8};            // ping-pong blur passes (each pass is one direction)

  // Tonemap / output
  float exposure{1.00f};
  float gamma{2.20f};

  // Extra film-ish effects
  float vignette{0.18f};         // 0..1
  float grain{0.035f};           // 0..1
  float chromaticAberration{0.0015f}; // UV offset magnitude (0 disables)

  // Screen-space "warp streaks" (cheap speed effect). 0 disables.
  float warp{0.0f};
};

// HDR scene render target + post-processing (bloom + tonemap) for the game.
// Intended to be lightweight and hackable.
class PostFX {
public:
  PostFX() = default;
  ~PostFX();

  PostFX(const PostFX&) = delete;
  PostFX& operator=(const PostFX&) = delete;

  bool init(std::string* outError = nullptr);

  // Recreate render targets if the viewport size changes.
  void ensureSize(int w, int h);

  // Bind the HDR scene framebuffer as the current draw target.
  // Call this before rendering the 3D scene.
  void beginScene(int w, int h) const;

  // Run post-processing and draw to the default framebuffer (screen).
  // Call this after the 3D scene is rendered.
  void present(int w, int h, const PostFXSettings& settings, float timeSeconds);

  int width() const { return w_; }
  int height() const { return h_; }

private:
  void destroy();
  void createOrResize(int w, int h);

  void drawFullscreen() const;

  int w_{0};
  int h_{0};
  int bloomW_{0};
  int bloomH_{0};

  unsigned int sceneFbo_{0};
  unsigned int sceneTex_{0};
  unsigned int depthRbo_{0};

  unsigned int brightFbo_{0};
  unsigned int brightTex_{0};

  unsigned int pingFbo_[2]{0, 0};
  unsigned int pingTex_[2]{0, 0};

  unsigned int vao_{0};
  unsigned int vbo_{0};

  ShaderProgram bright_{};
  ShaderProgram blur_{};
  ShaderProgram composite_{};
};

} // namespace stellar::render
