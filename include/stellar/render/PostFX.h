#pragma once

#include "stellar/render/FrameGraph.h"
#include "stellar/render/Shader.h"
#include "stellar/render/AutoExposure.h"

#include <array>
#include <string>
#include <string_view>

namespace stellar::render {

// HDR scene render target format.
//
// R11G11B10F is a packed floating-point format (no alpha) that uses roughly
// half the bandwidth of RGBA16F. On integrated / bandwidth-limited GPUs this
// can significantly improve FPS.
//
// If the chosen format is not supported as an FBO color attachment, PostFX will
// automatically fall back to RGBA16F at runtime.
enum class HdrBufferFormat : int {
  Rgba16f = 0,
  R11G11B10f = 1,
};

struct PostFXSettings {
  bool enabled{true};

  // ---------------------------------------------------------------------------
  // Performance: render-scale / dynamic resolution
  // ---------------------------------------------------------------------------
  // Render the 3D scene + post-FX at a lower resolution and upscale in the
  // composite pass. Values < 1.0 reduce fill-rate and can significantly improve
  // FPS on integrated GPUs.
  //
  // NOTE: the UI runs at full resolution (ImGui is drawn after post-FX).
  float renderScale{1.0f}; // 0.25 .. 1.0

  // HDR buffer format used for the internal scene target and HDR post-FX.
  // The packed R11G11B10F format is typically faster on bandwidth-limited GPUs.
  HdrBufferFormat hdrBufferFormat{HdrBufferFormat::R11G11B10f};

  // If enabled, PostFX will automatically adjust the internal renderScale to
  // target a desired frame rate using a smoothed CPU frame time estimate.
  bool  dynamicResolutionEnabled{false};
  float dynamicResolutionTargetFps{60.0f};
  float dynamicResolutionMinScale{0.55f};
  float dynamicResolutionMaxScale{1.0f};

  // How aggressively we change resolution when we miss the target.
  // 1.0 = default, higher reacts faster.
  float dynamicResolutionResponse{1.0f};

  // Dead-zone around the target (fraction). Prevents oscillation.
  float dynamicResolutionHysteresis{0.06f};

  // Quantize the scale to reduce realloc churn when the scale changes.
  // (Typical values: 0.01..0.05)
  float dynamicResolutionStep{0.02f};

  // Optional upscale sharpening (0 disables). Useful when renderScale < 1.
  float upsampleSharpen{0.0f}; // 0..1

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

  // Auto exposure (optional)
  //
  // When enabled, PostFX estimates the average scene luminance and computes a
  // smooth exposure multiplier that is applied on top of `exposure`.
  bool autoExposureEnabled{false};
  float autoExposureKey{0.18f};      // target middle-grey
  float autoExposureMin{0.25f};      // clamp multiplier
  float autoExposureMax{4.00f};
  float autoExposureSpeedUp{2.50f};  // 1/sec (brighten)
  float autoExposureSpeedDown{1.00f}; // 1/sec (darken)
  float autoExposureCenterWeight{0.65f}; // 0..1 (0 = uniform metering)
  int autoExposureMaxSize{256};      // cap luminance buffer to keep cost bounded
  bool autoExposureAsyncReadback{true}; // use async GPU readback (PBO) to avoid stalls

  // Extra film-ish effects
  float vignette{0.18f};         // 0..1
  float grain{0.035f};           // 0..1
  float chromaticAberration{0.0015f}; // UV offset magnitude (0 disables)

  // Retro / CRT style compositor (optional).
  //
  // This is a stylized alternative to the default composite path.
  // It performs optional pixelation (sample the HDR scene on a coarse grid),
  // then applies ordered dithering + color quantization and optional CRT-ish
  // embellishments.
  bool retroEnabled{false};

  // Pixel block size in screen pixels.
  // 1 = full resolution (no pixelation), 2 = half-res, etc.
  int retroPixelSize{1};

  // Quantization levels per channel.
  // 256 ~= "no quantization". Typical values: 8, 16, 32, 64.
  int retroColorSteps{32};

  // Ordered dithering strength (0 disables).
  float retroDitherStrength{0.65f};

  // CRT-ish extras (0 disables each)
  float retroScanlines{0.0f};   // 0..1
  float retroCurvature{0.0f};   // 0..1 barrel distortion
  float retroJitter{0.0f};      // 0..1 horizontal per-scanline jitter

  // Screen-space "warp streaks" (cheap speed effect). 0 disables.
  float warp{0.0f};

  // Hyperspace / jump tunnel overlay (screen-space). 0 disables.
  // This is intended to be driven by gameplay (FSD state), but can also be tuned manually.
  float hyperspace{0.0f};        // 0..1
  float hyperspaceTwist{0.45f};  // 0..1
  float hyperspaceDensity{0.65f}; // 0..1
  float hyperspaceNoise{0.35f};  // 0..1
  float hyperspaceIntensity{1.15f}; // multiplier
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

  // Apply settings that affect internal render targets (render-scale, etc.).
  // This should be called once per-frame *before* ensureSize/beginScene so the
  // HDR scene FBO matches the currently selected render scale.
  void configure(const PostFXSettings& settings);

  // Recreate render targets if the viewport size changes.
  void ensureSize(int w, int h);

  // Bind the HDR scene framebuffer as the current draw target.
  // Call this before rendering the 3D scene.
  void beginScene(int w, int h) const;

  // Run post-processing and draw to the default framebuffer (screen).
  // Call this after the 3D scene is rendered.
  void present(int w, int h, const PostFXSettings& settings, float timeSeconds);

  // --- Composite shader variants ------------------------------------------------
  //
  // The composite pass is where bloom is combined with the scene, tonemapping is
  // applied, and stylized screen-space effects are layered on top.
  //
  // For "procedural graphics engine" experiments we allow building shader
  // *variants* at runtime by injecting `#define` lines into the built-in composite
  // fragment shader.
  //
  // `defines` should be a string containing one or more `#define ...` lines.
  // It must NOT include a `#version` directive.
  bool buildCompositeVariant(std::string_view defines, std::string* outError = nullptr);

  // Restore the built-in composite fragment shader.
  bool resetCompositeShader(std::string* outError = nullptr);

  // Debug helpers for UI / logging.
  bool hasCustomComposite() const { return hasCustomComposite_; }
  const std::string& compositeFragmentSource() const { return compositeFsBuilt_; }

  // Auto-exposure debug values (only meaningful when auto-exposure is enabled).
  float lastAutoExposure() const { return autoExposure_.exposure; }
  float lastAvgLuminance() const { return autoExposure_.avgLuminance; }

  int width() const { return w_; }
  int height() const { return h_; }

  // Effective render scale used for the current frame (scene size / window size).
  float renderScale() const { return renderScaleUsed_; }

  // HDR format currently in use (may differ from settings if the driver lacks
  // support for the requested format).
  HdrBufferFormat hdrBufferFormatActive() const { return hdrBufferFormatActive_; }
  const char* hdrBufferFormatName() const {
    switch (hdrBufferFormatActive_) {
      case HdrBufferFormat::R11G11B10f: return "R11G11B10F";
      case HdrBufferFormat::Rgba16f:
      default: return "RGBA16F";
    }
  }

  // Access the internal FrameGraph (useful for debug UI and custom pipelines).
  FrameGraph& frameGraph() { return frameGraph_; }
  const FrameGraph& frameGraph() const { return frameGraph_; }

private:
  void destroy();
  void createOrResize(int w, int h);

  void drawFullscreen() const;

  int w_{0};
  int h_{0};

  // Dynamic render-scale state.
  // `renderScaleUsed_` is the scale that was applied to allocate the current
  // scene target. `dynamicScale_` is the DRS state that is updated in present()
  // and applied on the next frame via configure().
  float renderScaleUsed_{1.0f};     // scale used for the *current* scene target
  float dynamicScale_{1.0f};        // DRS state (applied next frame)
  float lastUserScale_{1.0f};       // last seen settings.renderScale
  float smoothedFrameMs_{16.67f};   // EMA of CPU frame time
  float drsCooldownSec_{0.0f};      // prevents resize thrash
  bool  drsWasEnabled_{false};

  // HDR render target format selection.
  HdrBufferFormat hdrBufferFormatRequested_{HdrBufferFormat::R11G11B10f};
  HdrBufferFormat hdrBufferFormatActive_{HdrBufferFormat::Rgba16f};
  bool forceRecreate_{false};

  // GL tex allocation params for the active HDR format.
  int hdrInternalFormat_{0};
  unsigned int hdrUploadFormat_{0};
  unsigned int hdrUploadType_{0};

  unsigned int sceneFbo_{0};
  unsigned int sceneTex_{0};
  unsigned int depthRbo_{0};

  // A 1x1 black texture for "null" inputs (e.g. bloom disabled).
  unsigned int blackTex_{0};

  // FrameGraph-managed post-processing resources.
  FrameGraph frameGraph_{};

  unsigned int vao_{0};
  unsigned int vbo_{0};

  ShaderProgram bright_{};
  ShaderProgram blur_{};
  ShaderProgram luminance_{};
  ShaderProgram reduce_{};
  ShaderProgram composite_{};

  // Auto exposure (CPU state; updated from GPU luminance measurement).
  AutoExposureState autoExposure_{};
  // Optional async auto-exposure readback (PBO + fence sync).
  // This avoids stalling the CPU every frame with glGetTexImage when auto-exposure is enabled.
  static constexpr int kAeReadbackFrames = 3;
  std::array<unsigned int, kAeReadbackFrames> aeReadbackPbos_{};
  std::array<void*, kAeReadbackFrames> aeReadbackFences_{};
  std::array<int, kAeReadbackFrames> aeReadbackTags_{};
  int  aeReadbackWrite_{0};
  int  aeReadbackTagCounter_{1};
  bool aeReadbackInited_{false};
  bool aeReadbackSupported_{false};

  float lastPresentTime_{0.0f};
  bool  haveLastPresentTime_{false};

  bool hasCustomComposite_{false};
  std::string compositeFsBuilt_{}; // last successfully built composite FS source
};

} // namespace stellar::render
