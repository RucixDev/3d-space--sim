#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/RenderTarget.h"
#include "stellar/render/ShaderToyGraph.h"
#include "stellar/render/ShaderToyParams.h"

#include <array>
#include <functional>
#include <string>

namespace stellar::game {

struct ProceduralShaderLabPassState {
  // Buffers can be disabled; Image is always rendered.
  bool enabled{true};

  // Resolution divisor for this pass (1 = full, 2 = half, 4 = quarter).
  int resolutionScale{1};

  // Channel routing (iChannel0..3).
  std::array<stellar::render::ShaderToyChannelSource, 4> channels{
    stellar::render::ShaderToyChannelSource::None,
    stellar::render::ShaderToyChannelSource::None,
    stellar::render::ShaderToyChannelSource::None,
    stellar::render::ShaderToyChannelSource::None,
  };

  bool dirty{false};
};

struct ProceduralShaderLabWindowState {
  bool open{false};

  // Preview
  bool previewInited{false};
  int previewResolution{512};
  stellar::render::RenderTarget2D previewTarget{};
  std::string previewInitError{};

  // Multi-pass graph (A..D buffers + Image).
  bool graphInited{false};
  stellar::render::ShaderToyGraph graph{};

  bool showGeneratedSource{false};

  // Editor
  int presetIndex{0};
  bool autoCompileOnPreset{true};

  // 0..3 = Buffer A..D, 4 = Image
  int editPassIndex{4};
  int previewPassIndex{4};

  // Per-pass settings
  static constexpr int kPassCount = 5;
  std::array<ProceduralShaderLabPassState, kPassCount> passes{};

  // NOTE: Fixed-size buffers keep ImGui usage simple.
  static constexpr int kCodeCapacity = 32768;
  std::array<std::array<char, kCodeCapacity>, kPassCount> code{};

  std::string fileError{};

  // File I/O
  char filePath[256]{"shaders/proc_shader_graph.stoy"};
  bool liveReload{false};
  bool liveReloadMissingOk{true};
  std::uint64_t liveReloadStamp{0};

  // Time
  bool useRealTime{true};
  float timeOverride{0.0f};
  float timeScale{1.0f};

  bool paused{false};
  bool stepFrame{false};
  float fixedStepSec{1.0f / 60.0f};

  // Internal time accumulator when useRealTime=true.
  bool timeInit{false};
  float timeAccumSec{0.0f};
  float lastRealTimeSec{0.0f};

  // iFrame counter (incremented when rendering advances).
  int frame{0};

  // Camera (orbit)
  float yawDeg{35.0f};
  float pitchDeg{15.0f};
  float distance{3.0f};
  bool autoOrbit{true};
  float orbitDegPerSec{8.0f};
  float fovYDeg{60.0f};

  // Inputs
  core::u64 seed{0xC0FFEEULL};

  // Parsed user parameters (`// @param ...`) shared across the graph.
  // Values are live-updated by the UI and sent as uniforms each frame.
  stellar::render::ShaderToyParamSet params{};

  // UI
  bool paramsOpen{true};
  bool paramsShowAdvanced{false};
  char paramsFilter[64]{""};

  // Mouse stored normalized in [0,1] relative to the preview image.
  float mouseU{0.0f};
  float mouseV{0.0f};
  float mouseDownU{0.0f};
  float mouseDownV{0.0f};
  bool mouseDown{false};

  // One-shot actions (handled inside the render section to keep GL state tidy).
  bool requestResetBuffers{false};
  bool requestCompileAll{false};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawProceduralShaderLabWindow(ProceduralShaderLabWindowState& st, float timeSec, const ToastFn& toast);

} // namespace stellar::game
