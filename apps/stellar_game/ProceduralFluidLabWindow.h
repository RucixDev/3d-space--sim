#pragma once

#include "stellar/core/Types.h"
#include "stellar/proc/FluidSim2D.h"
#include "stellar/render/Texture.h"

#include <functional>
#include <string>
#include <vector>

namespace stellar::game {

struct ProceduralFluidLabWindowState {
  bool open{false};

  // Simulation
  int gridSize{160};
  int iterations{25};
  bool paused{false};
  bool singleStep{false};
  bool autoDt{true};
  float fixedDt{1.0f / 60.0f};
  float maxDt{1.0f / 20.0f};

  // Visuals
  bool showVelocity{false};
  float displayExposure{0.08f};
  float velocityVizScale{0.015f};
  int previewPixels{520};

  // Brush
  float brushRadius01{0.04f};
  float brushDye{8.0f};
  float brushForce{40.0f};
  float brushColor[3]{0.85f, 0.55f, 0.15f};
  bool rightClickErases{true};

  // Procedural forcing seed
  core::u64 noiseSeed{0xBADC0FFEEull};

  // Export
  char exportPath[256]{"screenshots/fluid_lab.png"};
  bool exportFlipY{true};

  // Internal
  proc::FluidSim2D sim{};
  render::Texture2D tex{};
  std::vector<unsigned char> rgba{};
  float lastTimeSec{0.0f};
  bool mouseDown{false};
  float lastMouseU{0.0f};
  float lastMouseV{0.0f};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawProceduralFluidLabWindow(ProceduralFluidLabWindowState& st, float timeSec, const ToastFn& toast);

} // namespace stellar::game
