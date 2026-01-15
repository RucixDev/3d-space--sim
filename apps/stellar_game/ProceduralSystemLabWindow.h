#pragma once

#include "stellar/core/Types.h"

#include <functional>
#include <string>
#include <vector>
#include <cstdint>

namespace stellar {
namespace sim {
  class Universe;
  using SystemId = core::u64;
}

namespace game {

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

struct ProceduralSystemLabWindowState {
  bool open{false};
  bool followCurrentSystem{true};
  sim::SystemId systemId{0};

  // UI
  bool showPlanets{true};
  bool showMoons{true};
  bool showZones{true};

  // --- Procedural ring preview (render-only debug) ---
  bool showRings{false};
  double ringChanceMul{1.0};
  double ringOpacity{1.0};
  int ringPlanetIndex{0};
  bool ringUseAutoVariant{true};
  int ringVariant{0}; // 0..2
  int ringPreviewW{256};
  int ringPreviewH{80};
  bool ringPreviewAlphaOnly{false};

  // Cached preview image (CPU-side).
  core::u64 ringPreviewSeed{0};
  int ringPreviewWCache{0};
  int ringPreviewHCache{0};
  std::vector<std::uint8_t> ringPreviewRGBA; // w*h*4
  std::vector<float> ringRadialMean01;       // h values in [0..1]

  // --- Signals / mining sites ---
  bool showSignals{false};
  int resourceFieldCount{3};
  bool includeDailyDerelict{false};
  bool includeDistress{false};
  bool includeTrafficConvoys{false};

  // Time selection for deterministic daily signals.
  bool useCurrentTime{true};
  double timeDaysOverride{0.0};

  // Resource field selection / preview.
  int selectedResourceField{-1};
  bool showAsteroidScatter{true};
  int scatterMaxPoints{512};

  // Scatter visualization controls.
  bool scatterDensityHeatmap{true};
  int scatterHeatmapRes{48};
  bool scatterColorByDensity{true};
};

void drawProceduralSystemLabWindow(ProceduralSystemLabWindowState& state,
                                  sim::Universe& universe,
                                  sim::SystemId currentSystemId,
                                  double timeDays,
                                  float timeSec,
                                  const ToastFn& toast);

} // namespace game
} // namespace stellar
