#pragma once

#include "stellar/sim/LambertPlanner.h"

#include <imgui.h>

#include <string>
#include <vector>

namespace stellar::game {

// Small ImGui helper to visualize a LambertPorkchopResult as a classic
// “porkchop” heatmap (departure vs time-of-flight).

enum class PorkchopMetric : int {
  // Lower is better for all metrics.
  Score = 0,
  DepartDv = 1,
  ArrivalRelSpeed = 2,
  TotalDv = 3,
};

struct PorkchopPlotOptions {
  PorkchopMetric metric{PorkchopMetric::Score};

  // Gamma curve applied to the normalized value (1 = linear). Values < 1 make
  // small differences near “best” more visible.
  float gamma{0.65f};

  // Downsample caps. The plot will bin the full grid into <= maxCols x maxRows.
  int maxCols{220};
  int maxRows{160};

  bool showBestMarkers{true};
  bool showAxesLabels{true};
};

struct PorkchopPlotPick {
  bool clicked{false};
  double departAfterSec{0.0};
  double tofSec{0.0};
};

// Draw a heatmap from a completed porkchop search result.
//
// Returns a pick when the user clicks a cell in the plot.
PorkchopPlotPick drawPorkchopPlot(const sim::LambertPorkchopResult& res,
                                  const std::vector<sim::LambertTransferMetrics>& best,
                                  const PorkchopPlotOptions& opt,
                                  ImVec2 size = ImVec2(0.0f, 0.0f));

// Export the grid to a CSV file (one row per cell).
bool exportPorkchopCsv(const char* path,
                       const sim::LambertPorkchopResult& res,
                       std::string* err = nullptr);

} // namespace stellar::game
