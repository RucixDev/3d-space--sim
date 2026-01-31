#pragma once

#include "stellar/sim/Gravity.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/System.h"
#include "stellar/sim/TrajectoryAnalysis.h"
#include "stellar/sim/TrajectoryEvents.h"

#include <string>
#include <vector>

namespace stellar::game {

struct OrbitAnalyzerWindowState {
  bool open{false};

  // When enabled, multiply mu by gravityParams.scale so readouts match the game's
  // effective gravity strength.
  bool useGravityScale{true};

  // Planner targets (altitudes above the reference body's radius).
  float targetApoAltKm{200000.0f};
  float targetPeriAltKm{200000.0f};

  // Plane-change planner options.
  bool planeAlignForcePrograde{false};

  // Forecast controls.
  bool forecastAutoUpdate{false};
  float forecastAutoUpdateIntervalSec{1.0f};
  float forecastHorizonMin{60.0f};
  float forecastStepSec{10.0f};
  int forecastMaxSamples{900};
  int forecastRefineDepth{14};

  // Forecast cache.
  double forecastLastComputeTimeDays{-1.0};
  bool forecastValid{false};
  std::string forecastStatus{};
  std::vector<stellar::sim::TrajectorySample> forecastSamples{};
  stellar::sim::TrajectoryAnalysisResult forecastAnalysis{};
  std::vector<stellar::sim::DominantBodyTransition> forecastDominantTransitions{};
};

struct OrbitAnalyzerBindings {
  float* maneuverNodeTimeSec{};
  float* maneuverNodeDirX{};
  float* maneuverNodeDirY{};
  float* maneuverNodeDirZ{};

  float* maneuverNodeDvMs{};

  bool* maneuverNodeEnabled{};
  int* maneuverNodeRefBodyChoice{};
};

void drawOrbitAnalyzerWindow(OrbitAnalyzerWindowState& state,
                            const OrbitAnalyzerBindings& bindings,
                            const stellar::sim::StarSystem& sys,
                            double timeDays,
                            const stellar::sim::Ship& ship,
                            const stellar::sim::GravityParams& gravityParams);

} // namespace stellar::game
