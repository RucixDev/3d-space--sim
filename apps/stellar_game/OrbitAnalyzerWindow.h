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
  // Optional shared reference body selection used by the trajectory preview and
  // maneuver node RTN frame.
  // -1 = auto/dominant, 0 = star, 1..N = planet index+1
  int* refBodyChoice{};

  // Optional maneuver node bindings (seconds-from-now + RTN delta-v in m/s).
  // If any are missing, planner actions are disabled.
  bool* maneuverNodeEnabled{};
  float* maneuverNodeTimeSec{};
  float* dvAlongMS{};
  float* dvNormalMS{};
  float* dvRadialMS{};
};

void drawOrbitAnalyzerWindow(OrbitAnalyzerWindowState& state,
                            const OrbitAnalyzerBindings& bindings,
                            const stellar::sim::StarSystem& sys,
                            double timeDays,
                            const stellar::sim::Ship& ship,
                            const stellar::sim::GravityParams& gravityParams);

} // namespace stellar::game
