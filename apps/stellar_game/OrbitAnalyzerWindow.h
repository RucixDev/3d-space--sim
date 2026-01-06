#pragma once

#include "stellar/sim/Gravity.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/System.h"

namespace stellar::game {

struct OrbitAnalyzerWindowState {
  bool open{false};

  // When enabled, multiply mu by gravityParams.scale so readouts match the game's
  // "effective" gravity strength.
  bool useGravityScale{true};

  // Planner targets (altitudes above the reference body's radius).
  float targetApoAltKm{200000.0f};
  float targetPeriAltKm{200000.0f};
};

struct OrbitAnalyzerBindings {
  // Optional: share the trajectory preview reference body selection.
  // -1 = auto (dominant), 0 = star, 1..n = planet index + 1
  int* refBodyChoice{nullptr};

  // Maneuver node controls to write into the existing planner UI.
  bool* maneuverNodeEnabled{nullptr};
  float* maneuverNodeTimeSec{nullptr};
  float* dvAlongMS{nullptr};
  float* dvNormalMS{nullptr};
  float* dvRadialMS{nullptr};
};

void drawOrbitAnalyzerWindow(OrbitAnalyzerWindowState& state,
                            const stellar::sim::StarSystem& sys,
                            double timeDays,
                            const stellar::sim::Ship& ship,
                            const stellar::sim::GravityParams& gravityParams,
                            OrbitAnalyzerBindings bindings);

} // namespace stellar::game
