#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/Gravity.h"
#include "stellar/sim/System.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Trajectory predictor (headless)
// -----------------------------------------------------------------------------
//
// A lightweight path propagator for in-system debugging and gameplay.
//
// - Uses RK4 integration so the preview stays stable even with larger step sizes.
// - Supports a single instantaneous maneuver node (delta-v) at a requested time.
// - Deterministic + headless; safe to call from tools/tests.
//
// Units:
//   - position: km
//   - velocity: km/s
//   - time: seconds (relative to the start), plus startTimeDays for planet ephemerides.

struct ManeuverNode {
  // Time after the start of the prediction when the burn is applied.
  double timeSec{0.0};

  // Instantaneous delta-v applied in world space.
  math::Vec3d deltaVKmS{0, 0, 0};
};

struct TrajectorySample {
  double tSec{0.0};
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
};

struct TrajectoryPredictParams {
  // Total prediction time horizon (seconds).
  double horizonSec{1800.0};

  // Integrator step size (seconds). The final step is shortened to end exactly at horizonSec.
  double stepSec{2.0};

  // Safety cap for output size.
  int maxSamples{2000};

  // Include gravity from the provided system (star + planets) using GravityParams.
  bool includeGravity{true};
  GravityParams gravity{};
};

// Predict a future trajectory using 4th-order Runge-Kutta integration.
//
// If `node` is provided and node->timeSec is within the prediction horizon, an instantaneous
// delta-v is applied at that exact time (the integrator step is split as needed).
std::vector<TrajectorySample> predictTrajectoryRK4(const StarSystem& sys,
                                                   double startTimeDays,
                                                   const math::Vec3d& startPosKm,
                                                   const math::Vec3d& startVelKmS,
                                                   const TrajectoryPredictParams& params,
                                                   const ManeuverNode* node = nullptr);

} // namespace stellar::sim
