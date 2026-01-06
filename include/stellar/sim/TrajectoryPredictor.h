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
// - Deterministic + headless; safe to call from tools/tests.
// - Supports a single instantaneous maneuver node (delta-v) at a requested time.
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

  // Integrator step size (seconds).
  //
  // For predictTrajectoryRK4(): this is the fixed integrator step.
  // For predictTrajectoryRK45Adaptive(): this is the *output sampling* interval.
  //   (internal adaptive steps are chosen automatically).
  double stepSec{2.0};

  // Safety cap for output size.
  int maxSamples{2000};

  // Include gravity from the provided system (star + planets) using GravityParams.
  bool includeGravity{true};
  GravityParams gravity{};

  // ---------------------------------------------------------------------------
  // Adaptive stepping controls (RK45)
  // ---------------------------------------------------------------------------
  //
  // Used only by predictTrajectoryRK45Adaptive(). These are intentionally opt-in
  // so existing gameplay behavior stays stable.

  // Smallest / largest internal step (seconds).
  double minStepSec{0.05};
  double maxStepSec{30.0};

  // Error tolerances for adaptive step control.
  //
  // A step is accepted if the estimated local truncation error is <= 1.0 when
  // normalized by the tolerance scales.
  double relTol{1e-6};
  double absTolPosKm{1e-3};   // 1 meter
  double absTolVelKmS{1e-6};  // 1 mm/s
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

// Predict a future trajectory using an adaptive Dormand–Prince RK5(4) (RK45) integrator.
//
// Differences vs RK4:
// - Uses an embedded 5th/4th order pair to estimate local error and adapt the internal timestep.
// - `params.stepSec` controls output sampling interval (trajectory samples are emitted at that
//   cadence, plus a sample at the maneuver node time if provided).
// - Internal integration steps are clamped to [minStepSec, maxStepSec].
//
// This is useful for long-horizon previews: you can request coarse output samples (stepSec)
// while still keeping accuracy during fast dynamics (close passes) via smaller internal steps.
std::vector<TrajectorySample> predictTrajectoryRK45Adaptive(const StarSystem& sys,
                                                           double startTimeDays,
                                                           const math::Vec3d& startPosKm,
                                                           const math::Vec3d& startVelKmS,
                                                           const TrajectoryPredictParams& params,
                                                           const ManeuverNode* node = nullptr);

} // namespace stellar::sim
