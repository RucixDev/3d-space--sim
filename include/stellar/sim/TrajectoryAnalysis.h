#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/Gravity.h"
#include "stellar/sim/TrajectoryPredictor.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// TrajectoryAnalysis (headless)
// -----------------------------------------------------------------------------
//
// Utility routines to extract actionable information from a discrete trajectory
// (usually emitted by TrajectoryPredictor).
//
// Primary use cases:
//  - Tools/UI: closest approach readouts, collision warnings, encounter timing.
//  - Gameplay: preflight checks (impact within horizon), mission logic triggers.
//
// Design constraints:
//  - Deterministic + headless.
//  - Cheap: O(N samples * M bodies), where M is star + planets.
//
// Units:
//  - position: km
//  - velocity: km/s
//  - time: seconds (relative to the prediction start)

struct TrajectoryApproach {
  bool valid{false};
  GravityBody body{};

  // Time of closest approach (seconds after prediction start).
  double tSec{0.0};

  // Distance to body center at closest approach (km).
  double distanceKm{0.0};

  // Altitude above body surface (km). May be negative if inside the radius.
  double altitudeKm{0.0};

  // Relative state at closest approach (body-centric).
  math::Vec3d relPosKm{0, 0, 0};
  math::Vec3d relVelKmS{0, 0, 0};
};

struct TrajectoryImpact {
  bool valid{false};
  GravityBody body{};

  // Time of first surface intersection (seconds after prediction start).
  double tSec{0.0};

  // Approximate ship state at impact time (world/inertial frame).
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
};

struct TrajectoryAnalysisResult {
  // Closest approach for every considered body, in gather order.
  std::vector<TrajectoryApproach> closestApproachByBody{};

  // Smallest distance/altitude across all bodies.
  TrajectoryApproach closestOverall{};

  // Earliest detected impact with any body (if any).
  TrajectoryImpact firstImpact{};
};

// Analyze a discrete trajectory against system gravity bodies.
//
// Notes:
//  - Body motion is accounted for by evaluating ephemerides at the sample times.
//  - Closest approaches are computed using segment-min distance (not just sample
//    points), improving accuracy for coarse sampling.
TrajectoryAnalysisResult analyzeTrajectory(const StarSystem& sys,
                                          double startTimeDays,
                                          const std::vector<TrajectorySample>& samples,
                                          const GravityParams& params = {});

} // namespace stellar::sim
