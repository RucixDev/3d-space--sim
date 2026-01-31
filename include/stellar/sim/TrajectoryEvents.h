#pragma once

#include "stellar/sim/Gravity.h"
#include "stellar/sim/TrajectoryPredictor.h"

#include <vector>

namespace stellar::sim {

struct DominantBodyTransition {
  // Seconds since prediction start.
  double tSec{0.0};

  GravityBody from{};
  GravityBody to{};
};

struct DominantBodyTransitionParams {
  // Maximum recursion depth used to refine transition times.
  int refineDepth{14};

  // Minimum separation between adjacent reported transitions (seconds).
  double minSeparationSec{0.0};
};

// Detect times when the dominant gravity body changes along a discrete trajectory.
//
// Dominance is defined as the gravity body that contributes the largest |a|
// at the query point (see dominantGravityBody()).
std::vector<DominantBodyTransition> detectDominantBodyTransitions(const StarSystem& sys,
                                                                 double startTimeDays,
                                                                 const std::vector<TrajectorySample>& samples,
                                                                 const GravityParams& gravityParams = {},
                                                                 DominantBodyTransitionParams params = {});

} // namespace stellar::sim
