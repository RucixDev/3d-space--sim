#pragma once

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/core/Types.h"
#include "stellar/math/Math.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"

#include <string>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Time Trials (Gameplay-friendly but deterministic / sim-side)
// -----------------------------------------------------------------------------
//
// A lightweight "gate course" system meant to support moment-to-moment 3D gameplay
// loops (racing / training) without coupling to rendering or UI.
//
// Gates are represented as a disk (plane + radius) in world space.
// A gate is considered "passed" when the ship crosses the plane from the
// negative side to the positive side of the gate normal AND the intersection
// point lies within the gate radius.
//

struct TimeTrialGate {
  math::Vec3d posKm{};     // gate center in world-space kilometers
  math::Vec3d normal{};    // unit normal; positive side is the "forward" direction
  double radiusKm{2500.0}; // pass radius
};

struct TimeTrialCourse {
  core::u64 key{0};
  core::u64 seed{0};
  std::string name;
  std::vector<TimeTrialGate> gates;
};

struct TimeTrialCourseParams {
  int gateCount{12};
  double gateRadiusKm{2500.0};

  // Course shape around an anchor.
  double courseRadiusKm{50000.0};
  double courseHeightKm{15000.0};
  double jitterKm{8000.0};
  int loops{1};

  // If true, the last gate points toward the first gate.
  bool closedLoop{false};
};

// Generate a "slalom" course around an anchor frame.
// The anchor orientation defines the local axes used for the course:
//   - local +X : right
//   - local +Y : up
//   - local +Z : forward
TimeTrialCourse generateTimeTrialCourseStationSlalomKm(const math::Vec3d& anchorPosKm,
                                                       const math::Quatd& anchorOrient,
                                                       core::u64 seed,
                                                       const TimeTrialCourseParams& p);

// True if the segment [prevPosKm, posKm] crosses the gate plane in the forward
// direction (negative -> positive halfspace) and the intersection lies within
// the gate radius.
//
// velKmS is optional (can be {0,0,0}); when non-zero, it is used as an extra
// robustness check to ensure the motion is generally aligned with the gate normal.
// Optional outTHit is the normalized segment parameter of the plane intersection
// clamped to [0,1] (i.e., 0 = prevPosKm, 1 = posKm). This is useful for
// sub-frame timing (e.g., more accurate lap times when using variable dt).
bool timeTrialGatePassed(const TimeTrialGate& g,
                         const math::Vec3d& prevPosKm,
                         const math::Vec3d& posKm,
                         const math::Vec3d& velKmS,
                         double* outTHit = nullptr);

// Utility: quantize a position to a stable integer hash (helps course keys).
core::u64 hashPosKmQuantized(const math::Vec3d& pKm, double quantumKm);

} // namespace stellar::sim
