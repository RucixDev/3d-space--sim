#pragma once

#include "stellar/math/Vec3.h"

#include <optional>

namespace stellar::sim {

// Utilities for simple projectile lead / intercept math.
//
// We model a projectile that travels at constant speed `projectileSpeedKmS` in world space,
// starting at the shooter's current position and inheriting the shooter's velocity.
// This matches the game's projectile spawning rule:
//   projectileVel = shooterVel + dir * projectileSpeed.
//
// The intercept time solves:
//   |(targetPos - shooterPos) + (targetVel - shooterVel) * t| = projectileSpeed * t
// for t > 0.

// Returns the smallest positive intercept time (>= minTimeSec) if a valid solution exists.
std::optional<double> solveInterceptTimeSec(
  const math::Vec3d& shooterPosKm,
  const math::Vec3d& shooterVelKmS,
  const math::Vec3d& targetPosKm,
  const math::Vec3d& targetVelKmS,
  double projectileSpeedKmS,
  double minTimeSec = 1.0e-5
);

struct LeadSolution {
  double tSec{0.0};
  math::Vec3d leadPointKm{0, 0, 0};
  math::Vec3d aimDirWorld{0, 0, 1};
};

// Convenience wrapper that returns intercept time + lead point + unit aim direction.
// If maxTimeSec is finite, the solution must satisfy tSec <= maxTimeSec.
std::optional<LeadSolution> solveProjectileLead(
  const math::Vec3d& shooterPosKm,
  const math::Vec3d& shooterVelKmS,
  const math::Vec3d& targetPosKm,
  const math::Vec3d& targetVelKmS,
  double projectileSpeedKmS,
  double maxTimeSec = 1.0e9,
  double minTimeSec = 1.0e-5
);

} // namespace stellar::sim
