#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/ProximityField.h"

#include <cstddef>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// CollisionWarning (headless helpers)
// -----------------------------------------------------------------------------
//
// A lightweight, deterministic helper for cockpit-level collision prediction.
//
// Inputs:
//  - a ProximityFieldKm (typically populated with nearby asteroids)
//  - the ship's current world-space position/velocity
//  - an assumed maximum deceleration (km/s^2)
//
// Output:
//  - earliest predicted impact along the current velocity vector within a horizon
//  - stopping distance/time estimates under the provided decel model
//  - a simple hazard scalar + discrete warning level (advisory/caution/danger)
//
// This is intentionally conservative: it assumes the ship continues on its current
// velocity vector, and it uses a constant deceleration model for the stop estimate.

struct CollisionWarningParams {
  // Max lookahead window for impact prediction.
  double horizonSec{40.0};

  // Inflate obstacles by this amount (ship radius + safety margin).
  double padKm{0.08};

  // Ignore warnings below this speed.
  double minSpeedKmS{0.05};

  // Warning thresholds (seconds to impact).
  double cautionTtiSec{14.0};
  double dangerTtiSec{5.0};

  // If true, also raise severity when the predicted impact distance is within the
  // ship's stopping distance (under maxDecelKmS2).
  bool useStopDistance{true};

  // Extra margin applied to the stopping distance.
  // stopDistWithMargin = stopDist * (1 + stopMarginFactor)
  double stopMarginFactor{0.10};
};

enum class CollisionWarningLevel : int {
  None = 0,
  Advisory = 1,
  Caution = 2,
  Danger = 3,
};

struct CollisionWarningResult {
  bool hasImpact{false};

  // ID returned by ProximityFieldKm (obstacle index).
  int obstacleId{-1};

  // Time-to-impact along linear motion (seconds).
  double ttiSec{0.0};

  // Distance along the velocity ray to the first hit point (km).
  double impactDistKm{0.0};

  // Predicted hit point under constant velocity.
  math::Vec3d impactPointKm{0, 0, 0};

  // Snapshot values used to compute hazard.
  double speedKmS{0.0};
  double maxDecelKmS2{0.0};

  double stopTimeSec{0.0};
  double stopDistKm{0.0};
  double stopDistWithMarginKm{0.0};

  // Positive => can stop before the inflated obstacle along the current velocity.
  double marginKm{0.0};
  bool canStopBeforeImpact{true};

  // 0..1, where 1 indicates imminent impact (<= danger threshold).
  double hazard01{0.0};
  CollisionWarningLevel level{CollisionWarningLevel::None};
};

// Constant-deceleration estimates.
double stoppingTimeSec(double speedKmS, double decelKmS2);
double stoppingDistanceKm(double speedKmS, double decelKmS2);

// Compute a collision warning result for the earliest impact along the velocity
// vector within params.horizonSec.
CollisionWarningResult computeCollisionWarning(const ProximityFieldKm& field,
                                              const math::Vec3d& posKm,
                                              const math::Vec3d& velKmS,
                                              double maxDecelKmS2,
                                              const CollisionWarningParams& params = {},
                                              int ignoreId = -1);

} // namespace stellar::sim
