#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/Combat.h"

#include <cstddef>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// MissileDefense (headless helpers)
// -----------------------------------------------------------------------------
//
// Utilities for detecting inbound missiles in a deterministic, frame-rate
// agnostic way. This is used by the game app to give NPCs basic defensive
// behaviors (countermeasures + evasive thrust jinks) without baking the logic
// directly into the big game loop.

struct MissileThreatParams {
  // Cosine of the maximum "off boresight" angle for considering a missile
  // inbound. 1.0 means perfectly aligned, 0.0 means 90 degrees.
  double minApproachCos{0.2};

  // Ignore cases where relative closing speed is extremely low.
  double minClosingKmS{0.02};

  // Ignore missiles farther than this distance.
  double maxConsiderDistKm{250000.0};
};

struct MissileThreatSummary {
  bool inbound{false};
  MissileSeekerType seeker{MissileSeekerType::Heat};

  double distKm{0.0};
  double closingKmS{0.0};
  double ttiSec{0.0};
  double approachCos{0.0};

  std::size_t missileIndex{0};
  core::u64 shooterId{0};
  bool fromPlayer{false};
};

// Recommended evasion direction for an inbound missile (deterministic).
//
// The goal is to provide a *direction only*; the game can decide how strongly
// to apply thrust based on pilot skill, ship performance, and time-to-impact.
struct MissileEvasionPlan {
  bool valid{false};

  // Unit vector (world space) pointing in a good direction for applying
  // lateral thrust to increase the predicted miss distance.
  math::Vec3d dirWorld{0, 0, 0};

  // Time of closest approach under constant relative velocity (seconds).
  double tClosestSec{0.0};

  // Predicted miss distance at closest approach if the target does nothing (km).
  double missDistanceKm{0.0};

  // Relative closing speed along the line of sight (km/s). Positive means inbound.
  double closingKmS{0.0};
};

struct MissileEvasionParams {
  // Ignore cases where relative speed is extremely low (fallback direction is used).
  double minRelSpeedKmS{0.02};

  // Treat the closest-approach offset as degenerate below this threshold (km).
  double minMissVecKm{1e-3};

  // If true, project the output onto the plane perpendicular to the current
  // line-of-sight (produces a more "lateral jink" that tends to increase LOS rate).
  bool enforceLateralToLos{true};
};

// Compute a suggested evasion direction against a specific missile/target pair.
//
// seed:
//   Any stable seed (hash(universeSeed, npcId, timePhase, etc.)). Only used to break
//   ties in degenerate head-on cases.
MissileEvasionPlan planMissileEvasion(const Missile& missile,
                                     const math::Vec3d& targetPosKm,
                                     const math::Vec3d& targetVelKmS,
                                     core::u64 seed,
                                     const MissileEvasionParams& params = {});

// Find the nearest (minimum time-to-impact) inbound missile tracking the given
// target, according to a simple relative kinematic estimate.
MissileThreatSummary nearestInboundMissile(const Missile* missiles,
                                          std::size_t missileCount,
                                          CombatTargetKind targetKind,
                                          core::u64 targetId,
                                          const math::Vec3d& targetPosKm,
                                          const math::Vec3d& targetVelKmS,
                                          const MissileThreatParams& params = {});

}  // namespace stellar::sim