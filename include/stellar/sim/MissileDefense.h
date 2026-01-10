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
