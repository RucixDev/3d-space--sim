#include "stellar/sim/MissileDefense.h"

#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace stellar::sim {

static math::Vec3d safeNormalized(const math::Vec3d& v, const math::Vec3d& fallback) {
  if (v.lengthSq() < 1e-12) return fallback;
  return v.normalized();
}

MissileThreatSummary nearestInboundMissile(const Missile* missiles,
                                          std::size_t missileCount,
                                          CombatTargetKind targetKind,
                                          core::u64 targetId,
                                          const math::Vec3d& targetPosKm,
                                          const math::Vec3d& targetVelKmS,
                                          const MissileThreatParams& params) {
  MissileThreatSummary best{};
  double bestTtiSec = std::numeric_limits<double>::infinity();

  if (!missiles || missileCount == 0) return best;

  const double minCos = std::clamp(params.minApproachCos, -1.0, 1.0);
  const double minClosing = std::max(0.0, params.minClosingKmS);
  const double maxDist = std::max(0.0, params.maxConsiderDistKm);

  for (std::size_t i = 0; i < missileCount; ++i) {
    const Missile& m = missiles[i];

    if (!m.hasTarget) continue;
    if (m.targetKind != targetKind) continue;
    if (m.targetId != targetId) continue;
    if (m.ttlSimSec <= 0.0) continue;

    const math::Vec3d toTarget = targetPosKm - m.posKm;
    const double distKm = toTarget.length();
    if (distKm > maxDist) continue;

    const math::Vec3d toDir = safeNormalized(toTarget, math::Vec3d{0, 0, 1});
    const math::Vec3d mvDir = safeNormalized(m.velKmS, math::Vec3d{0, 0, 1});
    const double approachCos = math::dot(mvDir, toDir);
    if (approachCos < minCos) continue;

    // Relative closing speed along the line-of-sight.
    const math::Vec3d relVel = m.velKmS - targetVelKmS;
    const double closingKmS = math::dot(relVel, toDir);
    if (closingKmS <= minClosing) continue;

    const double ttiSec = (distKm > 1e-9) ? (distKm / std::max(closingKmS, 1e-9)) : 0.0;
    if (ttiSec < bestTtiSec) {
      bestTtiSec = ttiSec;

      best.inbound = true;
      best.seeker = m.seeker;
      best.distKm = distKm;
      best.closingKmS = closingKmS;
      best.ttiSec = ttiSec;
      best.approachCos = approachCos;
      best.missileIndex = i;
      best.shooterId = m.shooterId;
      best.fromPlayer = m.fromPlayer;
    }
  }

  return best;
}

static math::Vec3d unitPerpFromSeed(const math::Vec3d& axisIn, core::u64 seed) {
  math::Vec3d axis = axisIn;
  if (axis.lengthSq() < 1e-12) axis = {0, 0, 1};
  axis = axis.normalized();

  math::Vec3d b = math::cross(axis, math::Vec3d{0, 1, 0});
  if (b.lengthSq() < 1e-12) b = math::cross(axis, math::Vec3d{1, 0, 0});
  if (b.lengthSq() < 1e-12) b = math::cross(axis, math::Vec3d{0, 0, 1});
  if (b.lengthSq() < 1e-12) return {1, 0, 0};
  b = b.normalized();

  math::Vec3d c = math::cross(axis, b);
  if (c.lengthSq() < 1e-12) return b;
  c = c.normalized();

  core::SplitMix64 rng(seed);
  const double ang = rng.range(0.0, 2.0 * std::numbers::pi);
  const double cs = std::cos(ang);
  const double sn = std::sin(ang);
  return safeNormalized(b * cs + c * sn, b);
}

MissileEvasionPlan planMissileEvasion(const Missile& missile,
                                     const math::Vec3d& targetPosKm,
                                     const math::Vec3d& targetVelKmS,
                                     core::u64 seed,
                                     const MissileEvasionParams& params) {
  MissileEvasionPlan out{};

  // Relative motion: missile in the target frame.
  const math::Vec3d r = missile.posKm - targetPosKm;
  const math::Vec3d v = missile.velKmS - targetVelKmS;

  const double distSq = r.lengthSq();
  const double dist = std::sqrt(std::max(0.0, distSq));
  const math::Vec3d los = safeNormalized(r, math::Vec3d{0, 0, 1});

  const double vSq = v.lengthSq();
  const double vLen = std::sqrt(std::max(0.0, vSq));

  // Closing speed along LOS (positive when inbound).
  const double closing = -math::dot(v, los);
  out.closingKmS = closing;

  const double minRel = std::max(0.0, params.minRelSpeedKmS);
  if (!(vSq > 1e-18) || vLen < minRel || dist < 1e-9) {
    // Degenerate: pick any direction perpendicular to LOS.
    out.dirWorld = unitPerpFromSeed(los, seed);
    out.valid = (out.dirWorld.lengthSq() > 1e-12);
    return out;
  }

  if (closing <= 0.0) {
    // Not inbound (or grazing away).
    return out;
  }

  // Time of closest approach under constant relative velocity (clamped to now..inf).
  double tClosest = -math::dot(r, v) / vSq;
  if (!std::isfinite(tClosest)) tClosest = 0.0;
  if (tClosest < 0.0) tClosest = 0.0;
  out.tClosestSec = tClosest;

  // Predicted relative offset at closest approach.
  math::Vec3d missVec = r + v * tClosest;
  // Numerically ensure perpendicular to v (especially if tClosest was clamped).
  missVec = missVec - v * (math::dot(missVec, v) / vSq);

  const double missDist = missVec.length();
  out.missDistanceKm = missDist;

  math::Vec3d dir{0, 0, 0};
  if (missDist > std::max(0.0, params.minMissVecKm)) {
    // Move away from the predicted closest-approach point.
    dir = (-missVec) / missDist;
  } else {
    // Head-on: any direction perpendicular to relative velocity works; break ties via seed.
    dir = unitPerpFromSeed(v, seed);
  }

  if (params.enforceLateralToLos) {
    // Project into the plane perpendicular to LOS (lateral jink).
    dir = dir - los * math::dot(dir, los);
    dir = safeNormalized(dir, unitPerpFromSeed(los, seed ^ 0x9e3779b97f4a7c15ull));
  } else {
    dir = safeNormalized(dir, unitPerpFromSeed(los, seed ^ 0x9e3779b97f4a7c15ull));
  }

  out.dirWorld = dir;
  out.valid = (dir.lengthSq() > 1e-12);
  return out;
}

}  // namespace stellar::sim
