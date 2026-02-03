#include "stellar/sim/MissileDefense.h"

#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace stellar::sim {

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

    const math::Vec3d toDir = math::safeNormalized(toTarget, math::Vec3d{0, 0, 1}, 1e-12);
    const math::Vec3d mvDir = math::safeNormalized(m.velKmS, math::Vec3d{0, 0, 1}, 1e-12);
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
  return math::safeNormalized(b * cs + c * sn, b, 1e-12);
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
  const math::Vec3d los = math::safeNormalized(r, math::Vec3d{0, 0, 1}, 1e-12);

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

  // Optional: radar "beaming" bias.
  //
  // If the inbound missile is a radar seeker that models a doppler notch, bias the
  // jink direction toward rotating the target's velocity into the plane perpendicular
  // to the current line-of-sight. This increases LOS rate / reduces range-rate in a
  // deterministic way and synergizes with the notch mechanic.
  if (params.enableRadarBeaming && missile.seeker == MissileSeekerType::Radar) {
    const double notch = std::max(0.0, missile.radarDopplerNotchKmS);
    const double blend = std::clamp(params.radarBeamBlend, 0.0, 1.0);
    const double engageMul = std::max(0.0, params.radarBeamEngageNotchMultiple);
    if (notch > 0.0 && blend > 0.0) {
      // LOS from missile to target.
      const math::Vec3d toDir = -los;
      const double vrKmS = math::dot(targetVelKmS - missile.velKmS, toDir);
      const double absVr = std::fabs(vrKmS);

      if (absVr > notch * engageMul) {
        // Desired velocity direction: current velocity projected into the LOS-perpendicular plane.
        const double vT = targetVelKmS.length();

        math::Vec3d vPerp = targetVelKmS - toDir * math::dot(targetVelKmS, toDir);
        if (vPerp.lengthSq() < 1e-12) {
          // Degenerate (already LOS-aligned or nearly zero speed): pick any perpendicular direction.
          const math::Vec3d u = unitPerpFromSeed(toDir, seed ^ 0x6d4f2e31b7a3c0d5ull);
          vPerp = (vT > 1e-9) ? (u * vT) : u;
        } else {
          vPerp = vPerp.normalized() * std::max(0.0, vT);
        }

        const math::Vec3d deltaV = vPerp - targetVelKmS;
        const math::Vec3d fallback = unitPerpFromSeed(toDir, seed ^ 0x9e3779b97f4a7c15ull);
        const math::Vec3d beamDir = math::safeNormalized(deltaV, fallback, 1e-12);

        // Engagement weight ramps up as we exceed the notch threshold.
        const double w = std::clamp((absVr - notch) / std::max(absVr, 1e-9), 0.0, 1.0) * blend;

        dir = math::safeNormalized(dir * (1.0 - w) + beamDir * w, dir, 1e-12);
      }
    }
  }

  if (params.enforceLateralToLos) {
    // Project into the plane perpendicular to LOS (lateral jink).
    dir = dir - los * math::dot(dir, los);
    dir = math::safeNormalized(dir, unitPerpFromSeed(los, seed ^ 0x9e3779b97f4a7c15ull), 1e-12);
  } else {
    dir = math::safeNormalized(dir, unitPerpFromSeed(los, seed ^ 0x9e3779b97f4a7c15ull), 1e-12);
  }

  out.dirWorld = dir;
  out.valid = (dir.lengthSq() > 1e-12);
  return out;
}

}  // namespace stellar::sim
