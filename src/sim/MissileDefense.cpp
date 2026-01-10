#include "stellar/sim/MissileDefense.h"

#include <algorithm>

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

}  // namespace stellar::sim