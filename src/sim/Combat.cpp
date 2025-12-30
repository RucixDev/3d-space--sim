#include "stellar/sim/Combat.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

bool raySphereIntersectKm(const math::Vec3d& originKm,
                          const math::Vec3d& dirNormalized,
                          const math::Vec3d& centerKm,
                          double radiusKm,
                          double& outTEnterKm) {
  radiusKm = std::max(0.0, radiusKm);
  const math::Vec3d oc = centerKm - originKm;
  const double tProj = math::dot(oc, dirNormalized);
  const double dist2 = oc.lengthSq();
  const double d2 = dist2 - tProj * tProj;
  const double r2 = radiusKm * radiusKm;
  if (d2 > r2) return false;

  const double thc = std::sqrt(std::max(0.0, r2 - d2));
  double tEnter = tProj - thc;
  // If we start inside the sphere, tEnter can be negative; clamp to 0.
  if (tEnter < 0.0) tEnter = 0.0;
  outTEnterKm = tEnter;
  return true;
}

RaycastHit raycastNearestSphereKm(const math::Vec3d& originKm,
                                  const math::Vec3d& dirNormalized,
                                  double maxRangeKm,
                                  const SphereTarget* targets,
                                  std::size_t targetCount) {
  RaycastHit best{};
  best.hit = false;
  best.tKm = std::max(0.0, maxRangeKm);
  best.pointKm = originKm + dirNormalized * best.tKm;

  if (!targets || targetCount == 0) return best;

  const double maxR = std::max(0.0, maxRangeKm);
  double bestT = maxR;

  for (std::size_t i = 0; i < targetCount; ++i) {
    const SphereTarget& t = targets[i];
    const math::Vec3d toC = t.centerKm - originKm;
    const double dist2 = toC.lengthSq();
    if (dist2 < 1e-12) continue;

    // Aim cone filter (optional).
    if (t.minAimCos > -0.5) {
      const double dist = std::sqrt(dist2);
      const double aim = math::dot(dirNormalized, toC / dist);
      if (aim < t.minAimCos) continue;
    }

    // Cheap early reject: center projection outside [0, maxRange + radius].
    const double tProj = math::dot(toC, dirNormalized);
    if (tProj < -t.radiusKm) continue;
    if (tProj > maxR + t.radiusKm) continue;

    double tEnter = 0.0;
    if (!raySphereIntersectKm(originKm, dirNormalized, t.centerKm, t.radiusKm, tEnter)) continue;
    if (tEnter > maxR) continue;

    if (tEnter < bestT) {
      bestT = tEnter;
      best.hit = true;
      best.tKm = tEnter;
      best.pointKm = originKm + dirNormalized * tEnter;
      best.kind = t.kind;
      best.index = t.index;
      best.id = t.id;
    }
  }

  if (!best.hit) {
    best.tKm = maxR;
    best.pointKm = originKm + dirNormalized * maxR;
  }
  return best;
}

bool segmentHitsSphereKm(const math::Vec3d& aKm,
                         const math::Vec3d& bKm,
                         const math::Vec3d& centerKm,
                         double radiusKm) {
  const math::Vec3d ab = bKm - aKm;
  const double abLenSq = ab.lengthSq();
  const double r2 = radiusKm * radiusKm;
  if (abLenSq < 1e-12) {
    return (aKm - centerKm).lengthSq() <= r2;
  }

  const double t = std::clamp(math::dot(centerKm - aKm, ab) / abLenSq, 0.0, 1.0);
  const math::Vec3d closest = aKm + ab * t;
  return (closest - centerKm).lengthSq() <= r2;
}

static double segmentClosestT(const math::Vec3d& aKm,
                              const math::Vec3d& bKm,
                              const math::Vec3d& pKm) {
  const math::Vec3d ab = bKm - aKm;
  const double abLenSq = ab.lengthSq();
  if (abLenSq < 1e-12) return 0.0;
  return std::clamp(math::dot(pKm - aKm, ab) / abLenSq, 0.0, 1.0);
}

void stepProjectiles(std::vector<Projectile>& projectiles,
                     double dtSim,
                     const SphereTarget* targets,
                     std::size_t targetCount,
                     std::vector<ProjectileHit>& outHits) {
  if (dtSim <= 0.0) return;
  if (projectiles.empty()) return;

  outHits.clear();

  for (auto& p : projectiles) {
    if (p.ttlSimSec <= 0.0) continue;

    const math::Vec3d a = p.posKm;
    const math::Vec3d b = p.posKm + p.velKmS * dtSim;

    p.prevKm = a;
    p.posKm = b;
    p.ttlSimSec -= dtSim;

    if (!targets || targetCount == 0) continue;

    // Decide which target kinds this projectile can hit.
    const bool allowHitPlayer = !p.fromPlayer;
    const bool allowHitShips = true;
    const bool allowHitAsteroids = p.fromPlayer; // fun: player slugs can smack rocks

    for (std::size_t ti = 0; ti < targetCount; ++ti) {
      const SphereTarget& t = targets[ti];

      if (!allowHitPlayer && t.kind == CombatTargetKind::Player) continue;
      if (!allowHitShips && t.kind == CombatTargetKind::Ship) continue;
      if (!allowHitAsteroids && t.kind == CombatTargetKind::Asteroid) continue;

      // Avoid immediate self-hits when the shooter has a real id.
      if (p.shooterId != 0 && t.id != 0 && t.id == p.shooterId) continue;

      const double hitRadiusKm = std::max(0.0, t.radiusKm) + std::max(0.0, p.radiusKm);
      if (!segmentHitsSphereKm(a, b, t.centerKm, hitRadiusKm)) continue;

      const double tt = segmentClosestT(a, b, t.centerKm);
      const math::Vec3d hitPoint = a + (b - a) * tt;

      ProjectileHit h{};
      h.kind = t.kind;
      h.targetIndex = t.index;
      h.targetId = t.id;
      h.pointKm = hitPoint;
      h.dmg = p.dmg;
      h.fromPlayer = p.fromPlayer;
      h.shooterId = p.shooterId;
      outHits.push_back(h);

      p.ttlSimSec = 0.0;
      break;
    }
  }
}

double weaponHeatDelta(const WeaponDef& w, int distributorMk) {
  // Better distributors run cooler.
  const int dMk = std::clamp(distributorMk, 1, 3);
  const double heatFactor = std::clamp(1.08 - 0.08 * (double)dMk, 0.80, 1.10);
  return std::max(0.0, w.heatPerShot * heatFactor);
}

FireResult tryFireWeapon(Ship& shooter,
                         WeaponType weapon,
                         double cooldownSimSec,
                         int distributorMk,
                         core::u64 shooterId,
                         bool fromPlayer,
                         const SphereTarget* beamTargets,
                         std::size_t beamTargetCount) {
  FireResult out{};
  out.fired = false;
  out.newCooldownSimSec = cooldownSimSec;

  if (cooldownSimSec > 0.0) return out;

  const WeaponDef& w = weaponDef(weapon);
  out.fired = true;
  out.newCooldownSimSec = w.cooldownSimSec;
  out.heatDelta = weaponHeatDelta(w, distributorMk);

  if (w.beam) {
    const math::Vec3d origin = shooter.positionKm();
    math::Vec3d dir = shooter.forward();
    if (dir.lengthSq() < 1e-12) dir = {0, 0, 1};
    dir = dir.normalized();

    const double rangeKm = std::max(0.0, w.rangeKm);
    const RaycastHit hit = raycastNearestSphereKm(origin, dir, rangeKm, beamTargets, beamTargetCount);

    out.hasBeam = true;
    out.beam.aKm = origin;
    out.beam.bKm = hit.pointKm;
    out.beam.r = w.r;
    out.beam.g = w.g;
    out.beam.b = w.b;

    out.hit = hit.hit;
    if (hit.hit) {
      out.hitKind = hit.kind;
      out.hitIndex = hit.index;
      out.hitId = hit.id;
      out.hitPointKm = hit.pointKm;
    }
    return out;
  }

  // Projectile weapon.
  const double muzzleSpeedKmS = std::max(1e-6, w.projSpeedKmS);
  const double rangeKm = std::max(0.0, w.rangeKm);
  const double ttlSim = rangeKm / muzzleSpeedKmS;

  math::Vec3d fwd = shooter.forward();
  if (fwd.lengthSq() < 1e-12) fwd = {0, 0, 1};
  fwd = fwd.normalized();

  const math::Vec3d spawnKm = shooter.positionKm() + fwd * 400.0;

  Projectile p{};
  p.prevKm = spawnKm;
  p.posKm = spawnKm;
  p.velKmS = shooter.velocityKmS() + fwd * muzzleSpeedKmS;
  p.r = w.r; p.g = w.g; p.b = w.b;
  p.ttlSimSec = ttlSim;
  p.radiusKm = (weapon == WeaponType::Railgun) ? 520.0 : 700.0;
  p.dmg = w.dmg;
  p.fromPlayer = fromPlayer;
  p.shooterId = shooterId;

  out.hasProjectile = true;
  out.projectile = p;

  // Small recoil impulse (matches historical tuning in stellar_game).
  const double recoil = (weapon == WeaponType::Railgun) ? 0.003 : 0.002;
  shooter.setVelocityKmS(shooter.velocityKmS() - fwd * recoil);

  return out;
}

} // namespace stellar::sim
