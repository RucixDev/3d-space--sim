#include "stellar/sim/Combat.h"

#include "stellar/sim/Ballistics.h"

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

static math::Vec3d rotateTowards(const math::Vec3d& fromDir,
                                 const math::Vec3d& toDir,
                                 double maxAngleRad) {
  const math::Vec3d f = fromDir.normalized();
  const math::Vec3d t = toDir.normalized();
  if (f.lengthSq() < 1e-12) return t;
  if (t.lengthSq() < 1e-12) return f;

  const double c = std::clamp(math::dot(f, t), -1.0, 1.0);
  const double ang = std::acos(c);
  if (ang <= 1e-9) return t;
  if (maxAngleRad <= 0.0) return f;
  if (ang <= maxAngleRad) return t;

  math::Vec3d axis = math::cross(f, t);
  const double al2 = axis.lengthSq();
  if (al2 < 1e-18) {
    // Nearly opposite; just keep current direction.
    return f;
  }
  axis = axis / std::sqrt(al2);

  const double ca = std::cos(maxAngleRad);
  const double sa = std::sin(maxAngleRad);

  // Rodrigues rotation formula.
  const math::Vec3d v = f;
  const math::Vec3d v2 = math::cross(axis, v);
  return (v * ca + v2 * sa).normalized();
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

void stepMissiles(std::vector<Missile>& missiles,
                  double dtSim,
                  const SphereTarget* targets,
                  std::size_t targetCount,
                  std::vector<MissileDetonation>& outDetonations,
                  std::vector<MissileHit>& outHits) {
  if (dtSim <= 0.0) return;
  if (missiles.empty()) return;

  outDetonations.clear();
  outHits.clear();

  for (auto& m : missiles) {
    if (m.ttlSimSec <= 0.0) continue;

    const math::Vec3d a = m.posKm;

    // --- Guidance ---
    math::Vec3d v = m.velKmS;
    double speed = v.length();
    if (speed < 1e-6) {
      speed = 1.0;
      v = {0, 0, 1};
    }

    math::Vec3d desiredDir = v.normalized();
    const SphereTarget* tgt = nullptr;

    // "Lock" score used to compare against decoys (inverse-square falloff).
    double targetScore = 0.0;

    if (m.hasTarget && targets && targetCount > 0) {
      for (std::size_t i = 0; i < targetCount; ++i) {
        const SphereTarget& t = targets[i];
        if (t.kind != m.targetKind) continue;
        if (t.id != m.targetId) continue;
        tgt = &t;
        break;
      }

      if (tgt) {
        // Lead solve (missile treated as constant-speed projectile).
        const auto lead = solveProjectileLead(m.posKm,
                                              /*shooterVelKmS=*/{0, 0, 0},
                                              tgt->centerKm,
                                              tgt->velKmS,
                                              speed,
                                              /*maxTimeSec=*/1.0e6,
                                              /*minTimeSec=*/1.0e-3);
        if (lead && lead->aimDirWorld.lengthSq() > 1e-12) {
          desiredDir = lead->aimDirWorld;
        } else {
          const math::Vec3d to = tgt->centerKm - m.posKm;
          if (to.lengthSq() > 1e-12) desiredDir = to.normalized();
        }

        const math::Vec3d toTgt = tgt->centerKm - m.posKm;
        targetScore = 1.0 / (toTgt.lengthSq() + 1.0e-9);
      }

      // --- Countermeasure / decoy attraction ---
      const SphereTarget* bestDecoy = nullptr;
      double bestDecoyScore = 0.0;

      if (targetScore > 0.0) {
        const math::Vec3d fwd = v.normalized();
        const double fovCos = std::clamp(m.seekerFovCos, -1.0, 1.0);

        for (std::size_t i = 0; i < targetCount; ++i) {
          const SphereTarget& t = targets[i];
          if (t.kind != CombatTargetKind::Decoy) continue;

          const double strength = (m.seeker == MissileSeekerType::Radar) ? t.decoyRadar : t.decoyHeat;
          if (strength <= 0.0) continue;

          const math::Vec3d to = t.centerKm - m.posKm;
          const double distSq = to.lengthSq();
          if (distSq < 1e-12) continue;

          const math::Vec3d toDir = to.normalized();
          const double cosAng = math::dot(fwd, toDir);
          if (cosAng < fovCos) continue;

          const double score = (strength * std::max(0.0, cosAng)) / (distSq + 1.0e-9);
          if (score > bestDecoyScore) {
            bestDecoyScore = score;
            bestDecoy = &t;
          }
        }
      }

      const double resist = std::max(0.0, m.decoyResistance);
      if (bestDecoy && bestDecoyScore > targetScore * resist) {
        const auto lead = solveProjectileLead(m.posKm,
                                              /*shooterVelKmS=*/{0, 0, 0},
                                              bestDecoy->centerKm,
                                              bestDecoy->velKmS,
                                              speed,
                                              /*maxTimeSec=*/1.0e6,
                                              /*minTimeSec=*/1.0e-3);
        if (lead && lead->aimDirWorld.lengthSq() > 1e-12) {
          desiredDir = lead->aimDirWorld;
        } else {
          const math::Vec3d to = bestDecoy->centerKm - m.posKm;
          if (to.lengthSq() > 1e-12) desiredDir = to.normalized();
        }
      }
    }


    const double maxTurn = std::max(0.0, m.turnRateRadS) * dtSim;
    const math::Vec3d newDir = rotateTowards(v.normalized(), desiredDir, maxTurn);
    m.velKmS = newDir * speed;

    const math::Vec3d b = m.posKm + m.velKmS * dtSim;

    m.prevKm = a;
    m.posKm = b;
    m.ttlSimSec -= dtSim;

    // --- Collision & detonation ---
    bool detonated = false;
    math::Vec3d detPoint = b;

    if (targets && targetCount > 0) {
      const bool allowHitPlayer = !m.fromPlayer;
      const bool allowHitShips = true;
      const bool allowHitAsteroids = true;

      double bestT = 2.0; // [0,1]

      for (std::size_t ti = 0; ti < targetCount; ++ti) {
        const SphereTarget& t = targets[ti];

        if (!allowHitPlayer && t.kind == CombatTargetKind::Player) continue;
        if (!allowHitShips && t.kind == CombatTargetKind::Ship) continue;
        if (!allowHitAsteroids && t.kind == CombatTargetKind::Asteroid) continue;

        // Avoid immediate self-hits when the shooter has a real id.
        if (m.shooterId != 0 && t.id != 0 && t.id == m.shooterId) continue;

        const double hitRadiusKm = std::max(0.0, t.radiusKm) + std::max(0.0, m.radiusKm);
        if (!segmentHitsSphereKm(a, b, t.centerKm, hitRadiusKm)) continue;

        const double tt = segmentClosestT(a, b, t.centerKm);
        if (tt < bestT) {
          bestT = tt;
          detonated = true;
          detPoint = a + (b - a) * tt;
        }
      }
    }

    // Expired missiles disappear without a detonation (keeps background clutter down).
    if (!detonated) {
      if (m.ttlSimSec <= 0.0) {
        m.ttlSimSec = 0.0;
      }
      continue;
    }

    // Detonate.
    MissileDetonation det{};
    det.pointKm = detPoint;
    det.blastRadiusKm = std::max(0.0, m.blastRadiusKm);
    det.baseDmg = std::max(0.0, m.dmg);
    det.fromPlayer = m.fromPlayer;
    det.shooterId = m.shooterId;
    outDetonations.push_back(det);

    // Splash hits.
    if (targets && targetCount > 0 && det.blastRadiusKm > 1e-6 && det.baseDmg > 1e-9) {
      const bool allowHitPlayer = !m.fromPlayer;
      const bool allowHitShips = true;
      const bool allowHitAsteroids = true;

      for (std::size_t ti = 0; ti < targetCount; ++ti) {
        const SphereTarget& t = targets[ti];

        if (!allowHitPlayer && t.kind == CombatTargetKind::Player) continue;
        if (!allowHitShips && t.kind == CombatTargetKind::Ship) continue;
        if (!allowHitAsteroids && t.kind == CombatTargetKind::Asteroid) continue;

        // Avoid self splash when shooter has a real id.
        if (m.shooterId != 0 && t.id != 0 && t.id == m.shooterId) continue;

        const double dist = (t.centerKm - detPoint).length();
        const double surface = std::max(0.0, dist - std::max(0.0, t.radiusKm));
        const double k = std::clamp(1.0 - surface / det.blastRadiusKm, 0.0, 1.0);
        if (k <= 0.0) continue;

        MissileHit h{};
        h.kind = t.kind;
        h.targetIndex = t.index;
        h.targetId = t.id;
        h.pointKm = detPoint;
        h.dmg = det.baseDmg * k;
        h.fromPlayer = m.fromPlayer;
        h.shooterId = m.shooterId;
        outHits.push_back(h);
      }
    }

    m.ttlSimSec = 0.0;
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

  // Guided weapon (missile).
  if (w.guided) {
    const double muzzleSpeedKmS = std::max(1e-6, w.projSpeedKmS);
    const double rangeKm = std::max(0.0, w.rangeKm);
    const double ttlSim = rangeKm / muzzleSpeedKmS;

    math::Vec3d fwd = shooter.forward();
    if (fwd.lengthSq() < 1e-12) fwd = {0, 0, 1};
    fwd = fwd.normalized();

    // Lock onto a target in front of the shooter.
    bool hasLock = false;
    CombatTargetKind lockKind = CombatTargetKind::Ship;
    core::u64 lockId = 0;

    if (beamTargets && beamTargetCount > 0) {
      double bestAim = -1.0;
      double bestDist = 1.0e30;
      const math::Vec3d origin = shooter.positionKm();
      const math::Vec3d dir = fwd;

      for (std::size_t i = 0; i < beamTargetCount; ++i) {
        const SphereTarget& t = beamTargets[i];
        if (t.kind != CombatTargetKind::Ship && t.kind != CombatTargetKind::Player) continue;
        if (fromPlayer && t.kind == CombatTargetKind::Player) continue;
        if (shooterId != 0 && t.id != 0 && t.id == shooterId) continue;

        const math::Vec3d to = t.centerKm - origin;
        const double dist2 = to.lengthSq();
        if (dist2 < 1e-12) continue;
        const double dist = std::sqrt(dist2);
        if (dist > rangeKm) continue;

        const double aim = math::dot(dir, to / dist);
        // Don't lock behind.
        if (aim <= 0.0) continue;

        const double req = (t.minAimCos > -0.5) ? t.minAimCos : -1.0;
        if (aim < req) continue;

        // Prefer highest aim, then nearest.
        if (aim > bestAim + 1e-6 || (std::abs(aim - bestAim) <= 1e-6 && dist < bestDist)) {
          bestAim = aim;
          bestDist = dist;
          hasLock = true;
          lockKind = t.kind;
          lockId = t.id;
        }
      }
    }

    const math::Vec3d spawnKm = shooter.positionKm() + fwd * 520.0;

    Missile m{};
    m.prevKm = spawnKm;
    m.posKm = spawnKm;
    m.velKmS = shooter.velocityKmS() + fwd * muzzleSpeedKmS;
    m.r = w.r;
    m.g = w.g;
    m.b = w.b;
    m.ttlSimSec = ttlSim;
    m.radiusKm = 760.0;
    m.dmg = w.dmg;
    m.blastRadiusKm = std::max(0.0, w.blastRadiusKm);
    m.turnRateRadS = std::max(0.0, w.turnRateRadS);

    // Seeker tuning (decoys + field-of-view).
    // Heat seekers can be lured by flares; radar seekers by chaff.
    // These parameters intentionally remain simple and deterministic.
    if (weapon == WeaponType::RadarMissile) {
      m.seeker = MissileSeekerType::Radar;
      // Narrower seeker cone (radians -> cosine).
      m.seekerFovCos = std::cos(0.80);
      // Slightly resistant to chaff, but can still be fooled.
      m.decoyResistance = 0.92;
    } else {
      // Default guided weapon: heat seeker.
      m.seeker = MissileSeekerType::Heat;
      // Wider cone than radar.
      m.seekerFovCos = std::cos(0.95);
      // Easier to decoy with flares.
      m.decoyResistance = 0.80;
    }
    m.fromPlayer = fromPlayer;
    m.shooterId = shooterId;

    m.hasTarget = hasLock;
    m.targetKind = lockKind;
    m.targetId = lockId;

    out.hasMissile = true;
    out.missile = m;

    // Minimal recoil.
    shooter.setVelocityKmS(shooter.velocityKmS() - fwd * 0.0015);

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
