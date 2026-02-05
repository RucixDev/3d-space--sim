#include "stellar/sim/ThreatAvoidance.h"

#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace stellar::sim {

static double clamp01(double x) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  return x;
}

static math::Vec3d unitPerpFromSeed(const math::Vec3d& axisIn, core::u64 seed) {
  math::Vec3d axis = axisIn;
  if (axis.lengthSq() < 1e-12) axis = {0, 0, 1};
  axis = axis.normalized();

  // Build a deterministic orthonormal basis around `axis`.
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

static math::Vec3d collisionJinkDirWorld(const Ship& ship,
                                        const ProximityFieldKm& obstacles,
                                        const CollisionWarningResult& c,
                                        core::u64 seed) {
  // Base axis is current velocity direction (or forward if nearly stationary).
  math::Vec3d vDir = ship.velocityKmS();
  if (vDir.lengthSq() < 1e-18) vDir = ship.forward();
  vDir = math::safeNormalized(vDir, math::Vec3d{0, 0, 1}, 1e-12);

  const SphereObstacleKm* obs = obstacles.obstacle(c.obstacleId);
  if (!obs) {
    return unitPerpFromSeed(vDir, seed ^ 0x8b1d4e7a3c2f91a5ull);
  }

  // Take the obstacle-center offset, then project into the plane perpendicular
  // to the current velocity direction. Accelerating along this direction tends
  // to increase miss distance without increasing closure along the current path.
  const math::Vec3d rel = ship.positionKm() - obs->centerKm;
  math::Vec3d lat = rel - vDir * math::dot(rel, vDir);
  if (lat.lengthSq() < 1e-12) {
    lat = unitPerpFromSeed(vDir, seed ^ 0x1d6d3d92d61f4b5bull);
  }
  lat = math::safeNormalized(lat, unitPerpFromSeed(vDir, seed ^ 0x9e3779b97f4a7c15ull), 1e-12);

  // Ensure we push away from the obstacle (should already be true for the above
  // projection, but be defensive).
  if (math::dot(lat, rel) < 0.0) lat = lat * -1.0;
  return lat;
}

ThreatAvoidanceResult computeThreatAvoidance(const Ship& ship,
                                            const ProximityFieldKm* obstacles,
                                            double maxDecelKmS2,
                                            const Missile* missiles,
                                            std::size_t missileCount,
                                            CombatTargetKind targetKind,
                                            core::u64 targetId,
                                            core::u64 seed,
                                            const ThreatAvoidanceParams& params,
                                            int ignoreObstacleId) {
  ThreatAvoidanceResult out{};

  const math::Vec3d pos = ship.positionKm();
  const math::Vec3d vel = ship.velocityKmS();
  const double speed = vel.length();

  // ------- Collision prediction -------
  math::Vec3d collisionDir{0, 0, 0};
  double wCollision = 0.0;

  if (params.collisionEnable && obstacles && !obstacles->empty() &&
      speed >= std::max(0.0, params.collisionMinSpeedForJinkKmS)) {
    out.collision = computeCollisionWarning(*obstacles, pos, vel, maxDecelKmS2, params.collision,
                                            ignoreObstacleId);

    if (out.collision.hasImpact && out.collision.hazard01 >= std::max(0.0, params.collisionEngageHazard01)) {
      out.collisionActive = true;
      collisionDir = collisionJinkDirWorld(ship, *obstacles, out.collision, seed);
      wCollision = std::max(0.0, params.collisionWeight) * clamp01(out.collision.hazard01);
    }
  }

  // ------- Missile threat -------
  math::Vec3d missileDir{0, 0, 0};
  double wMissile = 0.0;

  if (params.missileEnable && missiles && missileCount > 0) {
    out.missileThreat = nearestInboundMissile(missiles, missileCount, targetKind, targetId, pos, vel,
                                              params.missileThreat);

    const double engageTti = std::max(0.0, params.missileEngageTtiSec);
    if (out.missileThreat.inbound && engageTti > 1e-9 && out.missileThreat.ttiSec <= engageTti) {
      const std::size_t idx = out.missileThreat.missileIndex;
      if (idx < missileCount) {
        const Missile& m = missiles[idx];
        const core::u64 planSeed = seed ^ (static_cast<core::u64>(idx) * 0x9e3779b97f4a7c15ull);
        out.missilePlan = planMissileEvasion(m, pos, vel, planSeed, params.missileEvasion);
        if (out.missilePlan.valid) {
          out.missileActive = true;
          missileDir = out.missilePlan.dirWorld;

          const double ramp = clamp01((engageTti - std::max(0.0, out.missileThreat.ttiSec)) / engageTti);
          wMissile = std::max(0.0, params.missileWeight) * ramp;
        }
      }
    }
  }

  // ------- Blend -------
  math::Vec3d sum{0, 0, 0};
  double wMax = 0.0;

  if (out.collisionActive) {
    sum += collisionDir * wCollision;
    wMax = std::max(wMax, wCollision);
  }
  if (out.missileActive) {
    sum += missileDir * wMissile;
    wMax = std::max(wMax, wMissile);
  }

  if (sum.lengthSq() <= 1e-12) {
    return out;
  }

  math::Vec3d dir = sum.normalized();

  // Prefer classic lateral jinks when we have meaningful velocity.
  if (params.preferLateralJink && speed > 1e-6) {
    const math::Vec3d vDir = math::safeNormalized(vel, ship.forward(), 1e-12);
    dir = dir - vDir * math::dot(dir, vDir);
    dir = math::safeNormalized(dir, unitPerpFromSeed(vDir, seed ^ 0x7f4a7c159e3779b9ull), 1e-12);
  }

  out.dirWorld = dir;
  out.thrust01 = std::clamp(wMax, 0.0, std::max(0.0, params.maxThrust01));

  if (out.thrust01 > 1e-12) {
    out.active = true;
    out.input.dampers = true;
    out.input.thrustLocal = ship.orientation().conjugate().rotate(dir) * out.thrust01;

    if (params.collisionEnable && out.collisionActive) {
      const double brakeTh = std::max(0.0, params.collisionBrakeEngageHazard01);
      out.input.brake = (brakeTh > 0.0) && (out.collision.hazard01 >= brakeTh);
    }

    if (params.allowBoost) {
      const bool boostForCollision = out.collisionActive && (out.collision.hazard01 >= std::max(0.0, params.boostEngageHazard01));
      const bool boostForMissile = out.missileActive && (out.missileThreat.ttiSec <= std::max(0.0, params.boostEngageMissileTtiSec));
      out.input.boost = boostForCollision || boostForMissile;
    }
  }

  return out;
}

} // namespace stellar::sim
