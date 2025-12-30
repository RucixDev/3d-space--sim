#include "stellar/sim/FlightController.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static math::Vec3d clampMagnitude(const math::Vec3d& v, double maxLen) {
  const double len = v.length();
  if (len <= maxLen || len <= 1e-12) return v;
  return v * (maxLen / len);
}

static math::Vec3d clampComponents(const math::Vec3d& v, double lo, double hi) {
  return {
    std::clamp(v.x, lo, hi),
    std::clamp(v.y, lo, hi),
    std::clamp(v.z, lo, hi),
  };
}

static double safeSqrt(double x) {
  return (x <= 0.0) ? 0.0 : std::sqrt(x);
}

FlightControlOutput approachTarget(const Ship& ship,
                                  const math::Vec3d& targetPosKm,
                                  const math::Vec3d& targetVelKmS,
                                  const FlightControlParams& params,
                                  const AttitudeControlParams& attitude,
                                  const math::Vec3d& desiredForwardWorld,
                                  const math::Vec3d* desiredUpWorld) {
  FlightControlOutput out{};

  out.input.dampers = params.dampers;

  // --- Translation ---
  const math::Vec3d to = targetPosKm - ship.positionKm();
  const double dist = to.length();
  out.distKm = dist;

  math::Vec3d toN = (dist > 1e-9) ? (to / dist) : math::Vec3d{0, 0, 1};

  const double desiredDist = std::max(0.0, params.desiredDistKm);
  const double rem = std::max(0.0, dist - desiredDist);

  const double baseCap = std::max(1e-12, ship.maxLinearAccelKmS2() * params.accelScale);
  const double boostCap = std::max(baseCap, ship.maxLinearAccelBoostKmS2() * params.accelScale);

  // Speed profile: proportional approach with a conservative stopping-distance clamp.
  double desiredSpeed = std::min(params.maxSpeedKmS, params.speedGain * rem);
  const double stopSpeed = safeSqrt(2.0 * baseCap * rem);
  desiredSpeed = std::min(desiredSpeed, stopSpeed);

  // If we are deep inside the desired distance, gently back away.
  if (desiredDist > 1e-6 && dist < desiredDist * params.backoffFrac) {
    desiredSpeed = -std::abs(params.backoffSpeedKmS);
  }

  out.desiredVelKmS = targetVelKmS + toN * desiredSpeed;
  const math::Vec3d dv = out.desiredVelKmS - ship.velocityKmS();

  math::Vec3d aWanted = dv * params.velGain;

  double capUsed = baseCap;
  bool useBoost = false;
  if (params.allowBoost) {
    const double want = aWanted.length();
    if (want > baseCap * 0.98 && boostCap > baseCap * 1.01) {
      capUsed = boostCap;
      useBoost = true;
    }
  }

  out.usedBoost = useBoost;
  out.input.boost = useBoost;
  out.input.brake = false;

  const math::Vec3d aCmd = clampMagnitude(aWanted, capUsed);
  const math::Vec3d thrustWorld = (capUsed > 1e-12) ? (aCmd / capUsed) : math::Vec3d{0, 0, 0};
  const math::Vec3d thrustLocal = ship.orientation().conjugate().rotate(thrustWorld);
  out.input.thrustLocal = clampComponents(thrustLocal, -1.0, 1.0);

  // --- Attitude ---
  math::Vec3d fwdW = desiredForwardWorld;
  if (fwdW.lengthSq() <= 1e-12) {
    fwdW = ship.forward();
  }
  fwdW = fwdW.normalized();

  const math::Vec3d desiredFwdLocal = ship.orientation().conjugate().rotate(fwdW);
  const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
  const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);

  out.input.torqueLocal.x = std::clamp(pitchErr * attitude.faceGain, -1.0, 1.0);
  out.input.torqueLocal.y = std::clamp(yawErr * attitude.faceGain, -1.0, 1.0);
  out.input.torqueLocal.z = 0.0;

  if (attitude.alignUp && desiredUpWorld) {
    math::Vec3d upW = *desiredUpWorld;
    if (upW.lengthSq() > 1e-12) upW = upW.normalized();

    const math::Vec3d desiredUpLocal = ship.orientation().conjugate().rotate(upW);
    const double rollErr = std::atan2(desiredUpLocal.x, desiredUpLocal.y);
    out.input.torqueLocal.z = std::clamp(rollErr * attitude.rollGain, -1.0, 1.0);
  }

  return out;
}

FlightControlOutput chaseTarget(const Ship& ship,
                                const math::Vec3d& targetPosKm,
                                const math::Vec3d& targetVelKmS,
                                const FlightControlParams& params,
                                const AttitudeControlParams& attitude) {
  math::Vec3d to = (targetPosKm - ship.positionKm());
  math::Vec3d toN = to.normalized();
  if (toN.lengthSq() <= 1e-12) {
    toN = ship.forward().normalized();
  }
  return approachTarget(ship, targetPosKm, targetVelKmS, params, attitude, toN, nullptr);
}

} // namespace stellar::sim
