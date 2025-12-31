#include "stellar/sim/FlightController.h"

#include "stellar/sim/Ballistics.h"

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

// Build a quaternion from an orthonormal basis.
// The basis maps ship-local axes to world axes:
//   local +X -> rightW, local +Y -> upW, local +Z -> fwdW.
static math::Quatd quatFromBasis(const math::Vec3d& rightW,
                                 const math::Vec3d& upW,
                                 const math::Vec3d& fwdW) {
  // Column-major rotation matrix.
  const double m00 = rightW.x;
  const double m01 = upW.x;
  const double m02 = fwdW.x;
  const double m10 = rightW.y;
  const double m11 = upW.y;
  const double m12 = fwdW.y;
  const double m20 = rightW.z;
  const double m21 = upW.z;
  const double m22 = fwdW.z;

  const double trace = m00 + m11 + m22;
  math::Quatd q{1, 0, 0, 0};

  if (trace > 0.0) {
    const double s = std::sqrt(trace + 1.0) * 2.0;
    q.w = 0.25 * s;
    q.x = (m21 - m12) / s;
    q.y = (m02 - m20) / s;
    q.z = (m10 - m01) / s;
  } else if ((m00 > m11) && (m00 > m22)) {
    const double s = std::sqrt(1.0 + m00 - m11 - m22) * 2.0;
    q.w = (m21 - m12) / s;
    q.x = 0.25 * s;
    q.y = (m01 + m10) / s;
    q.z = (m02 + m20) / s;
  } else if (m11 > m22) {
    const double s = std::sqrt(1.0 + m11 - m00 - m22) * 2.0;
    q.w = (m02 - m20) / s;
    q.x = (m01 + m10) / s;
    q.y = 0.25 * s;
    q.z = (m12 + m21) / s;
  } else {
    const double s = std::sqrt(1.0 + m22 - m00 - m11) * 2.0;
    q.w = (m10 - m01) / s;
    q.x = (m02 + m20) / s;
    q.y = (m12 + m21) / s;
    q.z = 0.25 * s;
  }

  return q.normalized();
}

static math::Vec3d attitudeTorqueLocal(const Ship& ship,
                                             const math::Vec3d& desiredForwardWorld,
                                             const AttitudeControlParams& attitude,
                                             const math::Vec3d* desiredUpWorld) {
  math::Vec3d fwdW = desiredForwardWorld;
  if (fwdW.lengthSq() <= 1e-12) {
    fwdW = ship.forward();
  }
  fwdW = fwdW.normalized();

  // When we care about roll/up alignment, drive a single 3D orientation error
  // instead of decomposing yaw/pitch/roll. This avoids 180° flip ambiguity.
  if (attitude.alignUp && desiredUpWorld) {
    math::Vec3d upRefW = *desiredUpWorld;
    if (upRefW.lengthSq() <= 1e-12) {
      upRefW = {0, 1, 0};
    }
    upRefW = upRefW.normalized();

    math::Vec3d rightW = math::cross(upRefW, fwdW);
    if (rightW.lengthSq() <= 1e-12) {
      // Desired up is parallel to forward; pick a fallback reference up.
      upRefW = (std::abs(fwdW.y) < 0.9) ? math::Vec3d{0, 1, 0} : math::Vec3d{1, 0, 0};
      rightW = math::cross(upRefW, fwdW);
    }
    rightW = rightW.normalized();
    const math::Vec3d upW = math::cross(fwdW, rightW).normalized();

    const math::Quatd qDes = quatFromBasis(rightW, upW, fwdW);

    // Body-frame orientation error.
    math::Quatd qErr = ship.orientation().conjugate() * qDes;
    if (qErr.w < 0.0) {
      qErr = math::Quatd{-qErr.w, -qErr.x, -qErr.y, -qErr.z};
    }

    // For small angles, 2*imag(qErr) ~= axis*angle. Works well as a control error.
    const math::Vec3d err = math::Vec3d{qErr.x, qErr.y, qErr.z} * 2.0;
    return {
      std::clamp(err.x * attitude.faceGain, -1.0, 1.0),
      std::clamp(err.y * attitude.faceGain, -1.0, 1.0),
      std::clamp(err.z * attitude.rollGain, -1.0, 1.0),
    };
  }

  // Default: align forward only (roll is free).
  const math::Vec3d desiredFwdLocal = ship.orientation().conjugate().rotate(fwdW);
  const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
  const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);

  return {
    std::clamp(pitchErr * attitude.faceGain, -1.0, 1.0),
    std::clamp(yawErr * attitude.faceGain, -1.0, 1.0),
    0.0,
  };
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
  // Use a quaternion-based controller when aligning both forward + up (avoids
  // 180° roll ambiguity during docking, and converges more reliably than the
  // yaw/pitch + roll decomposition).
  out.input.torqueLocal = attitudeTorqueLocal(ship, desiredForwardWorld, attitude, desiredUpWorld);

  return out;
}

FlightControlOutput approachTargetIntercept(const Ship& ship,
                                           const math::Vec3d& targetPosKm,
                                           const math::Vec3d& targetVelKmS,
                                           const FlightControlParams& params,
                                           const AttitudeControlParams& attitude,
                                           const math::Vec3d& desiredForwardWorld,
                                           const InterceptCourseParams& intercept,
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

  // Intercept-course lead direction.
  //
  // Use a constant-speed intercept solve (same math as projectile lead) to find an
  // approach direction that "cuts off" lateral motion. This is only engaged when:
  //   - enabled
  //   - we are closing (desiredSpeed > 0)
  //   - desiredSpeed is large enough to avoid noisy solves near the target
  math::Vec3d approachDir = toN;
  if (intercept.enabled && desiredSpeed > 0.0 && desiredSpeed >= intercept.minSpeedKmS) {
    const double solveSpeed = std::max(1e-6,
                                      intercept.useMaxSpeedForSolve ? params.maxSpeedKmS : desiredSpeed);

    // Solve assuming the pursuer can travel at `solveSpeed` in world space. This is an
    // approximation (ships are acceleration-limited), but works well for guidance.
    const auto lead = solveProjectileLead(ship.positionKm(),
                                          /*shooterVelKmS=*/{0,0,0},
                                          targetPosKm,
                                          targetVelKmS,
                                          solveSpeed,
                                          /*maxTimeSec=*/std::max(0.0, intercept.maxLeadTimeSec),
                                          /*minTimeSec=*/1.0e-3);
    if (lead) {
      approachDir = lead->aimDirWorld;
      if (approachDir.lengthSq() < 1e-12) approachDir = toN;
    }
  }

  out.desiredVelKmS = targetVelKmS + approachDir * desiredSpeed;
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
  out.input.torqueLocal = attitudeTorqueLocal(ship, desiredForwardWorld, attitude, desiredUpWorld);

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

FlightControlOutput chaseTargetIntercept(const Ship& ship,
                                         const math::Vec3d& targetPosKm,
                                         const math::Vec3d& targetVelKmS,
                                         const FlightControlParams& params,
                                         const AttitudeControlParams& attitude,
                                         const InterceptCourseParams& intercept) {
  math::Vec3d to = (targetPosKm - ship.positionKm());
  math::Vec3d toN = to.normalized();
  if (toN.lengthSq() <= 1e-12) {
    toN = ship.forward().normalized();
  }
  return approachTargetIntercept(ship, targetPosKm, targetVelKmS, params, attitude, toN, intercept, nullptr);
}

} // namespace stellar::sim
