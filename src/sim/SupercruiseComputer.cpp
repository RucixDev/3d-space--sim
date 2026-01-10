#include "stellar/sim/SupercruiseComputer.h"

#include "stellar/math/Math.h"
#include "stellar/sim/Ballistics.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static math::Vec3d clampMagnitude(const math::Vec3d& v, double maxLen) {
  const double len = v.length();
  if (len <= maxLen || len <= 1e-12) return v;
  return v * (maxLen / len);
}

static double safeAcos(double x) {
  return std::acos(std::clamp(x, -1.0, 1.0));
}

static math::Vec3d nlerpDir(const math::Vec3d& a, const math::Vec3d& b, double t) {
  const math::Vec3d v = a * (1.0 - t) + b * t;
  const double len = v.length();
  if (len <= 1e-12) return a;
  return v * (1.0 / len);
}

SupercruiseGuidanceResult guideSupercruise(const Ship& ship,
                                          const math::Vec3d& destPosKm,
                                          const math::Vec3d& destVelKmS,
                                          double dropRadiusKm,
                                          bool navAssistEnabled,
                                          bool dropRequested,
                                          bool interdicted,
                                          const SupercruiseParams& params) {
  SupercruiseGuidanceResult out{};
  out.recommendedMaxLinearAccelKmS2 = params.accelCapKmS2;
  out.recommendedMaxAngularAccelRadS2 = params.angularCapRadS2;

  // Default control flags.
  out.input.dampers = true;
  out.input.brake = false;
  out.input.boost = false;
  out.input.thrustLocal = {0, 0, 0};
  out.input.torqueLocal = {0, 0, 0};

  // --- Relative geometry ---
  const math::Vec3d rel = destPosKm - ship.positionKm();
  const double dist = rel.length();
  out.hud.distKm = dist;

  // Degenerate case: we're essentially at the destination.
  // In gameplay code this typically means "drop now" (not an emergency).
  if (dist < 1e-6) {
    out.hud.safeDropReady = true;
    // Respect interdiction gating (callers may prevent drops while interdicted).
    if (!interdicted) {
      out.dropNow = true;
      out.emergencyDrop = false;
    }
    // dirToDest defaults to forward.
    return out;
  }

  math::Vec3d dir = rel / dist;
  out.dirToDest = dir;

  const math::Vec3d vRel = ship.velocityKmS() - destVelKmS;
  const double closing = math::dot(vRel, dir);
  out.hud.closingKmS = closing;

  // Lateral relative speed (orthogonal to the line of sight).
  const math::Vec3d vLat = vRel - dir * closing;
  const double lateral = vLat.length();
  out.hud.lateralKmS = lateral;

  const double tta = (closing > 1e-3) ? (dist / closing) : 1e9;
  out.hud.ttaSec = tta;

  // Predicted miss distance under constant relative velocity.
  // rel(t) = rel - vRel * t, so closest approach occurs at:
  //   tCA = dot(rel, vRel) / |vRel|^2
  // (clamped to [0, +inf)).
  const double vRelSq = vRel.lengthSq();
  double missKm = dist;
  if (vRelSq > 1e-12) {
    const double tca = std::max(0.0, math::dot(rel, vRel) / vRelSq);
    const math::Vec3d relAt = rel - vRel * tca;
    missKm = relAt.length();
  }
  out.hud.missKm = missKm;

  // --- Safe-drop heuristic ---
  const double safeMin = std::max(0.0, params.safeTtaSec - params.safeWindowSlackSec);
  const double safeMax = params.safeTtaSec + params.safeWindowSlackSec;
  const bool lateralOk = (closing > params.minClosingKmS)
    ? (lateral <= std::max(0.0, closing) * std::max(0.0, params.maxLateralFrac))
    : true;
  const bool safeWindow = (dist <= dropRadiusKm)
    && (tta > safeMin) && (tta < safeMax)
    && (closing > params.minClosingKmS)
    && lateralOk;
  out.hud.safeDropReady = safeWindow;

  // --- Drop decisions ---
  // Manual drop always happens immediately (safe -> normal drop, otherwise emergency).
  if (!interdicted && dropRequested) {
    out.dropNow = true;
    out.emergencyDrop = !safeWindow;
    return out;
  }

  // Nav assist auto-drop when safe.
  if (!interdicted && navAssistEnabled && safeWindow) {
    out.dropNow = true;
    out.emergencyDrop = false;
    return out;
  }

  // --- Speed profile ---
  double desiredSpeed = 0.0;
  double speedLimit = params.maxSpeedKmS;

  // Precompute a corridor metric from the current relative velocity. This is used
  // by the optional corridor speed penalty, and is also exposed in the HUD.
  const double latFrac = (closing > 1e-6) ? (lateral / std::max(1e-6, closing)) : 0.0;

  if (navAssistEnabled) {
    // Classic heuristic wants v ~= dist / safeTtaSec, but that ignores finite acceleration.
    // We keep the heuristic as the baseline, then optionally cap it by a braking-distance
    // constraint that ensures we can slow down enough to reach the safe window.
    const double safeTta = std::max(0.1, params.safeTtaSec);
    const double base = std::clamp(dist / safeTta, params.assistSpeedMinKmS, params.maxSpeedKmS);

    desiredSpeed = base;

    if (params.useBrakingDistanceLimit) {
      const double dToDrop = std::max(0.0, dist - dropRadiusKm);
      const double vDrop = (dropRadiusKm > 0.0) ? (dropRadiusKm / safeTta) : 0.0;
      const double vAllow = std::sqrt(std::max(0.0, vDrop * vDrop + 2.0 * params.accelCapKmS2 * dToDrop));
      speedLimit = std::min(speedLimit, vAllow);
      desiredSpeed = std::min(desiredSpeed, speedLimit);
    }

    // Corridor-aware speed penalty: if we're sliding past the target (high
    // lateral vs closing), temporarily cap the speed. This gives the controller
    // more time to null lateral velocity before reaching the drop radius.
    double corridorFactor = 1.0;
    if (params.corridorSpeedPenalty && params.maxLateralFrac > 1e-6 && closing > params.minClosingKmS) {
      if (latFrac > params.maxLateralFrac * 1.02) {
        corridorFactor = params.maxLateralFrac / std::max(params.maxLateralFrac, latFrac);
        corridorFactor = std::clamp(corridorFactor, params.corridorPenaltyMinFactor, 1.0);
        corridorFactor = std::pow(corridorFactor, std::max(0.0, params.corridorPenaltyPower));
        speedLimit = std::clamp(speedLimit * corridorFactor, params.assistSpeedMinKmS, params.maxSpeedKmS);
        desiredSpeed = std::min(desiredSpeed, speedLimit);
      }
    }
    out.hud.corridorFactor = corridorFactor;
  } else {
    desiredSpeed = std::clamp(dist * params.manualSpeedDistGain, params.manualSpeedMinKmS, params.maxSpeedKmS);
  }

  out.hud.desiredSpeedKmS = desiredSpeed;
  out.hud.speedLimitKmS = speedLimit;

  // --- Intercept-course lead direction (translation only) ---
  //
  // For moving destinations, computing a lead direction helps prevent the
  // controller from constantly chasing the tail of the target's motion.
  math::Vec3d approachDir = dir;
  bool usedLead = false;
  double leadTimeSec = 0.0;

  if (params.interceptEnabled && desiredSpeed >= params.interceptMinSpeedKmS && destVelKmS.lengthSq() > 1e-12) {
    const double solveSpeed = std::max(1e-6,
      params.interceptUseMaxSpeedForSolve ? params.maxSpeedKmS : desiredSpeed);

    const auto lead = solveProjectileLead(ship.positionKm(),
                                         /*shooterVelKmS=*/{0, 0, 0},
                                         destPosKm,
                                         destVelKmS,
                                         solveSpeed,
                                         /*maxTimeSec=*/std::max(0.0, params.interceptMaxLeadTimeSec),
                                         /*minTimeSec=*/1.0e-3);

    if (lead && lead->aimDirWorld.lengthSq() > 1e-12) {
      approachDir = lead->aimDirWorld.normalized();
      usedLead = true;
      leadTimeSec = lead->tSec;

      // Clamp extreme lead angles.
      const double maxAng = math::degToRad(std::max(0.0, params.interceptMaxAngleDeg));
      if (maxAng > 1e-6) {
        const double ang = safeAcos(math::dot(approachDir, dir));
        if (ang > maxAng) {
          const double t = std::clamp(maxAng / std::max(1e-6, ang), 0.0, 1.0);
          approachDir = nlerpDir(dir, approachDir, t);
        }
      }
    }
  }

  out.hud.leadUsed = usedLead;
  out.hud.leadTimeSec = leadTimeSec;
  out.hud.leadAngleDeg = math::radToDeg(safeAcos(math::dot(approachDir, dir)));

  // --- Translational control ---
  const math::Vec3d desiredVel = destVelKmS + approachDir * desiredSpeed;
  const math::Vec3d dv = desiredVel - ship.velocityKmS();

  math::Vec3d thrustWorld{0, 0, 0};
  if (dv.lengthSq() > 1e-12) {
    const double tau = std::max(0.05, params.velTimeConstantSec);
    const double denom = std::max(1e-6, params.accelCapKmS2 * tau);
    thrustWorld = dv / denom;
    thrustWorld = clampMagnitude(thrustWorld, 1.0);
  }
  out.input.thrustLocal = ship.orientation().conjugate().rotate(thrustWorld);

  // --- Attitude control ---
  // When interdicted, callers may want the player to steer (so don't override torque).
  if (!interdicted) {
    const math::Vec3d desiredFwdLocal = ship.orientation().conjugate().rotate(dir);
    const double yawErr = std::atan2(desiredFwdLocal.x, desiredFwdLocal.z);
    const double pitchErr = -std::atan2(desiredFwdLocal.y, desiredFwdLocal.z);

    out.input.torqueLocal.x = std::clamp(pitchErr * params.faceGain, -1.0, 1.0);
    out.input.torqueLocal.y = std::clamp(yawErr * params.faceGain, -1.0, 1.0);
    out.input.torqueLocal.z = 0.0;
  }

  return out;
}

} // namespace stellar::sim
