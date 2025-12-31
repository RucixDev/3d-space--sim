#include "stellar/sim/SupercruiseComputer.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static math::Vec3d clampMagnitude(const math::Vec3d& v, double maxLen) {
  const double len = v.length();
  if (len <= maxLen || len <= 1e-12) return v;
  return v * (maxLen / len);
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
    out.dropNow = true;
    out.emergencyDrop = false;
    // dirToDest defaults to forward.
    return out;
  }

  math::Vec3d dir = rel / dist;
  out.dirToDest = dir;

  const math::Vec3d vRel = ship.velocityKmS() - destVelKmS;
  const double closing = math::dot(vRel, dir);
  out.hud.closingKmS = closing;

  const double tta = (closing > 1e-3) ? (dist / closing) : 1e9;
  out.hud.ttaSec = tta;

  // --- Safe-drop heuristic ---
  const double safeMin = std::max(0.0, params.safeTtaSec - params.safeWindowSlackSec);
  const double safeMax = params.safeTtaSec + params.safeWindowSlackSec;
  const bool safeWindow = (dist < dropRadiusKm) && (tta > safeMin) && (tta < safeMax) && (closing > params.minClosingKmS);
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
  } else {
    desiredSpeed = std::clamp(dist * params.manualSpeedDistGain, params.manualSpeedMinKmS, params.maxSpeedKmS);
  }

  out.hud.desiredSpeedKmS = desiredSpeed;
  out.hud.speedLimitKmS = speedLimit;

  // --- Translational control ---
  const math::Vec3d desiredVel = destVelKmS + dir * desiredSpeed;
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
