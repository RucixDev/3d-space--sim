#include "stellar/sim/NavAssistComputer.h"

#include "stellar/sim/ProximityField.h"

#include "stellar/core/Clamp.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

NavAssistComputer::NavAssistComputer() = default;

void NavAssistComputer::disengage() {
  mode_ = NavAssistMode::Off;
  desiredDistKm_ = 0.0;
  followDirInit_ = false;
  orbitInit_ = false;
  orbitSign_ = 1.0;
  avoidDirInit_ = false;
  bypassObstacleId_ = -1;
  bypassSideSign_ = 0;
}

void NavAssistComputer::engageApproach(double desiredDistKm) {
  mode_ = NavAssistMode::Approach;
  desiredDistKm_ = std::max(0.0, desiredDistKm);
  followDirInit_ = false;
  orbitInit_ = false;
  orbitSign_ = 1.0;
  avoidDirInit_ = false;
  bypassObstacleId_ = -1;
  bypassSideSign_ = 0;
}

void NavAssistComputer::engageMatchVelocity(const Ship& ship,
                                            const math::Vec3d& targetPosKm,
                                            double desiredDistOverrideKm) {
  mode_ = NavAssistMode::MatchVelocity;
  followDirInit_ = false;
  orbitInit_ = false;
  orbitSign_ = 1.0;
  avoidDirInit_ = false;
  bypassObstacleId_ = -1;
  bypassSideSign_ = 0;

  double d = desiredDistOverrideKm;
  if (d < 0.0) {
    d = (targetPosKm - ship.positionKm()).length();
  }

  d = std::clamp(d, params_.matchHoldDistMinKm, params_.matchHoldDistMaxKm);
  desiredDistKm_ = std::max(0.0, d);
}

void NavAssistComputer::engageFollow(const Ship& ship,
                                    const math::Vec3d& targetPosKm,
                                    const math::Vec3d& targetVelKmS,
                                    double desiredDistOverrideKm) {
  mode_ = NavAssistMode::Follow;
  orbitInit_ = false;
  orbitSign_ = 1.0;
  avoidDirInit_ = false;
  bypassObstacleId_ = -1;
  bypassSideSign_ = 0;

  double d = desiredDistOverrideKm;
  if (d < 0.0) {
    d = (targetPosKm - ship.positionKm()).length();
  }
  d = std::clamp(d, params_.matchHoldDistMinKm, params_.matchHoldDistMaxKm);
  desiredDistKm_ = std::max(0.0, d);

  // Initialize follow direction from target velocity if possible, otherwise
  // from line-of-sight.
  const double sp = targetVelKmS.length();
  if (sp >= std::max(0.0, params_.followMinTargetSpeedKmS)) {
    followDirUnit_ = (targetVelKmS / sp);
    followDirInit_ = true;
  } else {
    const auto rel = (targetPosKm - ship.positionKm());
    const auto dir = rel.normalized();
    if (dir.lengthSq() > 1e-12) {
      followDirUnit_ = dir;
      followDirInit_ = true;
    } else {
      followDirUnit_ = {0,0,1};
      followDirInit_ = false;
    }
  }
}

void NavAssistComputer::engageOrbit(const Ship& ship,
                                     const math::Vec3d& targetPosKm,
                                     const math::Vec3d& targetVelKmS,
                                     double desiredDistOverrideKm) {
  mode_ = NavAssistMode::Orbit;
  followDirInit_ = false;
  avoidDirInit_ = false;
  bypassObstacleId_ = -1;
  bypassSideSign_ = 0;

  double d = desiredDistOverrideKm;
  if (d < 0.0) {
    d = (targetPosKm - ship.positionKm()).length();
  }

  d = std::clamp(d, params_.matchHoldDistMinKm, params_.matchHoldDistMaxKm);
  desiredDistKm_ = std::max(0.0, d);

  // Choose a stable orbit plane normal.
  // Prefer the ship's current up vector, but fall back to a world axis.
  const math::Vec3d worldUp{0, 1, 0};
  const math::Vec3d worldFwd{0, 0, 1};

  math::Vec3d radial = (ship.positionKm() - targetPosKm).normalized();
  if (radial.lengthSq() <= 1e-12) {
    radial = ship.forward().normalized();
    if (radial.lengthSq() <= 1e-12) radial = math::Vec3d{1, 0, 0};
  }

  math::Vec3d n = ship.up().normalized();
  if (n.lengthSq() <= 1e-12) n = worldUp;

  // Avoid a near-parallel normal (degenerate tangent).
  if (std::abs(math::dot(n, radial)) > 0.97) {
    n = (std::abs(math::dot(worldUp, radial)) < 0.97) ? worldUp : worldFwd;
  }

  orbitNormalUnit_ = n.normalized();
  if (orbitNormalUnit_.lengthSq() <= 1e-12) orbitNormalUnit_ = worldUp;
  orbitInit_ = true;

  // Decide orbit direction from current tangential motion (relative to target).
  const math::Vec3d relVel = ship.velocityKmS() - targetVelKmS;
  math::Vec3d tangent = math::cross(orbitNormalUnit_, radial);
  if (tangent.lengthSq() <= 1e-12) tangent = math::cross(worldUp, radial);
  tangent = tangent.normalized();
  const double proj = math::dot(relVel, tangent);
  orbitSign_ = (proj < 0.0) ? -1.0 : 1.0;
}

double NavAssistComputer::speedGainFromRange(double maxSpeedKmS, double slowDownRangeKm) {
  maxSpeedKmS = std::max(0.0, maxSpeedKmS);
  slowDownRangeKm = std::max(1e-6, slowDownRangeKm);
  return maxSpeedKmS / slowDownRangeKm;
}

FlightControlParams NavAssistComputer::makeFlightParams(double desiredDistKm) const {
  FlightControlParams fp{};
  fp.desiredDistKm = std::max(0.0, desiredDistKm);
  fp.accelScale = std::max(0.0, params_.accelScale);
  fp.dampers = params_.dampers;

  if (mode_ == NavAssistMode::MatchVelocity) {
    fp.maxSpeedKmS = std::max(0.0, params_.matchMaxSpeedKmS);
    fp.speedGain = speedGainFromRange(fp.maxSpeedKmS, params_.matchSlowDownRangeKm);
    fp.velGain = std::max(0.0, params_.matchVelGain);
    fp.allowBoost = params_.matchAllowBoost;
  } else if (mode_ == NavAssistMode::Follow) {
    fp.maxSpeedKmS = std::max(0.0, params_.followMaxSpeedKmS);
    fp.speedGain = speedGainFromRange(fp.maxSpeedKmS, params_.followSlowDownRangeKm);
    fp.velGain = std::max(0.0, params_.followVelGain);
    fp.allowBoost = params_.followAllowBoost;
  } else if (mode_ == NavAssistMode::Orbit) {
    fp.maxSpeedKmS = std::max(0.0, params_.orbitMaxSpeedKmS);
    fp.speedGain = speedGainFromRange(fp.maxSpeedKmS, params_.orbitSlowDownRangeKm);
    fp.velGain = std::max(0.0, params_.orbitVelGain);
    fp.allowBoost = params_.orbitAllowBoost;
  } else {
    // Default to Approach tuning for all other modes.
    fp.maxSpeedKmS = std::max(0.0, params_.approachMaxSpeedKmS);
    fp.speedGain = speedGainFromRange(fp.maxSpeedKmS, params_.approachSlowDownRangeKm);
    fp.velGain = std::max(0.0, params_.approachVelGain);
    fp.allowBoost = params_.approachAllowBoost;
  }

  return fp;
}

AttitudeControlParams NavAssistComputer::makeAttitudeParams() const {
  AttitudeControlParams ap{};
  ap.faceGain = std::max(0.0, params_.faceGain);
  ap.rollGain = std::max(0.0, params_.rollGain);
  ap.alignUp = false;
  return ap;
}

InterceptCourseParams NavAssistComputer::makeInterceptParams() const {
  InterceptCourseParams ic{};
  ic.enabled = params_.interceptEnabled;
  ic.maxLeadTimeSec = std::max(0.0, params_.interceptMaxLeadTimeSec);
  ic.minSpeedKmS = std::max(0.0, params_.interceptMinSpeedKmS);
  ic.useMaxSpeedForSolve = true;
  return ic;
}

NavAssistResult NavAssistComputer::update(const Ship& ship,
                                         const math::Vec3d& targetPosKm,
                                         const math::Vec3d& targetVelKmS,
                                         double dtSimSec) {
  return update(ship, targetPosKm, targetVelKmS, dtSimSec, nullptr, -1);
}

NavAssistResult NavAssistComputer::update(const Ship& ship,
                                         const math::Vec3d& targetPosKm,
                                         const math::Vec3d& targetVelKmS,
                                         double dtSimSec,
                                         const ProximityFieldKm* obstaclesKm,
                                         int ignoreObstacleId) {
  NavAssistResult res{};
  res.mode = mode_;
  res.desiredDistKm = desiredDistKm_;

  if (mode_ == NavAssistMode::Off) {
    return res;
  }

  // Optional: compute an avoidance-biased translation target while keeping the
  // facing direction locked onto the true target.
  auto computeAvoidedTarget = [&](const math::Vec3d& goalPosKm,
                                  double desiredSpeedKmS,
                                  const math::Vec3d& facingDirUnit,
                                  const math::Vec3d& fallbackDirUnit,
                                  int ignoreId) -> math::Vec3d {
    (void)facingDirUnit;
    math::Vec3d desiredDir = (goalPosKm - ship.positionKm()).normalized();
    if (desiredDir.lengthSq() <= 1e-12) {
      desiredDir = fallbackDirUnit.lengthSq() > 1e-12 ? fallbackDirUnit : math::Vec3d{0,0,1};
    }

    const bool useAvoid = (obstaclesKm != nullptr) && params_.avoidanceEnabled;


    // If avoidance isn't active, clear any bypass preference state.
    if (!useAvoid) {
      bypassObstacleId_ = -1;
      bypassSideSign_ = 0;
    } else if (params_.avoidanceTangentBypassEnabled) {
      TangentBypassParamsKm bp{};
      bp.enabled = params_.avoidanceTangentBypassEnabled;
      bp.shipRadiusKm = std::max(0.0, params_.avoidanceShipRadiusKm);
      bp.padKm = std::max(0.0, params_.avoidancePadKm);
      bp.nearMissExtraKm = std::max(0.0, params_.avoidanceNearMissExtraKm);
      bp.aheadClearanceMult = std::max(0.0, params_.avoidanceBypassAheadClearanceMult);
      bp.aheadExtraKm = std::max(0.0, params_.avoidanceBypassAheadExtraKm);
      bp.maxWaypointDistKm = std::max(0.0, params_.avoidanceBypassMaxWaypointDistKm);

      const int preferredSide = (bypassObstacleId_ >= 0) ? bypassSideSign_ : 0;
      const auto br = planTangentBypassWaypoint(*obstaclesKm,
                                                ship.positionKm(),
                                                goalPosKm,
                                                bp,
                                                ignoreId,
                                                preferredSide,
                                                /*maxSphereTests=*/256);

      if (br.used) {
        bypassObstacleId_ = br.obstacleId;
        bypassSideSign_ = br.sideSign;

        res.avoidanceBypassUsed = true;
        res.avoidanceBypassWaypointKm = br.waypointKm;
        res.avoidanceBypassObstacleId = br.obstacleId;
        res.avoidanceBypassSideSign = br.sideSign;

        // Reset the smoothed avoidance direction to the bypass direction.
        avoidDirUnit_ = (br.waypointKm - ship.positionKm()).normalized();
        avoidDirInit_ = (avoidDirUnit_.lengthSq() > 1e-12);

        // The bypass waypoint is already validated against the inflated obstacle set,
        // so steer directly toward it this frame.
        return br.waypointKm;
      }

      // If a bypass couldn't be found, fall back to potential-field steering.
      bypassObstacleId_ = -1;
      bypassSideSign_ = 0;
    } else {
      bypassObstacleId_ = -1;
      bypassSideSign_ = 0;
    }


    math::Vec3d safeDir = desiredDir;
    if (useAvoid) {
      AvoidanceParamsKm ap{};
      ap.enabled = params_.avoidanceEnabled;
      ap.shipRadiusKm = std::max(0.0, params_.avoidanceShipRadiusKm);
      ap.padKm = std::max(0.0, params_.avoidancePadKm);
      ap.lookaheadTimeSec = std::max(0.0, params_.avoidanceLookaheadTimeSec);
      ap.lookaheadBaseKm = std::max(0.0, params_.avoidanceLookaheadBaseKm);
      ap.minSpeedForLookaheadKmS = std::max(0.0, params_.avoidanceMinSpeedForLookaheadKmS);
      ap.nearMissExtraKm = std::max(0.0, params_.avoidanceNearMissExtraKm);
      ap.strength = std::max(0.0, params_.avoidanceStrength);
      ap.maxSteerDeg = std::max(0.0, params_.avoidanceMaxSteerDeg);

      const auto ar = steerAvoidObstacles(*obstaclesKm,
                                          ship.positionKm(),
                                          ship.velocityKmS(),
                                          desiredDir,
                                          desiredSpeedKmS,
                                          ap,
                                          ignoreId);
      safeDir = ar.safeDirUnit;

      res.avoidanceSteering = ar.steering;
      res.avoidanceThreatId = ar.threatId;
      res.avoidanceThreatClearanceKm = ar.threatClearanceKm;
      res.avoidanceLookaheadKm = ar.lookaheadKm;
    }

    // Smooth the safe direction to avoid jitter when multiple obstacles compete.
    // If avoidance isn't active, track the true desired direction with no lag.
    if (!useAvoid) {
      avoidDirUnit_ = desiredDir;
      avoidDirInit_ = true;
    } else {
      const double tau = params_.avoidanceDirBlendTauSec;
      if (tau <= 0.0 || dtSimSec <= 0.0) {
        avoidDirUnit_ = safeDir.normalized();
        avoidDirInit_ = (avoidDirUnit_.lengthSq() > 1e-12);
      } else {
        const double alpha = 1.0 - std::exp(-dtSimSec / std::max(1e-6, tau));
        if (!avoidDirInit_) {
          avoidDirUnit_ = safeDir.normalized();
          avoidDirInit_ = (avoidDirUnit_.lengthSq() > 1e-12);
        } else {
          const auto blended = (avoidDirUnit_ * (1.0 - alpha)) + (safeDir * alpha);
          const auto n = blended.normalized();
          if (n.lengthSq() > 1e-12) avoidDirUnit_ = n;
        }
      }
    }

    const math::Vec3d dirFinal = (avoidDirInit_ && avoidDirUnit_.lengthSq() > 1e-12) ? avoidDirUnit_ : desiredDir;
    const double dist = (goalPosKm - ship.positionKm()).length();
    if (dist <= 1e-9) return goalPosKm;
    return ship.positionKm() + dirFinal * dist;
  };

  // Orbit mode: hold a standoff distance while moving tangentially around the target.
  // This is useful for scanning patterns, dogfighting, or keeping a target in view
  // while maintaining lateral motion.
  if (mode_ == NavAssistMode::Orbit) {
    const auto ap = makeAttitudeParams();
    const auto ic = makeInterceptParams();

    // Radial direction from target -> ship.
    const auto rel = (ship.positionKm() - targetPosKm);
    math::Vec3d radial = rel.normalized();
    if (radial.lengthSq() <= 1e-12) {
      radial = ship.forward().normalized();
      if (radial.lengthSq() <= 1e-12) radial = math::Vec3d{1, 0, 0};
    }

    // Ensure we have a sane orbit normal that is not near-parallel to the radial.
    const math::Vec3d worldUp{0, 1, 0};
    const math::Vec3d worldFwd{0, 0, 1};

    math::Vec3d normal = orbitNormalUnit_;
    if (!orbitInit_ || normal.lengthSq() <= 1e-12) {
      normal = ship.up().normalized();
      if (normal.lengthSq() <= 1e-12) normal = worldUp;
    } else {
      normal = normal.normalized();
    }

    if (std::abs(math::dot(normal, radial)) > 0.97) {
      normal = (std::abs(math::dot(worldUp, radial)) < 0.97) ? worldUp : worldFwd;
    }

    normal = normal.normalized();
    if (normal.lengthSq() <= 1e-12) normal = worldUp;
    orbitNormalUnit_ = normal;
    orbitInit_ = true;

    // Tangential direction (right-handed around `normal`).
    math::Vec3d tangent = math::cross(normal, radial);
    if (tangent.lengthSq() <= 1e-12) {
      const math::Vec3d alt = (std::abs(math::dot(worldUp, radial)) < 0.97) ? worldUp : worldFwd;
      tangent = math::cross(alt, radial);
    }

    math::Vec3d tangentDir = tangent.normalized();
    if (tangentDir.lengthSq() <= 1e-12) tangentDir = worldFwd;
    tangentDir = tangentDir * ((orbitSign_ < 0.0) ? -1.0 : 1.0);

    const double maxSpeed = std::max(0.0, params_.orbitMaxSpeedKmS);
    const double tangentialSpeed = std::clamp(std::max(0.0, params_.orbitTangentialSpeedKmS), 0.0, maxSpeed);

    const double leadTime = std::max(0.0, params_.orbitLeadTimeSec);
    const double leadCap = std::max(0.0, params_.orbitLeadMaxFrac) * std::max(0.0, desiredDistKm_);
    const double leadDist = std::clamp(tangentialSpeed * leadTime, 0.0, leadCap);

    const auto ringPointKm = targetPosKm + radial * desiredDistKm_;
    const auto orbitPointKm = ringPointKm + tangentDir * leadDist;
    const auto orbitVelKmS = targetVelKmS + tangentDir * tangentialSpeed;

    math::Vec3d desiredForward{0, 0, 1};
    if (params_.orbitFaceTarget) {
      const auto toTarget = (targetPosKm - ship.positionKm());
      desiredForward = toTarget.normalized();
      if (desiredForward.lengthSq() <= 1e-12) desiredForward = (-radial).normalized();
    } else {
      desiredForward = tangentDir;
    }

    if (desiredForward.lengthSq() <= 1e-12) desiredForward = worldFwd;

    // Optionally bias the translation target away from obstacles.
    const auto orbitPointAvoidKm = computeAvoidedTarget(orbitPointKm,
                                                        /*desiredSpeedKmS=*/maxSpeed,
                                                        desiredForward,
                                                        tangentDir,
                                                        ignoreObstacleId);

    auto fp = makeFlightParams(/*desiredDistKm=*/0.0);
    fp.desiredDistKm = 0.0;

    const auto out = approachTargetIntercept(ship, orbitPointAvoidKm, orbitVelKmS, fp, ap, desiredForward, ic);
    res.input = out.input;
    res.usedBoost = out.usedBoost;

    res.distKm = (ship.positionKm() - targetPosKm).length();
    res.relSpeedKmS = (ship.velocityKmS() - orbitVelKmS).length();

    const double distErr = std::abs(res.distKm - desiredDistKm_);
    const bool closeEnough = distErr <= std::max(0.0, params_.arriveDistEpsKm);
    const bool speedEnough = res.relSpeedKmS <= std::max(0.0, params_.arriveRelSpeedEpsKmS);
    res.arrived = closeEnough && speedEnough;

    if (params_.disengageOnArrive && res.arrived) {
      disengage();
      res.mode = mode_;
    }

    return res;
  }

  // Follow mode = chase a *moving point behind the target* (formation assist).
  // We still measure arrival against the true target distance and relative speed
  // so UI can show meaningful numbers.
  if (mode_ == NavAssistMode::Follow) {
    const auto ap = makeAttitudeParams();
    const auto ic = makeInterceptParams();

    // Choose a travel direction to place the follow point behind the target.
    math::Vec3d dir{0,0,0};
    const double sp = targetVelKmS.length();
    if (sp >= std::max(0.0, params_.followMinTargetSpeedKmS)) {
      dir = targetVelKmS / sp;
    } else {
      const auto rel = (targetPosKm - ship.positionKm());
      dir = rel.normalized();
    }

    if (dir.lengthSq() < 1e-12) {
      dir = followDirUnit_.lengthSq() > 1e-12 ? followDirUnit_ : math::Vec3d{0,0,1};
    }

    // Exponential smoothing to avoid jitter when target velocity is noisy.
    const double tau = std::max(1e-3, params_.followDirBlendTauSec);
    const double alpha = (dtSimSec > 0.0) ? (1.0 - std::exp(-dtSimSec / tau)) : 1.0;
    if (!followDirInit_) {
      followDirUnit_ = dir.normalized();
      followDirInit_ = true;
    } else {
      const auto blended = (followDirUnit_ * (1.0 - alpha)) + (dir * alpha);
      const auto n = blended.normalized();
      if (n.lengthSq() > 1e-12) followDirUnit_ = n;
    }

    const auto followPointKm = targetPosKm - followDirUnit_ * desiredDistKm_;

    // For formation flying it's nicer to *keep facing the target* while translating
    // toward a point behind it.
    math::Vec3d desiredForward{0,0,1};
    if (params_.followFaceTarget) {
      const auto toTarget = (targetPosKm - ship.positionKm());
      if (toTarget.lengthSq() > 1e-12) desiredForward = toTarget.normalized();
      else if (followDirUnit_.lengthSq() > 1e-12) desiredForward = followDirUnit_;
    } else {
      desiredForward = (followDirUnit_.lengthSq() > 1e-12) ? followDirUnit_ : math::Vec3d{0,0,1};
    }

    // Optionally bias the *translation* target away from obstacles.
    const auto followPointAvoidKm = computeAvoidedTarget(followPointKm,
                                                         /*desiredSpeedKmS=*/std::max(0.0, params_.followMaxSpeedKmS),
                                                         desiredForward,
                                                         followDirUnit_,
                                                         ignoreObstacleId);

    // Chase the follow point, but measure errors vs the real target.
    auto fp = makeFlightParams(/*desiredDistKm=*/0.0);
    fp.desiredDistKm = 0.0;

    const auto out = approachTargetIntercept(ship, followPointAvoidKm, targetVelKmS, fp, ap, desiredForward, ic);
    res.input = out.input;
    res.usedBoost = out.usedBoost;

    res.distKm = (ship.positionKm() - targetPosKm).length();
    res.relSpeedKmS = (ship.velocityKmS() - targetVelKmS).length();

    const double distErr = std::abs(res.distKm - desiredDistKm_);
    const bool closeEnough = distErr <= std::max(0.0, params_.arriveDistEpsKm);
    const bool speedEnough = res.relSpeedKmS <= std::max(0.0, params_.arriveRelSpeedEpsKmS);
    res.arrived = closeEnough && speedEnough;

    if (params_.disengageOnArrive && res.arrived) {
      disengage();
      res.mode = mode_;
    }

    return res;
  }

  // Compute guidance using the shared FlightController helpers.
  const auto fp = makeFlightParams(desiredDistKm_);
  const auto ap = makeAttitudeParams();
  const auto ic = makeInterceptParams();

  math::Vec3d desiredForward = (targetPosKm - ship.positionKm()).normalized();
  if (desiredForward.lengthSq() <= 1e-12) {
    desiredForward = ship.forward().normalized();
    if (desiredForward.lengthSq() <= 1e-12) desiredForward = math::Vec3d{0,0,1};
  }

  const auto targetAvoidKm = computeAvoidedTarget(targetPosKm,
                                                  /*desiredSpeedKmS=*/fp.maxSpeedKmS,
                                                  desiredForward,
                                                  desiredForward,
                                                  ignoreObstacleId);

  const auto out = approachTargetIntercept(ship, targetAvoidKm, targetVelKmS, fp, ap, desiredForward, ic);
  res.input = out.input;
  res.usedBoost = out.usedBoost;

  // Metrics for UI / tests.
  res.distKm = out.distKm;
  res.relSpeedKmS = (ship.velocityKmS() - targetVelKmS).length();

  const double distErr = std::abs(res.distKm - desiredDistKm_);
  const bool closeEnough = distErr <= std::max(0.0, params_.arriveDistEpsKm);
  const bool speedEnough = res.relSpeedKmS <= std::max(0.0, params_.arriveRelSpeedEpsKmS);
  res.arrived = closeEnough && speedEnough;

  if (params_.disengageOnArrive && res.arrived) {
    disengage();
    res.mode = mode_;
  }

  return res;
}

} // namespace stellar::sim
