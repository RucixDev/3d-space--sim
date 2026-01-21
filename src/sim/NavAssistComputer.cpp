#include "stellar/sim/NavAssistComputer.h"

#include "stellar/core/Clamp.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

NavAssistComputer::NavAssistComputer() = default;

void NavAssistComputer::disengage() {
  mode_ = NavAssistMode::Off;
  desiredDistKm_ = 0.0;
  followDirInit_ = false;
}

void NavAssistComputer::engageApproach(double desiredDistKm) {
  mode_ = NavAssistMode::Approach;
  desiredDistKm_ = std::max(0.0, desiredDistKm);
  followDirInit_ = false;
}

void NavAssistComputer::engageMatchVelocity(const Ship& ship,
                                            const math::Vec3d& targetPosKm,
                                            double desiredDistOverrideKm) {
  mode_ = NavAssistMode::MatchVelocity;
  followDirInit_ = false;

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
  NavAssistResult res{};
  res.mode = mode_;
  res.desiredDistKm = desiredDistKm_;

  if (mode_ == NavAssistMode::Off) {
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

    // Chase the follow point, but measure errors vs the real target.
    auto fp = makeFlightParams(/*desiredDistKm=*/0.0);
    fp.desiredDistKm = 0.0;

    const auto out = approachTargetIntercept(ship, followPointKm, targetVelKmS, fp, ap, desiredForward, ic);
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

  const auto out = chaseTargetIntercept(ship, targetPosKm, targetVelKmS, fp, ap, ic);
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
