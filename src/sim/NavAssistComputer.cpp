#include "stellar/sim/NavAssistComputer.h"

#include "stellar/core/Clamp.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

NavAssistComputer::NavAssistComputer() = default;

void NavAssistComputer::disengage() {
  mode_ = NavAssistMode::Off;
  desiredDistKm_ = 0.0;
}

void NavAssistComputer::engageApproach(double desiredDistKm) {
  mode_ = NavAssistMode::Approach;
  desiredDistKm_ = std::max(0.0, desiredDistKm);
}

void NavAssistComputer::engageMatchVelocity(const Ship& ship,
                                            const math::Vec3d& targetPosKm,
                                            double desiredDistOverrideKm) {
  mode_ = NavAssistMode::MatchVelocity;

  double d = desiredDistOverrideKm;
  if (d < 0.0) {
    d = (targetPosKm - ship.positionKm()).length();
  }

  d = std::clamp(d, params_.matchHoldDistMinKm, params_.matchHoldDistMaxKm);
  desiredDistKm_ = std::max(0.0, d);
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
                                         double /*dtSimSec*/) {
  NavAssistResult res{};
  res.mode = mode_;
  res.desiredDistKm = desiredDistKm_;

  if (mode_ == NavAssistMode::Off) {
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
