#include "stellar/sim/ManeuverComputer.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static constexpr double kPi = 3.14159265358979323846;

static double clamp01(double v) { return std::clamp(v, 0.0, 1.0); }

static double clampD(double v, double lo, double hi) { return std::clamp(v, lo, hi); }

static math::Vec3d clampComponents(const math::Vec3d& v, double lo, double hi) {
  return {
    std::clamp(v.x, lo, hi),
    std::clamp(v.y, lo, hi),
    std::clamp(v.z, lo, hi),
  };
}

static double radToDeg(double r) { return r * (180.0 / kPi); }

// Minimal "face a direction" controller.
// Produces a normalized torque input in body-local axes.
//
// This is intentionally simpler than the full FlightController attitude code; the
// maneuver computer only needs to align forward, not manage roll precisely.
static math::Vec3d torqueToFaceForward(const Ship& ship,
                                       const math::Vec3d& desiredForwardWorld,
                                       double faceGain) {
  math::Vec3d fwdW = desiredForwardWorld;
  if (fwdW.lengthSq() <= 1e-12) return {0, 0, 0};
  fwdW = fwdW.normalized();

  math::Vec3d curFwdW = ship.forward();
  if (curFwdW.lengthSq() <= 1e-12) curFwdW = {0, 0, 1};
  curFwdW = curFwdW.normalized();

  math::Vec3d axisW = math::cross(curFwdW, fwdW);
  const double sinA = axisW.length();
  const double cosA = clampD(math::dot(curFwdW, fwdW), -1.0, 1.0);
  const double ang = std::atan2(sinA, cosA); // [0,pi]

  if (sinA <= 1e-9 || ang <= 1e-9) return {0, 0, 0};
  axisW = axisW * (1.0 / sinA);

  // Convert the world-space axis into ship body-local coordinates.
  const math::Vec3d axisLocal = ship.orientation().conjugate().rotate(axisW);

  // PD-lite: use angle as error; rely on ship's built-in angular dampers for D term.
  const math::Vec3d cmd = axisLocal * (ang * faceGain);

  return clampComponents(cmd, -1.0, 1.0);
}

void ManeuverComputer::disengage() {
  phase_ = ManeuverComputerPhase::Off;
  plan_ = {};
  burnDirWorld_ = {0, 0, 1};
  dvTotalKmS_ = 0.0;
  velAtBurnStartKmS_ = {0, 0, 0};
}

void ManeuverComputer::engage(const Ship& ship, const ManeuverPlan& plan) {
  plan_ = plan;
  dvTotalKmS_ = plan_.deltaVWorldKmS.length();
  burnDirWorld_ = (dvTotalKmS_ > 1e-12) ? (plan_.deltaVWorldKmS / dvTotalKmS_) : math::Vec3d{0,0,1};

  velAtBurnStartKmS_ = ship.velocityKmS();

  if (dvTotalKmS_ <= 1e-9) {
    phase_ = ManeuverComputerPhase::Complete;
  } else {
    phase_ = ManeuverComputerPhase::Orient;
  }
}

ManeuverComputerOutput ManeuverComputer::update(const Ship& ship,
                                               double nowTimeDays,
                                               double dtSimSec,
                                               const ManeuverComputerParams& params) {
  ManeuverComputerOutput out{};
  out.phase = phase_;

  if (phase_ == ManeuverComputerPhase::Off) {
    return out;
  }

  const double dvTotal = dvTotalKmS_;
  out.dvTotalKmS = dvTotal;

  // Time to node
  out.timeToNodeSec = (plan_.nodeTimeDays - nowTimeDays) * 86400.0;

  // Attitude guidance always runs while active (even during burn).
  {
    // Angle error
    math::Vec3d curFwd = ship.forward();
    if (curFwd.lengthSq() <= 1e-12) curFwd = {0,0,1};
    curFwd = curFwd.normalized();
    const math::Vec3d desFwd = burnDirWorld_;
    const double cosA = clampD(math::dot(curFwd, desFwd), -1.0, 1.0);
    const double ang = std::acos(cosA);
    out.alignmentErrorDeg = radToDeg(ang);

    out.input.torqueLocal = torqueToFaceForward(ship, desFwd, params.faceGain);
    out.input.dampers = true;
  }

  // Completion / early out
  if (phase_ == ManeuverComputerPhase::Complete) {
    out.finished = true;
    out.phase = phase_;
    return out;
  }
  if (phase_ == ManeuverComputerPhase::Aborted) {
    out.aborted = true;
    out.phase = phase_;
    return out;
  }

  if (dvTotal <= 1e-9) {
    phase_ = ManeuverComputerPhase::Complete;
    out.phase = phase_;
    out.finished = true;
    return out;
  }

  // Timing estimates
  const double accelCap = std::max(1e-9, params.allowBoost ? ship.maxLinearAccelBoostKmS2()
                                                           : ship.maxLinearAccelKmS2());
  out.burnDurationSec = dvTotal / accelCap;
  out.burnStartLeadSec = out.burnDurationSec * 0.5 + params.extraLeadTimeSec;

  // Abort if we missed the node by too much and never started burning.
  if (phase_ == ManeuverComputerPhase::Orient &&
      params.abortAfterMissedSec > 0.0 &&
      out.timeToNodeSec < -params.abortAfterMissedSec) {
    phase_ = ManeuverComputerPhase::Aborted;
    out.phase = phase_;
    out.aborted = true;
    return out;
  }

  const double tolRad = params.alignToleranceDeg * (kPi / 180.0);

  // Start burn when centered-burn start time is reached and we're aligned.
  if (phase_ == ManeuverComputerPhase::Orient) {
    const bool startTime = (out.timeToNodeSec <= out.burnStartLeadSec);
    const bool aligned = (out.alignmentErrorDeg * (kPi / 180.0) <= tolRad);
    if (startTime && aligned) {
      phase_ = ManeuverComputerPhase::Burn;
      velAtBurnStartKmS_ = ship.velocityKmS();
    }
  }

  if (phase_ == ManeuverComputerPhase::Burn) {
    // Estimate achieved delta-v by projecting actual velocity change along burn direction.
    const double dvSoFar = math::dot(ship.velocityKmS() - velAtBurnStartKmS_, burnDirWorld_);
    const double dvRemaining = std::max(0.0, dvTotal - dvSoFar);
    out.dvRemainingKmS = dvRemaining;

    if (dvRemaining <= params.dvToleranceKmS) {
      phase_ = ManeuverComputerPhase::Complete;
      out.phase = phase_;
      out.finished = true;
      return out;
    }

    // Throttle down on the final step to avoid overshooting the target delta-v.
    const double stepDvCap = accelCap * std::max(1e-6, dtSimSec);
    const double throttle = clamp01(dvRemaining / stepDvCap);

    out.input.thrustLocal = {0, 0, throttle};
    out.input.boost = params.allowBoost;

    if (params.disableDampersDuringBurn) {
      out.input.dampers = false;
    }
  }

  out.phase = phase_;
  return out;
}

} // namespace stellar::sim
