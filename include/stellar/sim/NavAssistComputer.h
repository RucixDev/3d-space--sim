#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/FlightController.h"
#include "stellar/sim/Ship.h"

namespace stellar::sim {

// Lightweight flight assistance for normal-space maneuvering.
//
// This sits between raw manual thruster input and full autopilots like the
// DockingComputer. It exposes two useful modes:
//  - Approach: close to a target and hold at a configurable standoff distance.
//  - MatchVelocity: match the target's velocity while holding the initial
//                   separation as a soft constraint.
//
// Design constraints:
//  - Deterministic + headless (usable in tests / sandbox / server sims).
//  - Minimal state: the computer stores only the active mode and a desired
//    distance (captured on engage for MatchVelocity).
enum class NavAssistMode : core::u8 {
  Off = 0,
  Approach = 1,
  MatchVelocity = 2,
};

// Tunables for NavAssistComputer.
//
// Internally, the computer uses the existing FlightController helpers.
// The "slowDownRangeKm" knobs are converted into FlightControlParams::speedGain
// such that the controller reaches maxSpeed at roughly that remaining distance.
struct NavAssistParams {
  // --- Approach mode ---
  double approachMaxSpeedKmS{0.75};
  double approachSlowDownRangeKm{8000.0};
  double approachVelGain{1.9};
  bool approachAllowBoost{true};

  // --- Match velocity mode ---
  double matchMaxSpeedKmS{0.35};
  double matchSlowDownRangeKm{2500.0};
  double matchVelGain{2.2};
  bool matchAllowBoost{false};

  // --- Common ---
  double accelScale{1.0};
  bool dampers{true};

  // Attitude: face the target direction.
  double faceGain{1.7};
  double rollGain{1.6};

  // Use intercept-course guidance to reduce tail chasing.
  bool interceptEnabled{true};
  double interceptMaxLeadTimeSec{180.0};
  double interceptMinSpeedKmS{0.05};

  // Hold-distance clamp for MatchVelocity (captured at engage).
  double matchHoldDistMinKm{0.5};
  double matchHoldDistMaxKm{200000.0};

  // "Arrived" thresholds (used for UI / optional auto-disengage).
  double arriveDistEpsKm{50.0};
  double arriveRelSpeedEpsKmS{0.03};
  bool disengageOnArrive{false};
};

struct NavAssistResult {
  ShipInput input{};
  NavAssistMode mode{NavAssistMode::Off};

  // Debug/UI metrics.
  double distKm{0.0};
  double desiredDistKm{0.0};
  double relSpeedKmS{0.0};

  bool usedBoost{false};
  bool arrived{false};
};

class NavAssistComputer {
public:
  NavAssistComputer();

  const NavAssistParams& params() const { return params_; }
  void setParams(const NavAssistParams& p) { params_ = p; }

  bool active() const { return mode_ != NavAssistMode::Off; }
  NavAssistMode mode() const { return mode_; }

  double desiredDistKm() const { return desiredDistKm_; }

  void disengage();

  // Engage "approach" mode with a desired standoff distance.
  void engageApproach(double desiredDistKm);

  // Engage "match velocity" mode. By default, captures the current distance to
  // the target as the hold distance.
  //
  // desiredDistOverrideKm can be set to a non-negative value to force a specific
  // hold distance.
  void engageMatchVelocity(const Ship& ship,
                           const math::Vec3d& targetPosKm,
                           double desiredDistOverrideKm = -1.0);

  // Compute assisted inputs for this frame.
  //
  // dtSimSec is currently only used for future-proofing and is safe to pass 0.
  NavAssistResult update(const Ship& ship,
                         const math::Vec3d& targetPosKm,
                         const math::Vec3d& targetVelKmS,
                         double dtSimSec);

private:
  static double speedGainFromRange(double maxSpeedKmS, double slowDownRangeKm);
  FlightControlParams makeFlightParams(double desiredDistKm) const;
  AttitudeControlParams makeAttitudeParams() const;
  InterceptCourseParams makeInterceptParams() const;

  NavAssistParams params_{};
  NavAssistMode mode_{NavAssistMode::Off};
  double desiredDistKm_{0.0};
};

} // namespace stellar::sim
