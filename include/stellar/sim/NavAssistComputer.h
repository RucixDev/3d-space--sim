#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/FlightController.h"
#include "stellar/sim/Ship.h"

namespace stellar::sim {

class ProximityFieldKm;

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
  // Formation-like guidance: keep an offset behind the target based on its
  // current travel direction while matching velocity.
  Follow = 3,

  // Strafe/orbit guidance: maintain a standoff distance while moving
  // tangentially around the target in a stable plane.
  Orbit = 4,
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

  // --- Follow mode (formation) ---
  // Keep a point behind the target (based on target velocity direction) at a
  // configurable distance while matching velocity.
  double followMaxSpeedKmS{0.55};
  double followSlowDownRangeKm{5000.0};
  double followVelGain{2.0};
  bool followAllowBoost{false};

  // When the target isn't moving much, fall back to the line-of-sight direction.
  double followMinTargetSpeedKmS{0.02};

  // Smoothing for the travel direction used to compute the follow point.
  // Higher values = smoother but more lag.
  double followDirBlendTauSec{0.75};

  // If true, keep facing the *target* while translating toward the follow point.
  // If false, face the target's travel direction.
  bool followFaceTarget{true};

  // --- Orbit mode (strafe) ---
  // Maintain a ring around the target while moving tangentially in a stable plane.
  double orbitMaxSpeedKmS{0.55};
  double orbitSlowDownRangeKm{4500.0};
  double orbitVelGain{2.0};
  bool orbitAllowBoost{false};

  // Tangential speed around the target (added to target velocity).
  double orbitTangentialSpeedKmS{0.22};

  // Place the orbit aim point ahead along the tangent direction by:
  //   leadDist = orbitTangentialSpeedKmS * orbitLeadTimeSec,
  // clamped to [0, orbitLeadMaxFrac * desiredDistKm].
  double orbitLeadTimeSec{4.0};
  double orbitLeadMaxFrac{0.65};

  // If true, keep facing the target while orbiting; otherwise face the travel direction.
  bool orbitFaceTarget{true};

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

  // --- Optional local obstacle avoidance ---
  // When enabled (and an obstacle field is provided), nav assist will bias the
  // translation guidance away from nearby obstacles while keeping the ship's
  // facing direction focused on the true target.
  bool avoidanceEnabled{true};

  // Obstacles are inflated by (shipRadius + pad) for clearance checks.
  double avoidanceShipRadiusKm{0.03};
  double avoidancePadKm{0.05};

  // Lookahead distance = base + speed * lookaheadTimeSec.
  double avoidanceLookaheadTimeSec{8.0};
  double avoidanceLookaheadBaseKm{0.0};
  double avoidanceMinSpeedForLookaheadKmS{0.05};

  // Begin steering when clearance drops below this threshold.
  double avoidanceNearMissExtraKm{0.35};

  // Steering strength and maximum deviation from the desired direction.
  double avoidanceStrength{1.25};
  double avoidanceMaxSteerDeg{35.0};

  // Smooth the avoidance direction to prevent jitter.
  double avoidanceDirBlendTauSec{0.35};


  // If true, and the straight-line segment to the current translation goal is blocked,
  // compute a deterministic “tangent bypass” waypoint around the nearest obstacle
  // before falling back to potential-field steering.
  bool avoidanceTangentBypassEnabled{true};

  // Waypoint placement tuning for tangent bypass:
  //   tBeyond = tExit + clearanceKm * avoidanceBypassAheadClearanceMult + avoidanceBypassAheadExtraKm
  double avoidanceBypassAheadClearanceMult{1.0};
  double avoidanceBypassAheadExtraKm{0.0};

  // Optional clamp on distance from ship -> bypass waypoint (0 = no clamp).
  double avoidanceBypassMaxWaypointDistKm{0.0};

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

  // Optional avoidance debug.
  bool avoidanceSteering{false};
  int avoidanceThreatId{-1};
  double avoidanceThreatClearanceKm{0.0};
  double avoidanceLookaheadKm{0.0};


  bool avoidanceBypassUsed{false};
  math::Vec3d avoidanceBypassWaypointKm{0, 0, 0};
  int avoidanceBypassObstacleId{-1};
  int avoidanceBypassSideSign{0};

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

  // Engage "follow" mode: hold behind the target based on its travel direction.
  //
  // desiredDistOverrideKm can be set to a non-negative value to force a specific
  // follow distance.
  void engageFollow(const Ship& ship,
                    const math::Vec3d& targetPosKm,
                    const math::Vec3d& targetVelKmS,
                    double desiredDistOverrideKm = -1.0);

  // Engage "orbit" mode: maintain a standoff distance while translating tangentially
  // around the target in a stable plane.
  //
  // desiredDistOverrideKm can be set to a non-negative value to force a specific
  // orbit radius (otherwise uses current separation).
  void engageOrbit(const Ship& ship,
                   const math::Vec3d& targetPosKm,
                   const math::Vec3d& targetVelKmS,
                   double desiredDistOverrideKm = -1.0);

  // Compute assisted inputs for this frame.
  //
  // dtSimSec is currently only used for future-proofing and is safe to pass 0.
  NavAssistResult update(const Ship& ship,
                         const math::Vec3d& targetPosKm,
                         const math::Vec3d& targetVelKmS,
                         double dtSimSec);

  // Like update(), but optionally biases translation away from obstacles.
  //
  // ignoreObstacleId can be used to exclude an obstacle (e.g. the target asteroid)
  // from avoidance steering.
  NavAssistResult update(const Ship& ship,
                         const math::Vec3d& targetPosKm,
                         const math::Vec3d& targetVelKmS,
                         double dtSimSec,
                         const ProximityFieldKm* obstaclesKm,
                         int ignoreObstacleId = -1);

private:
  static double speedGainFromRange(double maxSpeedKmS, double slowDownRangeKm);
  FlightControlParams makeFlightParams(double desiredDistKm) const;
  AttitudeControlParams makeAttitudeParams() const;
  InterceptCourseParams makeInterceptParams() const;

  NavAssistParams params_{};
  NavAssistMode mode_{NavAssistMode::Off};
  double desiredDistKm_{0.0};

  // Internal state for Follow mode.
  math::Vec3d followDirUnit_{0, 0, 1};
  bool followDirInit_{false};

  // Internal state for Orbit mode.
  math::Vec3d orbitNormalUnit_{0, 1, 0};
  bool orbitInit_{false};
  double orbitSign_{1.0};

  // Internal state for avoidance direction smoothing.
  math::Vec3d avoidDirUnit_{0, 0, 1};
  bool avoidDirInit_{false};


  // Internal state for tangent-bypass waypoint selection.
  int bypassObstacleId_{-1};
  int bypassSideSign_{0};

};

} // namespace stellar::sim
