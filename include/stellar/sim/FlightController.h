#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/Ship.h"

namespace stellar::sim {

// A tiny, reusable flight-control helper that turns high-level goals
// ("approach this moving point", "face this direction") into ShipInput.
//
// The prototype previously had several ad-hoc controllers (player autopilot,
// NPC chase logic). This module centralizes the common math so behavior is
// consistent, tuneable, and unit-testable.

struct FlightControlParams {
  // Desired max closing speed toward the target (km/s).
  double maxSpeedKmS{0.25};

  // Gain converting distance-to-target into desired closing speed:
  //   desiredSpeed = min(maxSpeed, speedGain * max(dist - desiredDist, 0)).
  double speedGain{0.000004};

  // Gain converting velocity error into acceleration demand (1/s):
  //   aWanted = dv * velGain.
  double velGain{1.8};

  // Desired stand-off distance from the target (km). Use 0 to "arrive at" the point.
  double desiredDistKm{0.0};

  // Back-off behavior when we are well inside desiredDist (helps avoid sitting on top).
  double backoffFrac{0.60};     // inside desiredDist * frac, request a small negative closing speed
  double backoffSpeedKmS{0.10}; // magnitude when backing off

  // Scales ship acceleration caps (useful for different AI tiers).
  double accelScale{1.0};

  // If true, enable boost when the demanded acceleration exceeds the base cap.
  bool allowBoost{false};

  // Whether to enable dampers in the resulting ShipInput.
  bool dampers{true};
};

struct AttitudeControlParams {
  // Gains mapping angular error (radians) to normalized [-1,1] torque command.
  double faceGain{1.6};
  double rollGain{1.6};

  // If true, compute roll to align up vector (requires desiredUpWorld).
  bool alignUp{false};
};

// Optional intercept-course ("lead") tuning for chasing moving targets.
//
// When a target has significant lateral velocity, a naive pursuit controller can end up
// "tail chasing" and orbiting around it. Enabling intercept-course guidance computes a
// lead direction (using the same closed-form intercept solve as projectile lead) and
// uses that direction for the translation velocity setpoint.
//
// This is still an approximation (ships are acceleration-limited, not speed-limited),
// but it significantly improves chase behavior in practice without requiring a full
// guidance law.
struct InterceptCourseParams {
  // Master toggle.
  bool enabled{true};

  // Reject lead solutions that would require aiming too far into the future.
  // Large values can look "psychic" for very distant targets.
  double maxLeadTimeSec{120.0};

  // Only compute a lead solution when the commanded closing speed is at least this.
  // Prevents noisy/unstable lead directions when maneuvering slowly near a target.
  double minSpeedKmS{0.05};

  // If true, use params.maxSpeedKmS for the lead solve instead of the current desired
  // closing speed. This tends to yield more stable lead directions at the cost of a
  // small loss of near-range accuracy.
  bool useMaxSpeedForSolve{true};
};

struct FlightControlOutput {
  ShipInput input{};
  math::Vec3d desiredVelKmS{0,0,0};
  bool usedBoost{false};
  double distKm{0.0};
};

// Approach a moving target position/velocity while trying to maintain a standoff distance.
// Provide a desiredForwardWorld (what the ship's forward should point at), and optionally
// desiredUpWorld for roll alignment.
FlightControlOutput approachTarget(
  const Ship& ship,
  const math::Vec3d& targetPosKm,
  const math::Vec3d& targetVelKmS,
  const FlightControlParams& params,
  const AttitudeControlParams& attitude,
  const math::Vec3d& desiredForwardWorld,
  const math::Vec3d* desiredUpWorld = nullptr);

// Like approachTarget(), but uses an intercept-course lead direction for translation.
//
// The caller still chooses the facing direction (desiredForwardWorld). This allows
// AI to keep moving on an intercept course while aiming elsewhere (e.g. projectile lead).
FlightControlOutput approachTargetIntercept(
  const Ship& ship,
  const math::Vec3d& targetPosKm,
  const math::Vec3d& targetVelKmS,
  const FlightControlParams& params,
  const AttitudeControlParams& attitude,
  const math::Vec3d& desiredForwardWorld,
  const InterceptCourseParams& intercept,
  const math::Vec3d* desiredUpWorld = nullptr);

// Convenience wrapper: desired forward points toward the target.
FlightControlOutput chaseTarget(
  const Ship& ship,
  const math::Vec3d& targetPosKm,
  const math::Vec3d& targetVelKmS,
  const FlightControlParams& params,
  const AttitudeControlParams& attitude);

// Convenience wrapper: faces the target like chaseTarget(), but uses intercept-course
// translation guidance.
FlightControlOutput chaseTargetIntercept(
  const Ship& ship,
  const math::Vec3d& targetPosKm,
  const math::Vec3d& targetVelKmS,
  const FlightControlParams& params,
  const AttitudeControlParams& attitude,
  const InterceptCourseParams& intercept);

} // namespace stellar::sim
