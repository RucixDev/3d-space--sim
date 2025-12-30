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

// Convenience wrapper: desired forward points toward the target.
FlightControlOutput chaseTarget(
  const Ship& ship,
  const math::Vec3d& targetPosKm,
  const math::Vec3d& targetVelKmS,
  const FlightControlParams& params,
  const AttitudeControlParams& attitude);

} // namespace stellar::sim
