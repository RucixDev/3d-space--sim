#pragma once

#include "stellar/math/Vec3.h"

namespace stellar::sim {

// Input state for a simple 6DOF-ish ship model.
// Thrust values are expected to be in [-1, 1].
struct ShipControls {
  double thrustForward = 0.0; // + forward, - backward
  double thrustRight = 0.0;   // + right, - left
  double thrustUp = 0.0;      // + up, - down

  // Mouse/keyboard look deltas in degrees (applied this frame)
  double yawDeltaDeg = 0.0;
  double pitchDeltaDeg = 0.0;

  bool brake = false;
  bool boost = false;
};

// Lightweight player ship state.
// Units:
// - positionAU: AU
// - velocityAUPerDay: AU/day
struct Ship {
  stellar::math::Vec3d positionAU{0.0, 0.0, 5.0};
  stellar::math::Vec3d velocityAUPerDay{0.0, 0.0, 0.0};

  double yawDeg = 0.0;
  double pitchDeg = 0.0;

  // Tuning
  double maxAccelAUPerDay2 = 25.0;
  double boostAccelAUPerDay2 = 75.0;
  double maxSpeedAUPerDay = 8.0;
  double brakeStrengthPerSecond = 3.0;

  // Orientation basis vectors in world space
  stellar::math::Vec3d forward() const;
  stellar::math::Vec3d right() const;
  stellar::math::Vec3d up() const;

  // Integrate ship state.
  void step(const ShipControls& c, double dtSeconds);
};

} // namespace stellar::sim
