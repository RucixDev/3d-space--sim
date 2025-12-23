#include "stellar/sim/Ship.h"

#include "stellar/math/Math.h"

#include <algorithm>

namespace stellar::sim {
namespace {
constexpr double kSecondsPerDay = 86400.0;
}

stellar::math::Vec3d Ship::forward() const {
  const double yaw = stellar::math::degToRad(yawDeg);
  const double pitch = stellar::math::degToRad(pitchDeg);

  const double cy = std::cos(yaw);
  const double sy = std::sin(yaw);
  const double cp = std::cos(pitch);
  const double sp = std::sin(pitch);

  // Coordinate system:
  // - X right
  // - Y up
  // - Z forward
  return {sy * cp, sp, cy * cp};
}

stellar::math::Vec3d Ship::right() const {
  const double yaw = stellar::math::degToRad(yawDeg);
  const double cy = std::cos(yaw);
  const double sy = std::sin(yaw);
  return {cy, 0.0, -sy};
}

stellar::math::Vec3d Ship::up() const {
  // Ensure an orthonormal basis.
  const auto f = forward();
  const auto r = right();
  return stellar::math::cross(r, f).normalized();
}

void Ship::step(const ShipControls& c, double dtSeconds) {
  // Apply look
  yawDeg += c.yawDeltaDeg;
  pitchDeg += c.pitchDeltaDeg;

  // Keep pitch sane.
  pitchDeg = std::clamp(pitchDeg, -89.9, 89.9);

  const double dtDays = dtSeconds / kSecondsPerDay;

  // Thrust (local ship frame -> world)
  const double accel = c.boost ? boostAccelAUPerDay2 : maxAccelAUPerDay2;

  const auto f = forward();
  const auto r = right();
  const auto u = up();

  const stellar::math::Vec3d aWorldAUPerDay2 =
    (f * c.thrustForward + r * c.thrustRight + u * c.thrustUp) * accel;

  velocityAUPerDay += aWorldAUPerDay2 * dtDays;

  if (c.brake) {
    const double k = std::clamp(brakeStrengthPerSecond * dtSeconds, 0.0, 1.0);
    velocityAUPerDay *= (1.0 - k);
  }

  // Clamp speed
  const double speed = velocityAUPerDay.length();
  if (speed > maxSpeedAUPerDay && speed > 0.0) {
    velocityAUPerDay *= (maxSpeedAUPerDay / speed);
  }

  // Integrate position
  positionAU += velocityAUPerDay * dtDays;
}

} // namespace stellar::sim
