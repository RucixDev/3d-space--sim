#pragma once

#include "stellar/math/Mat4.h"
#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

namespace stellar::render {

class Camera {
public:
  stellar::math::Vec3d position{0.0, 0.0, -10.0};
  double yawDeg = 0.0;
  double pitchDeg = 0.0;

  double fovDeg = 70.0;
  double nearPlane = 0.01;
  double farPlane = 5000.0;

  stellar::math::Vec3d forward() const {
    const double yaw = stellar::math::degToRad(yawDeg);
    const double pitch = stellar::math::degToRad(pitchDeg);
    const double cy = std::cos(yaw);
    const double sy = std::sin(yaw);
    const double cp = std::cos(pitch);
    const double sp = std::sin(pitch);
    return {sy * cp, sp, cy * cp};
  }

  stellar::math::Vec3d right() const {
    const double yaw = stellar::math::degToRad(yawDeg);
    const double cy = std::cos(yaw);
    const double sy = std::sin(yaw);
    return {cy, 0.0, -sy};
  }

  stellar::math::Vec3d up() const {
    return stellar::math::cross(right(), forward()).normalized();
  }

  stellar::math::Mat4f view() const {
    const auto f = forward();
    return stellar::math::Mat4f::lookAt(position, position + f, {0.0, 1.0, 0.0});
  }

  stellar::math::Mat4f projection(float aspect) const {
    return stellar::math::Mat4f::perspective(
      static_cast<float>(stellar::math::degToRad(fovDeg)),
      aspect,
      static_cast<float>(nearPlane),
      static_cast<float>(farPlane));
  }

  stellar::math::Mat4f viewProj(float aspect) const {
    return projection(aspect) * view();
  }
};

} // namespace stellar::render
