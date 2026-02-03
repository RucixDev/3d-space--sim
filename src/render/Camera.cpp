#include "stellar/render/Camera.h"

namespace stellar::render {

void Camera::setPerspective(double fovYRad, double aspect, double zNear, double zFar) {
  proj_ = stellar::math::Mat4d::perspective(fovYRad, aspect, zNear, zFar);
}

stellar::math::Mat4d Camera::viewMatrix() const {
  // View = (T * R)^-1 where world = T(pos) * R(orient).
  // Using inverseRigid keeps the intent clear and centralizes the math.
  const auto world = stellar::math::Mat4d::rigidTransform(pos_, orient_);
  return world.inverseRigid();
}

} // namespace stellar::render
