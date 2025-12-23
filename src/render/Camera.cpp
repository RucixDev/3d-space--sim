#include "stellar/render/Camera.h"

namespace stellar::render {

void Camera::setPerspective(double fovYRad, double aspect, double zNear, double zFar) {
  proj_ = stellar::math::Mat4d::perspective(fovYRad, aspect, zNear, zFar);
}

stellar::math::Mat4d Camera::viewMatrix() const {
  // View = R^-1 * T^-1
  const stellar::math::Quatd inv = orient_.conjugate();
  const auto rot = stellar::math::Mat4d::rotation(inv);
  const auto tr = stellar::math::Mat4d::translation({-pos_.x, -pos_.y, -pos_.z});
  return rot * tr;
}

} // namespace stellar::render
