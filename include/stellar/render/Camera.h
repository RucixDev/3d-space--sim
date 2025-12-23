#pragma once

#include "stellar/math/Mat4.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"

namespace stellar::render {

class Camera {
public:
  void setPerspective(double fovYRad, double aspect, double zNear, double zFar);

  void setPosition(const stellar::math::Vec3d& p) { pos_ = p; }
  void setOrientation(const stellar::math::Quatd& q) { orient_ = q; }

  stellar::math::Vec3d position() const { return pos_; }
  stellar::math::Quatd orientation() const { return orient_; }

  stellar::math::Mat4d viewMatrix() const;
  stellar::math::Mat4d projectionMatrix() const { return proj_; }

private:
  stellar::math::Vec3d pos_{0,0,0};
  stellar::math::Quatd orient_{1,0,0,0};
  stellar::math::Mat4d proj_{stellar::math::Mat4d::identity()};
};

} // namespace stellar::render
