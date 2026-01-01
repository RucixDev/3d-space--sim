#pragma once
#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

#include <cmath>

namespace stellar::math {

struct Quatd {
  double w{1}, x{0}, y{0}, z{0};

  constexpr Quatd() = default;
  constexpr Quatd(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}

  static Quatd identity() { return {1,0,0,0}; }

  static Quatd fromAxisAngle(const Vec3d& axis, double radians) {
    const Vec3d n = axis.normalized();
    const double half = radians * 0.5;
    const double s = std::sin(half);
    return { std::cos(half), n.x*s, n.y*s, n.z*s };
  }

  static Quatd fromEulerYXZ(double yawRad, double pitchRad, double rollRad) {
    // yaw (Y), pitch (X), roll (Z) order
    const double cy = std::cos(yawRad*0.5),   sy = std::sin(yawRad*0.5);
    const double cx = std::cos(pitchRad*0.5), sx = std::sin(pitchRad*0.5);
    const double cz = std::cos(rollRad*0.5),  sz = std::sin(rollRad*0.5);

    Quatd qy{cy, 0, sy, 0};
    Quatd qx{cx, sx, 0, 0};
    Quatd qz{cz, 0, 0, sz};
    return (qy * qx * qz).normalized();
  }

  Quatd normalized(double eps=1e-12) const {
    const double n = std::sqrt(w*w + x*x + y*y + z*z);
    if (n < eps) return identity();
    return {w/n, x/n, y/n, z/n};
  }

  Quatd conjugate() const { return {w, -x, -y, -z}; }

  Quatd operator*(const Quatd& b) const {
    return {
      w*b.w - x*b.x - y*b.y - z*b.z,
      w*b.x + x*b.w + y*b.z - z*b.y,
      w*b.y - x*b.z + y*b.w + z*b.x,
      w*b.z + x*b.y - y*b.x + z*b.w
    };
  }

  Vec3d rotate(const Vec3d& v) const {
    // q * (0,v) * q^-1
    Quatd qv{0, v.x, v.y, v.z};
    Quatd r = (*this) * qv * conjugate();
    return {r.x, r.y, r.z};
  }

  // Convenience: rotate a vector by this quaternion.
  // Enables syntax like: Vec3d world = q * local;
  Vec3d operator*(const Vec3d& v) const { return rotate(v); }

  // Integrate by angular velocity (rad/day or rad/sec depending on caller units)
  Quatd integrateAngular(const Vec3d& omega, double dt) const {
    // dq/dt = 0.5 * q * (0, omega)
    Quatd o{0, omega.x, omega.y, omega.z};
    Quatd dq = (*this) * o;
    Quatd q2{ w + 0.5 * dq.w * dt,
              x + 0.5 * dq.x * dt,
              y + 0.5 * dq.y * dt,
              z + 0.5 * dq.z * dt };
    return q2.normalized();
  }
};

} // namespace stellar::math
