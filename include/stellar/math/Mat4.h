#pragma once
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/math/Math.h"

#include <cmath>

namespace stellar::math {

// Column-major 4x4 matrix (OpenGL-friendly).
struct Mat4d {
  double m[16]{};

  static Mat4d identity() {
    Mat4d r{};
    r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1.0;
    return r;
  }

  static Mat4d translation(const Vec3d& t) {
    Mat4d r = identity();
    r.m[12] = t.x;
    r.m[13] = t.y;
    r.m[14] = t.z;
    return r;
  }

  static Mat4d scale(double s) {
    Mat4d r{};
    r.m[0] = r.m[5] = r.m[10] = s;
    r.m[15] = 1.0;
    return r;
  }

  static Mat4d scale(const Vec3d& s) {
    Mat4d r{};
    r.m[0] = s.x;
    r.m[5] = s.y;
    r.m[10] = s.z;
    r.m[15] = 1.0;
    return r;
  }

  static Mat4d rotation(const Quatd& qIn) {
    const Quatd q = qIn.normalized();
    const double xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z;
    const double xy = q.x*q.y, xz = q.x*q.z, yz = q.y*q.z;
    const double wx = q.w*q.x, wy = q.w*q.y, wz = q.w*q.z;

    Mat4d r = identity();
    r.m[0]  = 1.0 - 2.0*(yy + zz);
    r.m[1]  = 2.0*(xy + wz);
    r.m[2]  = 2.0*(xz - wy);

    r.m[4]  = 2.0*(xy - wz);
    r.m[5]  = 1.0 - 2.0*(xx + zz);
    r.m[6]  = 2.0*(yz + wx);

    r.m[8]  = 2.0*(xz + wy);
    r.m[9]  = 2.0*(yz - wx);
    r.m[10] = 1.0 - 2.0*(xx + yy);
    return r;
  }

  static Mat4d perspective(double fovYRad, double aspect, double zNear, double zFar) {
    Mat4d r{};
    const double f = 1.0 / std::tan(fovYRad * 0.5);
    r.m[0] = f / aspect;
    r.m[5] = f;
    r.m[10] = (zFar + zNear) / (zNear - zFar);
    r.m[11] = -1.0;
    r.m[14] = (2.0 * zFar * zNear) / (zNear - zFar);
    return r;
  }

  static Mat4d lookAt(const Vec3d& eye, const Vec3d& center, const Vec3d& up) {
    const Vec3d f = (center - eye).normalized();
    const Vec3d s = cross(f, up.normalized()).normalized();
    const Vec3d u = cross(s, f);

    Mat4d r = identity();
    r.m[0] = s.x; r.m[4] = s.y; r.m[8]  = s.z;
    r.m[1] = u.x; r.m[5] = u.y; r.m[9]  = u.z;
    r.m[2] = -f.x; r.m[6] = -f.y; r.m[10] = -f.z;

    r.m[12] = -dot(s, eye);
    r.m[13] = -dot(u, eye);
    r.m[14] = dot(f, eye);
    return r;
  }

  Mat4d operator*(const Mat4d& b) const {
    Mat4d r{};
    for (int c = 0; c < 4; ++c) {
      for (int rrow = 0; rrow < 4; ++rrow) {
        r.m[c*4 + rrow] =
          m[0*4 + rrow] * b.m[c*4 + 0] +
          m[1*4 + rrow] * b.m[c*4 + 1] +
          m[2*4 + rrow] * b.m[c*4 + 2] +
          m[3*4 + rrow] * b.m[c*4 + 3];
      }
    }
    return r;
  }
};

} // namespace stellar::math
