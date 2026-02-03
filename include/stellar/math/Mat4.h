#pragma once
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/math/Vec4.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::math {

// Column-major 4x4 matrix (OpenGL-friendly).
struct Mat4d {
  double m[16]{};

  bool isFinite() const {
    for (double v : m) {
      if (!std::isfinite(v)) return false;
    }
    return true;
  }

  static Mat4d identity() {
    Mat4d r{};
    r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1.0;
    return r;
  }

  // Copy this matrix into a column-major float[16] array (OpenGL-friendly).
  void toFloat(float out[16]) const {
    for (int i = 0; i < 16; ++i) out[i] = static_cast<float>(m[i]);
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

  // Construct a rigid transform matrix representing world = T(t) * R(q).
  // This is a convenience wrapper that avoids repeated callsites doing
  // translation(t) * rotation(q).
  static Mat4d rigidTransform(const Vec3d& t, const Quatd& qIn) {
    const Quatd q = qIn.normalized();
    const double xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z;
    const double xy = q.x*q.y, xz = q.x*q.z, yz = q.y*q.z;
    const double wx = q.w*q.x, wy = q.w*q.y, wz = q.w*q.z;

    Mat4d r{};
    r.m[0]  = 1.0 - 2.0*(yy + zz);
    r.m[1]  = 2.0*(xy + wz);
    r.m[2]  = 2.0*(xz - wy);
    r.m[3]  = 0.0;

    r.m[4]  = 2.0*(xy - wz);
    r.m[5]  = 1.0 - 2.0*(xx + zz);
    r.m[6]  = 2.0*(yz + wx);
    r.m[7]  = 0.0;

    r.m[8]  = 2.0*(xz + wy);
    r.m[9]  = 2.0*(yz - wx);
    r.m[10] = 1.0 - 2.0*(xx + yy);
    r.m[11] = 0.0;

    r.m[12] = t.x;
    r.m[13] = t.y;
    r.m[14] = t.z;
    r.m[15] = 1.0;
    return r;
  }

  // Construct an affine TRS matrix representing world = T(t) * R(q) * S(s).
  // Scale is applied in local space (i.e., scales the columns of R).
  static Mat4d trs(const Vec3d& t, const Quatd& qIn, const Vec3d& s) {
    const Quatd q = qIn.normalized();
    const double xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z;
    const double xy = q.x*q.y, xz = q.x*q.z, yz = q.y*q.z;
    const double wx = q.w*q.x, wy = q.w*q.y, wz = q.w*q.z;

    // Rotation matrix (column-major) scaled by (s.x, s.y, s.z) per column.
    const double r00 = 1.0 - 2.0*(yy + zz);
    const double r01 = 2.0*(xy + wz);
    const double r02 = 2.0*(xz - wy);

    const double r10 = 2.0*(xy - wz);
    const double r11 = 1.0 - 2.0*(xx + zz);
    const double r12 = 2.0*(yz + wx);

    const double r20 = 2.0*(xz + wy);
    const double r21 = 2.0*(yz - wx);
    const double r22 = 1.0 - 2.0*(xx + yy);

    Mat4d r{};
    // Column 0 (X axis)
    r.m[0] = r00 * s.x;
    r.m[1] = r01 * s.x;
    r.m[2] = r02 * s.x;
    r.m[3] = 0.0;

    // Column 1 (Y axis)
    r.m[4] = r10 * s.y;
    r.m[5] = r11 * s.y;
    r.m[6] = r12 * s.y;
    r.m[7] = 0.0;

    // Column 2 (Z axis)
    r.m[8]  = r20 * s.z;
    r.m[9]  = r21 * s.z;
    r.m[10] = r22 * s.z;
    r.m[11] = 0.0;

    r.m[12] = t.x;
    r.m[13] = t.y;
    r.m[14] = t.z;
    r.m[15] = 1.0;
    return r;
  }

  static Mat4d trs(const Vec3d& t, const Quatd& qIn, double s) {
    return trs(t, qIn, {s, s, s});
  }

  Vec3d translationPart() const {
    return {m[12], m[13], m[14]};
  }

  // Extract the rotation as a quaternion. If the matrix contains scale, the
  // columns are normalized before converting.
  Quatd rotationPart() const {
    const Vec3d right = safeNormalized({m[0], m[1], m[2]}, {1.0, 0.0, 0.0}, 1e-18);
    const Vec3d up = safeNormalized({m[4], m[5], m[6]}, {0.0, 1.0, 0.0}, 1e-18);
    const Vec3d forward = safeNormalized({m[8], m[9], m[10]}, {0.0, 0.0, 1.0}, 1e-18);
    return Quatd::fromBasis(right, up, forward);
  }

  static Mat4d perspective(double fovYRad, double aspect, double zNear, double zFar) {
    // Defensive: avoid NaNs when callers feed bad camera parameters.
    if (!std::isfinite(fovYRad) || !std::isfinite(aspect) ||
        !std::isfinite(zNear) || !std::isfinite(zFar)) {
      return identity();
    }

    // Clamp aspect to a sane default.
    const double a = (aspect > 1e-12) ? aspect : 1.0;

    // Avoid tan(fov/2) blowing up at 0 or pi.
    const double fov = std::clamp(fovYRad, 1e-6, kPi - 1e-6);

    // Ensure a valid near/far pair.
    const double zn = std::max(1e-6, zNear);
    double zf = zFar;
    if (!(zf > zn + 1e-6)) {
      // Preserve ordering rather than emitting an invalid projection.
      zf = zn + 1.0;
    }

    Mat4d r{};
    const double f = 1.0 / std::tan(fov * 0.5);
    r.m[0] = f / a;
    r.m[5] = f;
    r.m[10] = (zf + zn) / (zn - zf);
    r.m[11] = -1.0;
    r.m[14] = (2.0 * zf * zn) / (zn - zf);
    return r;
  }

  static Mat4d lookAt(const Vec3d& eye, const Vec3d& center, const Vec3d& up) {
    // Robust lookAt that remains well-defined even when the caller provides:
    //  - eye == center (zero forward vector)
    //  - an up vector parallel to forward
    //
    // We choose stable fallbacks rather than returning an all-zero matrix.

    Vec3d f = safeNormalized(center - eye, {0.0, 0.0, -1.0}, 1e-18);
    Vec3d uHint = safeNormalized(up, {0.0, 1.0, 0.0}, 1e-18);

    // If up is nearly parallel to forward, pick a different up axis.
    if (std::abs(dot(uHint, f)) > 0.999) {
      uHint = (std::abs(f.y) < 0.999) ? Vec3d{0.0, 1.0, 0.0} : Vec3d{1.0, 0.0, 0.0};
    }

    auto fallbackRight = [&](const Vec3d& forward) {
      Vec3d a = (std::abs(forward.x) < 0.9) ? Vec3d{1.0, 0.0, 0.0} : Vec3d{0.0, 1.0, 0.0};
      Vec3d r = cross(forward, a);
      if (r.lengthSq() < 1e-18) {
        r = cross(forward, Vec3d{0.0, 0.0, 1.0});
      }
      return safeNormalized(r, {1.0, 0.0, 0.0}, 1e-18);
    };

    Vec3d s = cross(f, uHint);
    if (s.lengthSq() < 1e-18) {
      s = fallbackRight(f);
    } else {
      s = s.normalized();
    }

    Vec3d u2 = cross(s, f);
    u2 = safeNormalized(u2, {0.0, 1.0, 0.0}, 1e-18);

    // Re-orthonormalize right in case the up hint was slightly skew.
    s = safeNormalized(cross(f, u2), s, 1e-18);

    Mat4d r = identity();
    r.m[0] = s.x;  r.m[4] = s.y;  r.m[8]  = s.z;
    r.m[1] = u2.x; r.m[5] = u2.y; r.m[9]  = u2.z;
    r.m[2] = -f.x; r.m[6] = -f.y; r.m[10] = -f.z;

    r.m[12] = -dot(s, eye);
    r.m[13] = -dot(u2, eye);
    r.m[14] = dot(f, eye);
    return r;
  }

  Vec4d mul(const Vec4d& v) const {
    // Column-major matrix-vector multiply.
    return {
      m[0] * v.x + m[4] * v.y + m[8]  * v.z + m[12] * v.w,
      m[1] * v.x + m[5] * v.y + m[9]  * v.z + m[13] * v.w,
      m[2] * v.x + m[6] * v.y + m[10] * v.z + m[14] * v.w,
      m[3] * v.x + m[7] * v.y + m[11] * v.z + m[15] * v.w,
    };
  }

  Vec4d mulPoint(const Vec3d& p, double w = 1.0) const {
    return mul({p.x, p.y, p.z, w});
  }

  // Transform a point assuming an affine matrix (w stays 1).
  Vec3d transformPointAffine(const Vec3d& p) const {
    const auto v = mulPoint(p, 1.0);
    return {v.x, v.y, v.z};
  }

  // Transform a direction vector (ignores translation; w=0).
  Vec3d transformVector(const Vec3d& v) const {
    const auto r = mulPoint(v, 0.0);
    return {r.x, r.y, r.z};
  }

  Mat4d transposed() const {
    Mat4d r{};
    for (int c = 0; c < 4; ++c) {
      for (int row = 0; row < 4; ++row) {
        r.m[c * 4 + row] = m[row * 4 + c];
      }
    }
    return r;
  }

  // Invert a rigid transform (rotation + translation, no scale/shear).
  //
  // Assumes the matrix is of the form:
  //   [ R t ]
  //   [ 0 1 ]
  // where R is orthonormal. The inverse is:
  //   [ R^T -R^T t ]
  //   [ 0     1    ]
  Mat4d inverseRigid() const {
    if (!isFinite()) return identity();

    Mat4d r = identity();

    // Transpose the upper-left 3x3.
    r.m[0] = m[0];
    r.m[1] = m[4];
    r.m[2] = m[8];

    r.m[4] = m[1];
    r.m[5] = m[5];
    r.m[6] = m[9];

    r.m[8]  = m[2];
    r.m[9]  = m[6];
    r.m[10] = m[10];

    const Vec3d t{m[12], m[13], m[14]};

    // -R^T * t
    r.m[12] = -(r.m[0] * t.x + r.m[4] * t.y + r.m[8]  * t.z);
    r.m[13] = -(r.m[1] * t.x + r.m[5] * t.y + r.m[9]  * t.z);
    r.m[14] = -(r.m[2] * t.x + r.m[6] * t.y + r.m[10] * t.z);

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
