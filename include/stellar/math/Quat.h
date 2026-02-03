#pragma once
#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>

namespace stellar::math {

struct Quatd {
  double w{1}, x{0}, y{0}, z{0};

  constexpr Quatd() = default;
  constexpr Quatd(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}

  static Quatd identity() { return {1, 0, 0, 0}; }

  double normSq() const { return w * w + x * x + y * y + z * z; }

  // Quaternion dot product.
  static double dot(const Quatd& a, const Quatd& b) {
    return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
  }

  // Negation represents the same rotation (q and -q encode identical orientation).
  constexpr Quatd operator-() const { return {-w, -x, -y, -z}; }

  static Quatd fromAxisAngle(const Vec3d& axis, double radians) {
    // Defensive: if the axis is degenerate or inputs are non-finite, return identity.
    const double lsq = axis.lengthSq();
    if (!(lsq > 1e-24) || !std::isfinite(radians)) return identity();

    const Vec3d n = axis * (1.0 / std::sqrt(lsq));
    const double half = radians * 0.5;
    const double s = std::sin(half);
    return Quatd{std::cos(half), n.x * s, n.y * s, n.z * s}.normalized();
  }

  // Construct a quaternion that rotates vector `from` onto vector `to`.
  // If either vector is near zero, returns identity.
  // Construct a quaternion that rotates vector `from` onto vector `to`.
  //
  // Notes:
  // - If either vector is near zero, returns identity.
  // - Handles the 180-degree case (opposite vectors) by choosing an arbitrary perpendicular axis.
  // - Uses a robust "1+dot, cross" formulation so very small angles still produce the correct
  //   (near-identity) rotation rather than snapping to identity.
  static Quatd fromTo(const Vec3d& from, const Vec3d& to) {
    const double fromLenSq = from.lengthSq();
    const double toLenSq = to.lengthSq();
    if (!(fromLenSq > 1e-18) || !(toLenSq > 1e-18)) return identity();

    const Vec3d f = from * (1.0 / std::sqrt(fromLenSq));
    const Vec3d t = to * (1.0 / std::sqrt(toLenSq));

    const double c = std::clamp(math::dot(f, t), -1.0, 1.0);

    if (c < -0.999999) {
      // 180 degree rotation around any axis orthogonal to `from`.
      Vec3d axis = cross({1, 0, 0}, f);
      if (axis.lengthSq() < 1e-12) axis = cross({0, 1, 0}, f);
      return fromAxisAngle(axis, math::kPi);
    }

    const Vec3d axis = cross(f, t);

    // For all cases except the 180-degree flip, this is stable (including tiny angles).
    // q = [1 + dot(f,t), cross(f,t)] then normalized.
    return Quatd{1.0 + c, axis.x, axis.y, axis.z}.normalized();
  }

  // Construct from an orthonormal basis whose columns are (right, up, forward).
  static Quatd fromBasis(const Vec3d& right, const Vec3d& up, const Vec3d& forward) {
    const double m00 = right.x;
    const double m01 = up.x;
    const double m02 = forward.x;

    const double m10 = right.y;
    const double m11 = up.y;
    const double m12 = forward.y;

    const double m20 = right.z;
    const double m21 = up.z;
    const double m22 = forward.z;

    const double trace = m00 + m11 + m22;
    double qw = 1.0, qx = 0.0, qy = 0.0, qz = 0.0;

    if (trace > 0.0) {
      const double s = std::sqrt(trace + 1.0) * 2.0;
      qw = 0.25 * s;
      qx = (m21 - m12) / s;
      qy = (m02 - m20) / s;
      qz = (m10 - m01) / s;
    } else if (m00 > m11 && m00 > m22) {
      const double s = std::sqrt(1.0 + m00 - m11 - m22) * 2.0;
      qw = (m21 - m12) / s;
      qx = 0.25 * s;
      qy = (m01 + m10) / s;
      qz = (m02 + m20) / s;
    } else if (m11 > m22) {
      const double s = std::sqrt(1.0 + m11 - m00 - m22) * 2.0;
      qw = (m02 - m20) / s;
      qx = (m01 + m10) / s;
      qy = 0.25 * s;
      qz = (m12 + m21) / s;
    } else {
      const double s = std::sqrt(1.0 + m22 - m00 - m11) * 2.0;
      qw = (m10 - m01) / s;
      qx = (m02 + m20) / s;
      qy = (m12 + m21) / s;
      qz = 0.25 * s;
    }

    return Quatd{qw, qx, qy, qz}.normalized();
  }

  static Quatd lookRotation(const Vec3d& forwardDir, const Vec3d& upHint) {
    const Vec3d f = safeNormalized(forwardDir, {0.0, 0.0, 1.0});
    Vec3d u = safeNormalized(upHint, {0.0, 1.0, 0.0});

    // If up is nearly parallel to forward, pick a fallback up.
    if (std::abs(math::dot(u, f)) > 0.999) {
      u = (std::abs(f.y) < 0.999) ? Vec3d{0.0, 1.0, 0.0} : Vec3d{1.0, 0.0, 0.0};
    }

    Vec3d r = cross(u, f);
    r = safeNormalized(r, {1.0, 0.0, 0.0});
    const Vec3d u2 = cross(f, r);

    return fromBasis(r, u2, f);
  }

  static Quatd fromEulerYXZ(double yawRad, double pitchRad, double rollRad) {
    // yaw (Y), pitch (X), roll (Z) order
    const double cy = std::cos(yawRad * 0.5), sy = std::sin(yawRad * 0.5);
    const double cx = std::cos(pitchRad * 0.5), sx = std::sin(pitchRad * 0.5);
    const double cz = std::cos(rollRad * 0.5), sz = std::sin(rollRad * 0.5);

    Quatd qy{cy, 0, sy, 0};
    Quatd qx{cx, sx, 0, 0};
    Quatd qz{cz, 0, 0, sz};
    return (qy * qx * qz).normalized();
  }

  // Normalized linear interpolation with shortest-arc handling.
  static Quatd nlerp(const Quatd& a, const Quatd& bIn, double t) {
    t = std::clamp(t, 0.0, 1.0);

    Quatd b = bIn;
    if (dot(a, b) < 0.0) b = -b;

    Quatd q{
      a.w * (1.0 - t) + b.w * t,
      a.x * (1.0 - t) + b.x * t,
      a.y * (1.0 - t) + b.y * t,
      a.z * (1.0 - t) + b.z * t,
    };

    return q.normalized();
  }

  // Spherical linear interpolation with shortest-arc handling.
  // Falls back to nlerp when quaternions are nearly identical.
  static Quatd slerp(const Quatd& a, const Quatd& bIn, double t) {
    t = std::clamp(t, 0.0, 1.0);

    Quatd b = bIn;
    double d = dot(a, b);

    // Ensure shortest arc.
    if (d < 0.0) {
      b = -b;
      d = -d;
    }

    // If the angle is tiny, nlerp is stable and effectively identical.
    if (d > 0.9995) {
      return nlerp(a, b, t);
    }

    d = std::clamp(d, -1.0, 1.0);
    const double theta = std::acos(d);
    const double sinTheta = std::sin(theta);
    if (!(sinTheta > 1e-12)) {
      return nlerp(a, b, t);
    }

    const double wa = std::sin((1.0 - t) * theta) / sinTheta;
    const double wb = std::sin(t * theta) / sinTheta;

    Quatd q{
      a.w * wa + b.w * wb,
      a.x * wa + b.x * wb,
      a.y * wa + b.y * wb,
      a.z * wa + b.z * wb,
    };

    return q.normalized();
  }

  Quatd normalized(double eps = 1e-12) const {
    const double n = std::sqrt(normSq());
    if (!(n > eps)) return identity();
    return {w / n, x / n, y / n, z / n};
  }

  Quatd conjugate() const { return {w, -x, -y, -z}; }

  Quatd inverse(double eps = 1e-12) const {
    const double n2 = normSq();
    if (!(n2 > eps * eps)) return identity();
    return {w / n2, -x / n2, -y / n2, -z / n2};
  }

  Quatd operator*(const Quatd& b) const {
    return {
      w * b.w - x * b.x - y * b.y - z * b.z,
      w * b.x + x * b.w + y * b.z - z * b.y,
      w * b.y - x * b.z + y * b.w + z * b.x,
      w * b.z + x * b.y - y * b.x + z * b.w,
    };
  }

  Vec3d rotate(const Vec3d& v) const {
    // Rotate by q * v * q^-1.
    //
    // Implemented via a fast vector identity that is scale-invariant (works even
    // if this quaternion is not perfectly normalized):
    //   t  = 2 * cross(q.xyz, v) / |q|^2
    //   v' = v + w*t + cross(q.xyz, t)
    const double n2 = normSq();
    if (!(n2 > 1e-24)) return v;

    const Vec3d u{x, y, z};
    const Vec3d t = (2.0 / n2) * cross(u, v);
    return v + t * w + cross(u, t);
  }

  // Convenience: rotate a vector by this quaternion.
  // Enables syntax like: Vec3d world = q * local;
  Vec3d operator*(const Vec3d& v) const { return rotate(v); }

  // Integrate by angular velocity (rad/day or rad/sec depending on caller units)
  Quatd integrateAngular(const Vec3d& omega, double dt) const {
    // Exact integration for constant body-frame angular velocity:
    //   q(t+dt) = q(t) * exp(0.5 * omega * dt)
    // where omega is a pure-imag quaternion (0, wx, wy, wz).

    const double wx = omega.x;
    const double wy = omega.y;
    const double wz = omega.z;

    if (!std::isfinite(dt) || !std::isfinite(wx) || !std::isfinite(wy) || !std::isfinite(wz)) {
      return normalized();
    }

    const double w2 = wx * wx + wy * wy + wz * wz;
    if (!(w2 > 1e-24) || dt == 0.0) {
      // Preserve normalization even if callers feed slightly drifted quats.
      return normalized();
    }

    const double wmag = std::sqrt(w2);
    const double halfTheta = 0.5 * wmag * dt;

    // For the imaginary part we need: axis * sin(halfTheta)
    // where axis = omega / |omega|. Rewrite as omega * (sin(halfTheta)/|omega|).
    // This avoids a separate normalization and is stable for tiny angles.
    const double c = std::cos(halfTheta);
    double scale = 0.0;
    if (std::abs(halfTheta) < 1e-6) {
      // sin(x)/x ~ 1 - x^2/6 for small x.
      const double x2 = halfTheta * halfTheta;
      const double sinc = 1.0 - x2 / 6.0;
      scale = 0.5 * dt * sinc; // sin(halfTheta)/|omega| = (sin(x)/x) * (0.5*dt)
    } else {
      scale = std::sin(halfTheta) / wmag;
    }

    const Quatd dq{c, wx * scale, wy * scale, wz * scale};
    return ((*this) * dq).normalized();
  }
};

} // namespace stellar::math
