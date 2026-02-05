#pragma once

#include "stellar/math/Geometry.h"
#include "stellar/math/Mat4.h"

#include <array>
#include <algorithm>
#include <cmath>

namespace stellar::math {

// Plane represented as: n·x + d = 0.
//
// By convention for this header's helpers, points are considered "inside" a plane's
// half-space when:
//   signedDistance(p) >= 0
struct Plane3d {
  Vec3d n{0.0, 0.0, 1.0};
  double d{0.0};

  bool isFinite() const {
    return stellar::math::isFinite(n) && std::isfinite(d);
  }

  // Signed distance (in world units) when `n` is normalized.
  double signedDistance(const Vec3d& p) const {
    return dot(n, p) + d;
  }

  // Normalize plane coefficients so |n| == 1 (when possible).
  void normalize(double eps = 1e-12) {
    if (!isFinite()) return;
    const double lsq = n.lengthSq();
    if (!(lsq > eps * eps)) return;
    const double invL = 1.0 / std::sqrt(lsq);
    n *= invL;
    d *= invL;
  }

  static Plane3d fromPointNormal(const Vec3d& point, const Vec3d& normal) {
    Plane3d p{};
    p.n = safeNormalized(normal, {0.0, 0.0, 1.0}, 1e-18);
    p.d = -dot(p.n, point);
    // n is already normalized by safeNormalized.
    return p;
  }

  static Plane3d fromCoefficients(double a, double b, double c, double dIn) {
    Plane3d p{};
    p.n = {a, b, c};
    p.d = dIn;
    p.normalize();
    return p;
  }
};

enum class FrustumPlane : int {
  Left = 0,
  Right = 1,
  Bottom = 2,
  Top = 3,
  Near = 4,
  Far = 5,
};

// View frustum represented as 6 inward-facing planes.
struct Frustumd {
  std::array<Plane3d, 6> planes{};

  bool isFinite() const {
    for (const auto& p : planes) {
      if (!p.isFinite()) return false;
    }
    return true;
  }

  bool containsPoint(const Vec3d& p, double eps = 0.0) const {
    for (const auto& pl : planes) {
      if (pl.signedDistance(p) < -eps) return false;
    }
    return true;
  }

  bool intersectsSphere(const Vec3d& center, double radius, double eps = 0.0) const {
    if (!stellar::math::isFinite(center) || !std::isfinite(radius)) return false;
    if (radius < 0.0) return false;

    const double r = radius;
    for (const auto& pl : planes) {
      if (pl.signedDistance(center) < -(r + eps)) return false;
    }
    return true;
  }

  // Conservative AABB-frustum intersection using the "positive vertex" test.
  // Returns false only when the AABB is fully outside any plane.
  bool intersectsAabb(const Aabb3d& box, double eps = 0.0) const {
    if (!box.isFinite()) return false;

    for (const auto& pl : planes) {
      // Vertex with maximal signed distance in the direction of the plane normal.
      const Vec3d v{
        (pl.n.x >= 0.0) ? box.max.x : box.min.x,
        (pl.n.y >= 0.0) ? box.max.y : box.min.y,
        (pl.n.z >= 0.0) ? box.max.z : box.min.z,
      };
      if (pl.signedDistance(v) < -eps) return false;
    }

    return true;
  }
};

// Extract frustum planes from a clip matrix (projection * view).
//
// The returned planes are normalized when possible. This helper assumes an
// OpenGL-style clip space where -w <= (x,y,z) <= w.
inline Frustumd frustumFromClipMatrix(const Mat4d& clipFromWorld) {
  Frustumd f{};
  if (!clipFromWorld.isFinite()) return f;

  // clipFromWorld is column-major. Rows are:
  //   r0 = [m0,  m4,  m8,  m12]
  //   r1 = [m1,  m5,  m9,  m13]
  //   r2 = [m2,  m6,  m10, m14]
  //   r3 = [m3,  m7,  m11, m15]
  const double r0x = clipFromWorld.m[0];
  const double r0y = clipFromWorld.m[4];
  const double r0z = clipFromWorld.m[8];
  const double r0w = clipFromWorld.m[12];

  const double r1x = clipFromWorld.m[1];
  const double r1y = clipFromWorld.m[5];
  const double r1z = clipFromWorld.m[9];
  const double r1w = clipFromWorld.m[13];

  const double r2x = clipFromWorld.m[2];
  const double r2y = clipFromWorld.m[6];
  const double r2z = clipFromWorld.m[10];
  const double r2w = clipFromWorld.m[14];

  const double r3x = clipFromWorld.m[3];
  const double r3y = clipFromWorld.m[7];
  const double r3z = clipFromWorld.m[11];
  const double r3w = clipFromWorld.m[15];

  // Standard extraction: plane = row4 +/- row{i} (1-indexed).
  f.planes[static_cast<int>(FrustumPlane::Left)]   = Plane3d::fromCoefficients(r3x + r0x, r3y + r0y, r3z + r0z, r3w + r0w);
  f.planes[static_cast<int>(FrustumPlane::Right)]  = Plane3d::fromCoefficients(r3x - r0x, r3y - r0y, r3z - r0z, r3w - r0w);
  f.planes[static_cast<int>(FrustumPlane::Bottom)] = Plane3d::fromCoefficients(r3x + r1x, r3y + r1y, r3z + r1z, r3w + r1w);
  f.planes[static_cast<int>(FrustumPlane::Top)]    = Plane3d::fromCoefficients(r3x - r1x, r3y - r1y, r3z - r1z, r3w - r1w);
  f.planes[static_cast<int>(FrustumPlane::Near)]   = Plane3d::fromCoefficients(r3x + r2x, r3y + r2y, r3z + r2z, r3w + r2w);
  f.planes[static_cast<int>(FrustumPlane::Far)]    = Plane3d::fromCoefficients(r3x - r2x, r3y - r2y, r3z - r2z, r3w - r2w);

  return f;
}

inline Frustumd frustumFromViewProjection(const Mat4d& view, const Mat4d& projection) {
  return frustumFromClipMatrix(projection * view);
}

} // namespace stellar::math
