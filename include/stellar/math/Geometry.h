#pragma once

#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::math {

// Axis-aligned bounding box in 3D.
//
// Semantics:
// - min/max are inclusive.
// - This type is intentionally small and header-only, used across sim/render.
struct Aabb3d {
  // Default-constructed AABBs are empty (non-finite), so `expand()` works as expected.
  Vec3d min{
    std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::infinity()
  };
  Vec3d max{
    -std::numeric_limits<double>::infinity(),
    -std::numeric_limits<double>::infinity(),
    -std::numeric_limits<double>::infinity()
  };

  constexpr Aabb3d() = default;
  constexpr Aabb3d(const Vec3d& min_, const Vec3d& max_) : min(min_), max(max_) {}

  static Aabb3d fromMinMax(const Vec3d& a, const Vec3d& b) {
    return Aabb3d{
      {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)},
      {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)},
    };
  }

  static Aabb3d fromCenterExtents(const Vec3d& c, const Vec3d& eIn) {
    const Vec3d e{std::abs(eIn.x), std::abs(eIn.y), std::abs(eIn.z)};
    return Aabb3d{c - e, c + e};
  }

  bool isFinite() const {
    return stellar::math::isFinite(min) && stellar::math::isFinite(max);
  }

  Vec3d center() const { return (min + max) * 0.5; }
  Vec3d size() const { return max - min; }
  Vec3d extents() const { return (max - min) * 0.5; }

  bool contains(const Vec3d& p) const {
    if (!isFinite() || !stellar::math::isFinite(p)) return false;
    return (p.x >= min.x && p.y >= min.y && p.z >= min.z &&
            p.x <= max.x && p.y <= max.y && p.z <= max.z);
  }

  Vec3d clampPoint(const Vec3d& p) const {
    if (!isFinite() || !stellar::math::isFinite(p)) return p;
    return {
      std::clamp(p.x, min.x, max.x),
      std::clamp(p.y, min.y, max.y),
      std::clamp(p.z, min.z, max.z),
    };
  }

  // Squared distance from point to the box (0 for points inside).
  double distanceSqToPoint(const Vec3d& p) const {
    if (!isFinite() || !stellar::math::isFinite(p)) {
      return std::numeric_limits<double>::infinity();
    }

    const auto axisDist = [](double v, double lo, double hi) {
      if (v < lo) return lo - v;
      if (v > hi) return v - hi;
      return 0.0;
    };

    const double dx = axisDist(p.x, min.x, max.x);
    const double dy = axisDist(p.y, min.y, max.y);
    const double dz = axisDist(p.z, min.z, max.z);
    return dx * dx + dy * dy + dz * dz;
  }

  bool intersectsSphere(const Vec3d& c, double r) const {
    if (!isFinite() || !std::isfinite(r)) return false;
    if (r < 0.0) return false;
    return distanceSqToPoint(c) <= r * r;
  }


  bool intersectsAabb(const Aabb3d& b) const {
    if (!isFinite() || !b.isFinite()) return false;
    // Inclusive overlap.
    return !(b.min.x > max.x || b.max.x < min.x ||
             b.min.y > max.y || b.max.y < min.y ||
             b.min.z > max.z || b.max.z < min.z);
  }

  // Ray-AABB intersection (slab method).
  //
  // `dir` does not need to be normalized; distances are returned in world units.
  // If the ray origin lies inside the box, outTEnter is clamped to 0.
  //
  // Returns false when the box only intersects the *backwards* extension of the ray
  // (i.e., when outTExit < 0).
  bool rayIntersectionT(const Vec3d& origin,
                        const Vec3d& dir,
                        double& outTEnter,
                        double& outTExit,
                        double eps = 1e-12) const {
    outTEnter = 0.0;
    outTExit = 0.0;

    if (!isFinite() || !stellar::math::isFinite(origin) || !stellar::math::isFinite(dir)) return false;

    const Vec3d d = safeNormalized(dir, {0.0, 0.0, 0.0}, 1e-18);
    if (!(d.lengthSq() > 1e-18)) return false;

    double tMin = 0.0;
    double tMax = std::numeric_limits<double>::infinity();

    const auto slab = [&](double o, double da, double lo, double hi) {
      if (std::abs(da) <= eps) {
        // Parallel: must be inside slab.
        return (o >= lo && o <= hi);
      }
      const double inv = 1.0 / da;
      double t1 = (lo - o) * inv;
      double t2 = (hi - o) * inv;
      if (t1 > t2) std::swap(t1, t2);
      tMin = std::max(tMin, t1);
      tMax = std::min(tMax, t2);
      return tMax >= tMin;
    };

    if (!slab(origin.x, d.x, min.x, max.x)) return false;
    if (!slab(origin.y, d.y, min.y, max.y)) return false;
    if (!slab(origin.z, d.z, min.z, max.z)) return false;

    // Entire intersection interval behind the ray origin.
    if (tMax < 0.0) return false;

    outTEnter = std::max(0.0, tMin);
    outTExit = tMax;
    return true;
  }

  // Segment-AABB intersection (slab method).
  //
  // Returns true when the segment [a,b] intersects the box. If so:
  //  - outTEnter/outTExit are the normalized segment parameters in [0,1]
  //    for the entry/exit points (clamped to the segment domain).
  //
  // If the segment starts inside the box, outTEnter is clamped to 0.
  bool segmentIntersectionT(const Vec3d& a,
                            const Vec3d& b,
                            double& outTEnter,
                            double& outTExit,
                            double eps = 1e-12) const {
    outTEnter = 0.0;
    outTExit = 0.0;

    if (!isFinite() || !stellar::math::isFinite(a) || !stellar::math::isFinite(b)) return false;

    const Vec3d d = b - a;

    double tMin = 0.0;
    double tMax = 1.0;

    const auto slab = [&](double o, double da, double lo, double hi) {
      if (std::abs(da) <= eps) {
        // Parallel: must be inside slab.
        return (o >= lo && o <= hi);
      }
      const double inv = 1.0 / da;
      double t1 = (lo - o) * inv;
      double t2 = (hi - o) * inv;
      if (t1 > t2) std::swap(t1, t2);
      tMin = std::max(tMin, t1);
      tMax = std::min(tMax, t2);
      return tMax >= tMin;
    };

    if (!slab(a.x, d.x, min.x, max.x)) return false;
    if (!slab(a.y, d.y, min.y, max.y)) return false;
    if (!slab(a.z, d.z, min.z, max.z)) return false;

    // Reject if the intersection is entirely outside the segment domain.
    if (tMax < 0.0 || tMin > 1.0) return false;

    outTEnter = std::clamp(tMin, 0.0, 1.0);
    outTExit = std::clamp(tMax, 0.0, 1.0);
    return true;
  }

  void expand(const Vec3d& p) {
    if (!stellar::math::isFinite(p)) return;

    if (!isFinite()) {
      min = max = p;
      return;
    }

    min.x = std::min(min.x, p.x);
    min.y = std::min(min.y, p.y);
    min.z = std::min(min.z, p.z);

    max.x = std::max(max.x, p.x);
    max.y = std::max(max.y, p.y);
    max.z = std::max(max.z, p.z);
  }

  void expand(const Aabb3d& b) {
    if (!b.isFinite()) return;
    expand(b.min);
    expand(b.max);
  }
};

// Free-function helper to mirror common vector-math APIs.
inline double distanceSqToPoint(const Aabb3d& box, const Vec3d& p) {
  return box.distanceSqToPoint(p);
}

// Signed distance of a point to a plane.
//
// Plane is defined by a point on the plane and a (not necessarily unit-length) normal.
// If `planeNormal` is unit length, this returns a true signed distance in world units.
inline double planeSignedDistance(const Vec3d& p,
                                  const Vec3d& planePoint,
                                  const Vec3d& planeNormal) {
  if (!isFinite(p) || !isFinite(planePoint) || !isFinite(planeNormal)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return dot(p - planePoint, planeNormal);
}

// Segment-plane intersection.
//
// Returns true if the segment [a,b] intersects the plane defined by (planePoint, planeNormal).
// If so, `outT` is set to the normalized segment parameter in [0,1] where:
//   p = a + (b - a) * outT
//
// Notes:
// - The plane normal does not need to be normalized.
// - If the segment is (nearly) parallel to the plane, this returns false unless `a`
//   lies on the plane (in which case outT=0 and true is returned).
inline bool segmentPlaneIntersectionT(const Vec3d& a,
                                      const Vec3d& b,
                                      const Vec3d& planePoint,
                                      const Vec3d& planeNormal,
                                      double& outT,
                                      double eps = 1e-12) {
  outT = 0.0;

  if (!isFinite(a) || !isFinite(b) || !isFinite(planePoint) || !isFinite(planeNormal)) return false;
  if (!(planeNormal.lengthSq() > 1e-24)) return false;

  const Vec3d ab = b - a;
  const double denom = dot(ab, planeNormal);
  const double num = dot(planePoint - a, planeNormal);

  if (std::abs(denom) <= eps) {
    // Parallel: treat as hit only if the start point lies on the plane.
    if (std::abs(num) <= eps) {
      outT = 0.0;
      return true;
    }
    return false;
  }

  double t = num / denom;

  // Accept a small tolerance and clamp to [0,1] for numeric robustness.
  const double tEps = 1e-12;
  if (t < -tEps || t > 1.0 + tEps) return false;
  t = std::clamp(t, 0.0, 1.0);

  outT = t;
  return true;
}

// Convenience overload that also returns the intersection point.
inline bool segmentPlaneIntersection(const Vec3d& a,
                                     const Vec3d& b,
                                     const Vec3d& planePoint,
                                     const Vec3d& planeNormal,
                                     double& outT,
                                     Vec3d& outPoint,
                                     double eps = 1e-12) {
  if (!segmentPlaneIntersectionT(a, b, planePoint, planeNormal, outT, eps)) return false;
  outPoint = a + (b - a) * outT;
  return true;
}

// Return the clamped parametric t in [0,1] for the closest point on segment AB to point P.
//
// If the segment is degenerate (A == B), returns 0.
inline double segmentClosestT(const Vec3d& a, const Vec3d& b, const Vec3d& p) {
  if (!isFinite(a) || !isFinite(b) || !isFinite(p)) return 0.0;
  const Vec3d ab = b - a;
  const double abLenSq = ab.lengthSq();
  if (!(abLenSq > 1e-18)) return 0.0;
  return std::clamp(dot(p - a, ab) / abLenSq, 0.0, 1.0);
}

// Segment-sphere hit test.
//
// Returns true if any point on segment [A,B] lies within radius R of center C.
inline bool segmentHitsSphere(const Vec3d& a,
                              const Vec3d& b,
                              const Vec3d& center,
                              double radius) {
  if (!isFinite(a) || !isFinite(b) || !isFinite(center) || !std::isfinite(radius)) return false;
  if (radius < 0.0) return false;

  const Vec3d ab = b - a;
  const double abLenSq = ab.lengthSq();
  const double r2 = radius * radius;

  if (!(abLenSq > 1e-18)) {
    return (a - center).lengthSq() <= r2;
  }

  const double t = segmentClosestT(a, b, center);
  const Vec3d closest = a + ab * t;
  return (closest - center).lengthSq() <= r2;
}

// Segment-sphere intersection.
//
// Returns true if the segment [A,B] intersects the sphere. If so:
//  - outTEnter/outTExit are the normalized segment parameters (in [0,1]) of the
//    entry and exit points along the segment (clamped to the segment domain).
//
// Special cases:
//  - If A == B, this reduces to a point-in-sphere test with outTEnter=outTExit=0.
//  - If the segment starts inside the sphere, outTEnter is clamped to 0.
inline bool segmentSphereIntersectionT(const Vec3d& a,
                                       const Vec3d& b,
                                       const Vec3d& center,
                                       double radius,
                                       double& outTEnter,
                                       double& outTExit) {
  outTEnter = 0.0;
  outTExit = 0.0;

  if (!isFinite(a) || !isFinite(b) || !isFinite(center) || !std::isfinite(radius)) return false;
  if (radius < 0.0) return false;

  const Vec3d d = b - a;
  const double dd = d.lengthSq();
  const double r2 = radius * radius;

  // Degenerate segment: treat as point.
  if (!(dd > 1e-18)) {
    const double dist2 = (a - center).lengthSq();
    if (dist2 <= r2) {
      outTEnter = 0.0;
      outTExit = 0.0;
      return true;
    }
    return false;
  }

  // Solve: |(a + d*t) - center|^2 = r^2 for t in [0,1].
  const Vec3d m = a - center;

  // Quadratic coefficients:
  //   dd * t^2 + 2*(m·d) * t + (m·m - r^2) = 0
  const double bq = dot(m, d);
  const double cq = m.lengthSq() - r2;

  const double disc = bq * bq - dd * cq;
  if (disc < 0.0) return false;

  const double sdisc = std::sqrt(std::max(0.0, disc));
  double t0 = (-bq - sdisc) / dd;
  double t1 = (-bq + sdisc) / dd;
  if (t0 > t1) std::swap(t0, t1);

  // Reject if the entire intersection interval is outside the segment.
  if (t1 < 0.0 || t0 > 1.0) return false;

  outTEnter = std::clamp(t0, 0.0, 1.0);
  outTExit = std::clamp(t1, 0.0, 1.0);
  return true;
}

// Convenience overload returning only the entry parameter.
inline bool segmentSphereIntersectionEnterT(const Vec3d& a,
                                            const Vec3d& b,
                                            const Vec3d& center,
                                            double radius,
                                            double& outTEnter) {
  double tExit = 0.0;
  return segmentSphereIntersectionT(a, b, center, radius, outTEnter, tExit);
}

// Ray-sphere intersection.
//
// `dir` does not need to be normalized; distances are returned in world units.
// If the origin is inside the sphere, outTEnter is clamped to 0.
//
// Returns false when the sphere only intersects the *backwards* extension of the ray
// (i.e., when both intersection points are behind the origin).
inline bool raySphereIntersect(const Vec3d& origin,
                               const Vec3d& dir,
                               const Vec3d& center,
                               double radius,
                               double& outTEnter,
                               double& outTExit) {
  if (!isFinite(origin) || !isFinite(center) || !isFinite(dir) || !std::isfinite(radius)) return false;
  if (radius < 0.0) return false;

  const Vec3d d = safeNormalized(dir, {0.0, 0.0, 0.0}, 1e-18);
  if (!(d.lengthSq() > 1e-18)) return false;

  const Vec3d oc = center - origin;
  const double tProj = dot(oc, d);
  const double dist2 = oc.lengthSq();

  const double r2 = radius * radius;
  const double d2 = dist2 - tProj * tProj;
  if (d2 > r2) return false;

  const double thc = std::sqrt(std::max(0.0, r2 - d2));
  double tEnter = tProj - thc;
  double tExit = tProj + thc;

  // If the entire sphere is behind the ray origin, there is no forward hit.
  if (tExit < 0.0) return false;

  // If we start inside the sphere, clamp entry to 0.
  if (tEnter < 0.0) tEnter = 0.0;

  outTEnter = tEnter;
  outTExit = tExit;
  return true;
}

// Convenience overload returning only the entry distance.
inline bool raySphereIntersectEnter(const Vec3d& origin,
                                    const Vec3d& dir,
                                    const Vec3d& center,
                                    double radius,
                                    double& outTEnter) {
  double tExit = 0.0;
  return raySphereIntersect(origin, dir, center, radius, outTEnter, tExit);
}

} // namespace stellar::math
