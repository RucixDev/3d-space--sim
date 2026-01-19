#pragma once

#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace stellar::math {

// -----------------------------------------------------------------------------
// Catmull-Rom splines + arc-length reparameterization.
//
// These helpers are intentionally minimal (dependency-free) and tuned for
// animation/tooling use-cases (camera rigs, debug ghost playback, UI graphs).
// They are NOT intended as a fully-featured animation system.
// -----------------------------------------------------------------------------

// Uniform Catmull-Rom spline (u in [0,1]).
// Returns p1 at u=0 and p2 at u=1.
inline Vec3d catmullRomUniform(const Vec3d& p0,
                              const Vec3d& p1,
                              const Vec3d& p2,
                              const Vec3d& p3,
                              double u) {
  u = std::clamp(u, 0.0, 1.0);
  const double u2 = u * u;
  const double u3 = u2 * u;
  return (p1 * 2.0 + (p2 - p0) * u + (p0 * 2.0 - p1 * 5.0 + p2 * 4.0 - p3) * u2 +
          (-p0 + p1 * 3.0 - p2 * 3.0 + p3) * u3) *
         0.5;
}

inline Vec3d lerp(const Vec3d& a, const Vec3d& b, double t) {
  return a * (1.0 - t) + b * t;
}

// Lerp between points a/b defined at parameter times ta/tb evaluated at time t.
inline Vec3d lerpT(const Vec3d& a, const Vec3d& b, double ta, double tb, double t) {
  const double denom = tb - ta;
  if (std::abs(denom) < 1e-12) return a;
  const double u = (t - ta) / denom;
  return lerp(a, b, u);
}

// Centripetal Catmull-Rom (alpha=0.5 by default).
//
// Compared to uniform Catmull-Rom, the centripetal parameterization tends to
// reduce overshoot and self-intersections around sharp corners.
//
// Implementation notes:
//  - Uses the Barry & Goldman formulation (recursive lerps in parameter t-space).
//  - Falls back to uniform Catmull-Rom when the knot vector degenerates.
inline Vec3d catmullRomCentripetal(const Vec3d& p0,
                                  const Vec3d& p1,
                                  const Vec3d& p2,
                                  const Vec3d& p3,
                                  double u01,
                                  double alpha = 0.5) {
  u01 = std::clamp(u01, 0.0, 1.0);

  auto tj = [&](double ti, const Vec3d& a, const Vec3d& b) -> double {
    const double d = (b - a).length();
    // epsilon so repeated points don't trigger divide-by-zero in lerpT.
    const double v = std::pow(std::max(d, 1e-9), alpha);
    return ti + v;
  };

  const double t0 = 0.0;
  const double t1 = tj(t0, p0, p1);
  const double t2 = tj(t1, p1, p2);
  const double t3 = tj(t2, p2, p3);

  // Degenerate fallback.
  if (!(t1 > t0 + 1e-12 && t2 > t1 + 1e-12 && t3 > t2 + 1e-12)) {
    return catmullRomUniform(p0, p1, p2, p3, u01);
  }

  const double t = t1 + u01 * (t2 - t1);

  const Vec3d A1 = lerpT(p0, p1, t0, t1, t);
  const Vec3d A2 = lerpT(p1, p2, t1, t2, t);
  const Vec3d A3 = lerpT(p2, p3, t2, t3, t);

  const Vec3d B1 = lerpT(A1, A2, t0, t2, t);
  const Vec3d B2 = lerpT(A2, A3, t1, t3, t);

  const Vec3d C = lerpT(B1, B2, t1, t2, t);
  return C;
}

// -----------------------------------------------------------------------------
// Arc-length reparameterization.
//
// Given a curve evaluator f(u) with u in [0,1], we approximate its arc-length by
// sampling and build a table to invert s->u.
//
// This is useful when you want a constant-speed motion along a spline: drive a
// uniform u over time, then map it through reparamByArcLength().
// -----------------------------------------------------------------------------

struct ArcLengthTable {
  int samples{0};
  std::vector<double> s; // cumulative arc length; size = samples+1
  double total{0.0};
};

template <class Eval>
inline ArcLengthTable buildArcLengthTable(Eval&& eval, int samples = 32) {
  ArcLengthTable tbl{};
  tbl.samples = std::clamp(samples, 2, 512);
  tbl.s.resize(static_cast<std::size_t>(tbl.samples) + 1u);
  tbl.s[0] = 0.0;

  Vec3d prev = eval(0.0);
  double acc = 0.0;

  for (int i = 1; i <= tbl.samples; ++i) {
    const double u = static_cast<double>(i) / static_cast<double>(tbl.samples);
    const Vec3d p = eval(u);
    acc += (p - prev).length();
    tbl.s[static_cast<std::size_t>(i)] = acc;
    prev = p;
  }

  tbl.total = acc;
  return tbl;
}

// Invert arc length: given s in [0,total], return u in [0,1].
inline double invertArcLength(const ArcLengthTable& tbl, double s) {
  if (tbl.samples < 2 || tbl.s.empty() || !(tbl.total > 1e-12)) return 0.0;

  s = std::clamp(s, 0.0, tbl.total);

  // Find first element > s.
  auto it = std::upper_bound(tbl.s.begin(), tbl.s.end(), s);
  int hi = static_cast<int>(it - tbl.s.begin());
  hi = std::clamp(hi, 1, tbl.samples);

  const int lo = hi - 1;
  const double s0 = tbl.s[static_cast<std::size_t>(lo)];
  const double s1 = tbl.s[static_cast<std::size_t>(hi)];
  const double denom = std::max(1e-12, s1 - s0);
  const double t = (s - s0) / denom;

  const double u = (static_cast<double>(lo) + t) / static_cast<double>(tbl.samples);
  return std::clamp(u, 0.0, 1.0);
}

inline double reparamByArcLength(const ArcLengthTable& tbl, double u01) {
  u01 = std::clamp(u01, 0.0, 1.0);
  if (!(tbl.total > 1e-12)) return u01;
  return invertArcLength(tbl, u01 * tbl.total);
}

} // namespace stellar::math
