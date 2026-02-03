#pragma once

#include <algorithm>
#include <cmath>

namespace stellar::math {

constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kLn2 = 0.693147180559945309417232121458176568;

inline double degToRad(double deg) { return deg * (kPi / 180.0); }
inline double radToDeg(double rad) { return rad * (180.0 / kPi); }

template <class T>
inline T clamp(T v, T lo, T hi) {
  return std::max(lo, std::min(v, hi));
}

template <class T>
inline T lerp(T a, T b, double t) {
  return static_cast<T>(a + (b - a) * t);
}

// -----------------------------------------------------------------------------
// Half-life smoothing helpers
// -----------------------------------------------------------------------------
//
// We use a half-life parameterization because it is stable across dt and easy
// to reason about:
//   dt == halfLife -> value moves halfway toward the target.
//
// decayFactor: multiply a value by this to apply exponential decay toward 0.
// alpha:        interpolation weight for y += alpha * (x - y).

// Exponential decay factor in [0,1].
//
// - dt <= 0      => 1 (no decay)
// - halfLife<=0  => 0 for dt>0 (instant decay)
inline double halfLifeDecayFactor(double dt, double halfLife) {
  if (dt <= 0.0) return 1.0;
  if (halfLife <= 0.0) return 0.0;
  return std::exp(-kLn2 * dt / halfLife);
}

// Interpolation alpha in [0,1] for half-life exponential smoothing.
//
// Equivalent to:
//   y = lerp(y, x, alpha)
// and reaches half the gap after `halfLife` time.
inline double halfLifeAlpha(double dt, double halfLife) {
  return 1.0 - halfLifeDecayFactor(dt, halfLife);
}

} // namespace stellar::math
