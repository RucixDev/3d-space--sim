#pragma once

#include <algorithm>
#include <cmath>

namespace stellar::math {

constexpr double kPi = 3.1415926535897932384626433832795;

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

} // namespace stellar::math
