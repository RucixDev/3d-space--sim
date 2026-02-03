#pragma once
#include <algorithm>
#include <cmath>

namespace stellar::math {

struct Vec3d {
  double x{0}, y{0}, z{0};

  constexpr Vec3d() = default;
  constexpr Vec3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  Vec3d operator+(const Vec3d& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3d operator-(const Vec3d& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3d operator-() const { return {-x, -y, -z}; }
  Vec3d operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3d operator/(double s) const { return {x / s, y / s, z / s}; }

  Vec3d& operator+=(const Vec3d& o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  Vec3d& operator-=(const Vec3d& o) {
    x -= o.x;
    y -= o.y;
    z -= o.z;
    return *this;
  }
  Vec3d& operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
  }
  Vec3d& operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
  }

  double lengthSq() const { return x * x + y * y + z * z; }
  double length() const { return std::sqrt(lengthSq()); }

  Vec3d normalized(double eps = 1e-12) const {
    const double len = length();
    if (len < eps) return {0, 0, 0};
    return *this / len;
  }
};

// Allow scalar * vector (Vec3d already supports vector * scalar).
inline Vec3d operator*(double s, const Vec3d& v) { return v * s; }

inline bool isFinite(const Vec3d& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

inline double dot(const Vec3d& a, const Vec3d& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3d cross(const Vec3d& a, const Vec3d& b) {
  return {
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x,
  };
}

// Free-function helpers (mirrors common vector-math APIs).
inline double length(const Vec3d& v) { return v.length(); }

// Normalize a vector, returning `fallback` if the input is too small or non-finite.
//
// epsSq is a squared-length threshold.
inline Vec3d safeNormalized(const Vec3d& v, const Vec3d& fallback, double epsSq = 1e-18) {
  if (!isFinite(v)) return fallback;
  const double lsq = v.lengthSq();
  if (!(lsq > epsSq)) return fallback;
  return v * (1.0 / std::sqrt(lsq));
}

inline Vec3d clampComponents(const Vec3d& v, double lo, double hi) {
  return {
    std::clamp(v.x, lo, hi),
    std::clamp(v.y, lo, hi),
    std::clamp(v.z, lo, hi),
  };
}

// Clamp a vector's magnitude without changing direction.
//
// If the vector is near zero (<=eps), it is returned unchanged.
inline Vec3d clampMagnitude(const Vec3d& v, double maxLen, double eps = 1e-12) {
  if (!isFinite(v) || !std::isfinite(maxLen)) return v;
  if (maxLen <= 0.0) return {0, 0, 0};

  const double lsq = v.lengthSq();
  if (!(lsq > eps * eps)) return v;

  const double maxSq = maxLen * maxLen;
  if (lsq <= maxSq) return v;

  return v * (maxLen / std::sqrt(lsq));
}

} // namespace stellar::math
