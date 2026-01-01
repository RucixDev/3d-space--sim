#pragma once
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

  Vec3d& operator+=(const Vec3d& o) { x += o.x; y += o.y; z += o.z; return *this; }
  Vec3d& operator-=(const Vec3d& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
  Vec3d& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
  Vec3d& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }

  double lengthSq() const { return x*x + y*y + z*z; }
  double length() const { return std::sqrt(lengthSq()); }

  Vec3d normalized(double eps=1e-12) const {
    const double len = length();
    if (len < eps) return {0,0,0};
    return *this / len;
  }
};

// Allow scalar * vector (Vec3d already supports vector * scalar).
inline Vec3d operator*(double s, const Vec3d& v) { return v * s; }

inline double dot(const Vec3d& a, const Vec3d& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline Vec3d cross(const Vec3d& a, const Vec3d& b) {
  return {
    a.y*b.z - a.z*b.y,
    a.z*b.x - a.x*b.z,
    a.x*b.y - a.y*b.x
  };
}

} // namespace stellar::math
