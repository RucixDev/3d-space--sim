#pragma once

#include <cmath>

namespace stellar::math {

struct Vec4d {
  double x{0}, y{0}, z{0}, w{0};

  constexpr Vec4d() = default;
  constexpr Vec4d(double x_, double y_, double z_, double w_) : x(x_), y(y_), z(z_), w(w_) {}

  Vec4d operator+(const Vec4d& o) const { return {x + o.x, y + o.y, z + o.z, w + o.w}; }
  Vec4d operator-(const Vec4d& o) const { return {x - o.x, y - o.y, z - o.z, w - o.w}; }
  Vec4d operator*(double s) const { return {x * s, y * s, z * s, w * s}; }
  Vec4d operator/(double s) const { return {x / s, y / s, z / s, w / s}; }

  Vec4d& operator+=(const Vec4d& o) {
    x += o.x;
    y += o.y;
    z += o.z;
    w += o.w;
    return *this;
  }
  Vec4d& operator-=(const Vec4d& o) {
    x -= o.x;
    y -= o.y;
    z -= o.z;
    w -= o.w;
    return *this;
  }
  Vec4d& operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    w *= s;
    return *this;
  }
  Vec4d& operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    w /= s;
    return *this;
  }

  double lengthSq() const { return x * x + y * y + z * z + w * w; }
  double length() const { return std::sqrt(lengthSq()); }
};

inline Vec4d operator*(double s, const Vec4d& v) { return v * s; }

inline double dot(const Vec4d& a, const Vec4d& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

} // namespace stellar::math
