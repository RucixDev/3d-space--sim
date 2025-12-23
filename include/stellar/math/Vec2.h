#pragma once
#include <cmath>

namespace stellar::math {

struct Vec2d {
  double x{0}, y{0};

  constexpr Vec2d() = default;
  constexpr Vec2d(double x_, double y_) : x(x_), y(y_) {}

  Vec2d operator+(const Vec2d& o) const { return {x + o.x, y + o.y}; }
  Vec2d operator-(const Vec2d& o) const { return {x - o.x, y - o.y}; }
  Vec2d operator*(double s) const { return {x * s, y * s}; }
  Vec2d operator/(double s) const { return {x / s, y / s}; }

  Vec2d& operator+=(const Vec2d& o) { x += o.x; y += o.y; return *this; }
  Vec2d& operator-=(const Vec2d& o) { x -= o.x; y -= o.y; return *this; }
  Vec2d& operator*=(double s) { x *= s; y *= s; return *this; }
  Vec2d& operator/=(double s) { x /= s; y /= s; return *this; }

  double length() const { return std::sqrt(x*x + y*y); }
};

inline double dot(const Vec2d& a, const Vec2d& b) { return a.x*b.x + a.y*b.y; }

} // namespace stellar::math
