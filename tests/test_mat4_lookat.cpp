#include "stellar/math/Mat4.h"

#include <cmath>
#include <iostream>

static bool finiteMat(const stellar::math::Mat4d& m) {
  for (double v : m.m) {
    if (!std::isfinite(v)) return false;
  }
  return true;
}

static bool approx(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

static stellar::math::Vec3d col(const stellar::math::Mat4d& m, int c) {
  // Column-major.
  return {m.m[c * 4 + 0], m.m[c * 4 + 1], m.m[c * 4 + 2]};
}

int test_mat4_lookat() {
  int fails = 0;

  // Canonical camera: eye at origin looking down -Z with +Y up.
  {
    const auto v = stellar::math::Mat4d::lookAt({0, 0, 0}, {0, 0, -1}, {0, 1, 0});
    if (!finiteMat(v)) {
      std::cerr << "[test_mat4_lookat] canonical: matrix not finite\n";
      ++fails;
    }

    // Should be identity.
    const double eps = 1e-12;
    for (int i = 0; i < 16; ++i) {
      const double expected = (i == 0 || i == 5 || i == 10 || i == 15) ? 1.0 : 0.0;
      if (!approx(v.m[i], expected, eps)) {
        std::cerr << "[test_mat4_lookat] canonical: expected identity at i=" << i << "\n";
        ++fails;
        break;
      }
    }
  }

  // Degenerate forward vector: eye == center should still produce a valid view.
  {
    const stellar::math::Vec3d eye{1.0, 2.0, 3.0};
    const auto v = stellar::math::Mat4d::lookAt(eye, eye, {0, 1, 0});

    if (!finiteMat(v)) {
      std::cerr << "[test_mat4_lookat] eye==center: matrix not finite\n";
      ++fails;
    }

    const double eps = 1e-12;
    // Rotation part should match default look down -Z (i.e. identity basis).
    if (!approx(v.m[0], 1.0, eps) || !approx(v.m[5], 1.0, eps) || !approx(v.m[10], 1.0, eps)) {
      std::cerr << "[test_mat4_lookat] eye==center: expected default basis\n";
      ++fails;
    }

    // Translation should move the world by -eye.
    if (!approx(v.m[12], -eye.x, eps) || !approx(v.m[13], -eye.y, eps) || !approx(v.m[14], -eye.z, eps)) {
      std::cerr << "[test_mat4_lookat] eye==center: expected translation -eye\n";
      ++fails;
    }
  }

  // Degenerate up vector: up parallel to forward should still produce an orthonormal frame.
  {
    const auto v = stellar::math::Mat4d::lookAt({0, 0, 0}, {0, 1, 0}, {0, 1, 0});
    if (!finiteMat(v)) {
      std::cerr << "[test_mat4_lookat] up||forward: matrix not finite\n";
      ++fails;
    } else {
      const auto c0 = col(v, 0);
      const auto c1 = col(v, 1);
      const auto c2 = col(v, 2);

      const double eps = 1e-6;
      if (!approx(c0.length(), 1.0, eps) || !approx(c1.length(), 1.0, eps) || !approx(c2.length(), 1.0, eps)) {
        std::cerr << "[test_mat4_lookat] up||forward: basis not unit length\n";
        ++fails;
      }

      if (std::abs(stellar::math::dot(c0, c1)) > 1e-5 ||
          std::abs(stellar::math::dot(c0, c2)) > 1e-5 ||
          std::abs(stellar::math::dot(c1, c2)) > 1e-5) {
        std::cerr << "[test_mat4_lookat] up||forward: basis not orthogonal\n";
        ++fails;
      }
    }
  }

  return fails;
}
