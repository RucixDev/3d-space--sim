#include "stellar/math/Mat4.h"

#include <cmath>
#include <iostream>

static bool approx(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

static bool finiteMat(const stellar::math::Mat4d& m) {
  return m.isFinite();
}

static bool nearIdentity(const stellar::math::Mat4d& m, double eps) {
  for (int i = 0; i < 16; ++i) {
    const double expected = (i == 0 || i == 5 || i == 10 || i == 15) ? 1.0 : 0.0;
    if (!approx(m.m[i], expected, eps)) return false;
  }
  return true;
}

int test_mat4_inverse() {
  int fails = 0;

  // Pure translation.
  {
    const stellar::math::Mat4d t = stellar::math::Mat4d::translation({1.0, 2.0, 3.0});
    const auto inv = t.inverseRigid();
    const auto id = t * inv;
    if (!finiteMat(inv) || !finiteMat(id) || !nearIdentity(id, 1e-12)) {
      std::cerr << "[test_mat4_inverse] translation: inverseRigid failed\n";
      ++fails;
    }

    if (!approx(inv.m[12], -1.0, 1e-12) || !approx(inv.m[13], -2.0, 1e-12) || !approx(inv.m[14], -3.0, 1e-12)) {
      std::cerr << "[test_mat4_inverse] translation: expected -t in inverse\n";
      ++fails;
    }
  }

  // Rotation + translation.
  {
    const auto q = stellar::math::Quatd::fromAxisAngle({0.2, 1.0, 0.1}, stellar::math::degToRad(35.0));
    const auto m = stellar::math::Mat4d::translation({-4.0, 2.0, 1.0}) * stellar::math::Mat4d::rotation(q);
    const auto inv = m.inverseRigid();
    const auto id = m * inv;

    if (!finiteMat(inv) || !finiteMat(id) || !nearIdentity(id, 1e-10)) {
      std::cerr << "[test_mat4_inverse] rigid: m * inv != I\n";
      ++fails;
    }

    // Round-trip a point using homogeneous coordinates.
    const stellar::math::Vec3d p{3.0, -1.0, 7.0};
    const stellar::math::Vec4d wp = m.mulPoint(p, 1.0);
    const stellar::math::Vec4d lp = inv.mul(wp);
    const double eps = 1e-9;
    if (!approx(lp.w, 1.0, eps) ||
        !approx(lp.x, p.x, eps) || !approx(lp.y, p.y, eps) || !approx(lp.z, p.z, eps)) {
      std::cerr << "[test_mat4_inverse] rigid: point round-trip failed\n";
      ++fails;
    }
  }

  // For a pure rotation matrix, inverse == transpose.
  {
    const auto q = stellar::math::Quatd::fromEulerYXZ(stellar::math::degToRad(20.0),
                                                      stellar::math::degToRad(-15.0),
                                                      stellar::math::degToRad(7.0));
    const auto r = stellar::math::Mat4d::rotation(q);
    const auto rt = r.transposed();
    const auto id = r * rt;
    if (!nearIdentity(id, 1e-10)) {
      std::cerr << "[test_mat4_inverse] rotation: R*R^T != I\n";
      ++fails;
    }
  }

  return fails;
}
