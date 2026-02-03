#include "test_harness.h"

#include "stellar/math/Quat.h"

#include <cmath>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

static bool approxVec(const math::Vec3d& a, const math::Vec3d& b, double eps = 1e-12) {
  return approx(a.x, b.x, eps) && approx(a.y, b.y, eps) && approx(a.z, b.z, eps);
}

int test_quat_fromto() {
  int failures = 0;

  // Same direction => identity (or numerically equivalent).
  {
    const math::Quatd q = math::Quatd::fromTo({0.0, 0.0, 1.0}, {0.0, 0.0, 1.0});
    const math::Vec3d v = q.rotate({0.0, 0.0, 1.0});
    CHECK(approxVec(v, {0.0, 0.0, 1.0}, 1e-12));
  }

  // Nearly same direction should still produce the correct tiny rotation
  // (avoid snapping to identity).
  {
    const double theta = 1e-3; // ~0.057 degrees
    const math::Vec3d from{0.0, 0.0, 1.0};
    const math::Vec3d to{0.0, std::sin(theta), std::cos(theta)}; // already unit length

    const math::Quatd q = math::Quatd::fromTo(from, to);
    const math::Vec3d v = q.rotate(from);
    CHECK(approxVec(v, to, 1e-10));
  }

  // Opposite direction => 180-degree flip around some perpendicular axis.
  {
    const math::Quatd q = math::Quatd::fromTo({0.0, 0.0, 1.0}, {0.0, 0.0, -1.0});
    const math::Vec3d v = q.rotate({0.0, 0.0, 1.0});
    CHECK(approxVec(v, {0.0, 0.0, -1.0}, 1e-12));
  }

  // Simple orthogonal case.
  {
    const math::Quatd q = math::Quatd::fromTo({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
    const math::Vec3d v = q.rotate({1.0, 0.0, 0.0});
    CHECK(approxVec(v, {0.0, 1.0, 0.0}, 1e-12));
  }

  // Arbitrary directions.
  {
    const math::Vec3d from{1.0, 2.0, 3.0};
    const math::Vec3d to{-2.0, 1.0, 0.5};

    const math::Quatd q = math::Quatd::fromTo(from, to);
    const math::Vec3d v = q.rotate(from.normalized());
    CHECK(approxVec(v, to.normalized(), 1e-10));
  }

  return failures;
}
