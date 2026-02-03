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

int test_quat_rotate() {
  int failures = 0;

  // rotate() should be scale-invariant: q and (s*q) represent the same rotation.
  {
    const math::Quatd q = math::Quatd::fromAxisAngle({1.0, 2.0, 3.0}, 0.7);
    const math::Quatd qs{q.w * 3.0, q.x * 3.0, q.y * 3.0, q.z * 3.0};

    const math::Vec3d v{0.3, -0.5, 1.2};

    const math::Vec3d a = q.rotate(v);
    const math::Vec3d b = qs.rotate(v);

    CHECK(approxVec(a, b, 1e-12));
  }

  // fromAxisAngle should be defensive for degenerate axes.
  {
    const math::Quatd q = math::Quatd::fromAxisAngle({0.0, 0.0, 0.0}, 1.0);
    const math::Vec3d v{1.0, 2.0, 3.0};
    CHECK(approxVec(q.rotate(v), v, 1e-12));
  }

  return failures;
}
