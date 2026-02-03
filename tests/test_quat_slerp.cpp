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

int test_quat_slerp() {
  int failures = 0;

  const math::Quatd a = math::Quatd::identity();
  const math::Quatd b = math::Quatd::fromAxisAngle({0.0, 1.0, 0.0}, math::kPi * 0.5);

  // Midpoint should be a 45-degree yaw about +Y.
  {
    const math::Quatd mid = math::Quatd::slerp(a, b, 0.5);
    const math::Vec3d f = mid.rotate({0.0, 0.0, 1.0});
    const math::Vec3d expected{std::sin(math::kPi * 0.25), 0.0, std::cos(math::kPi * 0.25)};
    CHECK(approxVec(f, expected, 1e-12));
  }

  // Negating b represents the same rotation; slerp should still take the shortest arc.
  {
    const math::Quatd bNeg{-b.w, -b.x, -b.y, -b.z};
    const math::Quatd mid = math::Quatd::slerp(a, bNeg, 0.5);
    const math::Vec3d f = mid.rotate({0.0, 0.0, 1.0});
    const math::Vec3d expected{std::sin(math::kPi * 0.25), 0.0, std::cos(math::kPi * 0.25)};
    CHECK(approxVec(f, expected, 1e-12));
  }

  // Endpoints.
  {
    const math::Quatd q0 = math::Quatd::slerp(a, b, 0.0);
    const math::Quatd q1 = math::Quatd::slerp(a, b, 1.0);
    CHECK(approxVec(q0.rotate({0.0, 0.0, 1.0}), a.rotate({0.0, 0.0, 1.0}), 1e-12));
    CHECK(approxVec(q1.rotate({0.0, 0.0, 1.0}), b.rotate({0.0, 0.0, 1.0}), 1e-12));
  }

  return failures;
}
