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

int test_quat_integrate() {
  int failures = 0;

  // integrateAngular() should match an axis-angle rotation for constant body-frame omega,
  // and be dt-invariant (many small steps == one larger step).

  // ---- Known rotation: yaw about +Y ----
  {
    const math::Quatd q0 = math::Quatd::identity();
    const math::Vec3d omega{0.0, 1.0, 0.0}; // rad/s about +Y
    const double dt = 1.0;

    const math::Quatd q1 = q0.integrateAngular(omega, dt);

    const math::Vec3d f = q1.rotate({0.0, 0.0, 1.0});
    const math::Vec3d expected{std::sin(1.0), 0.0, std::cos(1.0)};

    CHECK(approxVec(f, expected, 1e-12));
  }

  // ---- dt invariance for arbitrary omega ----
  {
    const math::Quatd q0 = math::Quatd::identity();
    const math::Vec3d omega{0.2, -0.4, 1.0};

    const math::Quatd a = q0.integrateAngular(omega, 0.1);

    math::Quatd b = q0;
    for (int i = 0; i < 10; ++i) {
      b = b.integrateAngular(omega, 0.01);
    }

    const math::Vec3d fA = a.rotate({0.0, 0.0, 1.0});
    const math::Vec3d fB = b.rotate({0.0, 0.0, 1.0});

    CHECK(approxVec(fA, fB, 1e-12));
  }

  return failures;
}
