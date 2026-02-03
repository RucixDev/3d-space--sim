#include "stellar/math/Vec3.h"

#include <cmath>
#include <iostream>
#include <limits>

using namespace stellar;

static bool approxEq(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

static bool approxEqVec(const math::Vec3d& a, const math::Vec3d& b, double eps = 1e-12) {
  return approxEq(a.x, b.x, eps) && approxEq(a.y, b.y, eps) && approxEq(a.z, b.z, eps);
}

int test_vec3_utils() {
  int fails = 0;

  // safeNormalized: near-zero returns fallback
  {
    const math::Vec3d v{0, 0, 0};
    const math::Vec3d fb{1, 0, 0};
    const math::Vec3d n = math::safeNormalized(v, fb);
    if (!approxEqVec(n, fb)) {
      std::cerr << "[test_vec3_utils] safeNormalized(zero) expected fallback\n";
      ++fails;
    }
  }

  // safeNormalized: normal case returns unit length and direction
  {
    const math::Vec3d v{3, 4, 0};
    const math::Vec3d fb{1, 0, 0};
    const math::Vec3d n = math::safeNormalized(v, fb);
    if (!approxEq(n.length(), 1.0, 1e-12) || !approxEqVec(n, {0.6, 0.8, 0.0}, 1e-12)) {
      std::cerr << "[test_vec3_utils] safeNormalized(3,4,0) expected (0.6,0.8,0)\n";
      ++fails;
    }
  }

  // safeNormalized: non-finite returns fallback
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const math::Vec3d v{nan, 0.0, 1.0};
    const math::Vec3d fb{0, 1, 0};
    const math::Vec3d n = math::safeNormalized(v, fb);
    if (!approxEqVec(n, fb)) {
      std::cerr << "[test_vec3_utils] safeNormalized(NaN) expected fallback\n";
      ++fails;
    }
  }

  // clampComponents
  {
    const math::Vec3d v{2.0, -3.0, 0.5};
    const math::Vec3d c = math::clampComponents(v, -1.0, 1.0);
    if (!approxEqVec(c, {1.0, -1.0, 0.5})) {
      std::cerr << "[test_vec3_utils] clampComponents expected (1,-1,0.5)\n";
      ++fails;
    }
  }

  // clampMagnitude: already within limit
  {
    const math::Vec3d v{3.0, 4.0, 0.0}; // len=5
    const math::Vec3d c = math::clampMagnitude(v, 5.0);
    if (!approxEqVec(c, v)) {
      std::cerr << "[test_vec3_utils] clampMagnitude(len==max) should be unchanged\n";
      ++fails;
    }
  }

  // clampMagnitude: scales down
  {
    const math::Vec3d v{6.0, 8.0, 0.0}; // len=10
    const math::Vec3d c = math::clampMagnitude(v, 5.0);
    if (!approxEq(c.length(), 5.0, 1e-12) || !approxEqVec(c, {3.0, 4.0, 0.0}, 1e-12)) {
      std::cerr << "[test_vec3_utils] clampMagnitude expected scaled vector (3,4,0)\n";
      ++fails;
    }
  }

  // clampMagnitude: negative max => zero
  {
    const math::Vec3d v{1.0, 2.0, 3.0};
    const math::Vec3d c = math::clampMagnitude(v, -1.0);
    if (!approxEqVec(c, {0.0, 0.0, 0.0})) {
      std::cerr << "[test_vec3_utils] clampMagnitude(max<0) expected zero\n";
      ++fails;
    }
  }

  return fails;
}
