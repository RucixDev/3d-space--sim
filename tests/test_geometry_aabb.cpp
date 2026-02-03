#include "stellar/math/Geometry.h"

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

int test_geometry_aabb() {
  int fails = 0;

  // fromMinMax: orders per-axis
  {
    const math::Aabb3d b = math::Aabb3d::fromMinMax({1.0, -2.0, 5.0}, {-3.0, 4.0, 0.0});
    if (!approxEqVec(b.min, {-3.0, -2.0, 0.0}) || !approxEqVec(b.max, {1.0, 4.0, 5.0})) {
      std::cerr << "[test_geometry_aabb] fromMinMax expected ordered min/max\n";
      ++fails;
    }
  }

  // fromCenterExtents: abs(extents)
  {
    const math::Aabb3d b = math::Aabb3d::fromCenterExtents({1.0, 2.0, 3.0}, {-2.0, 0.5, -4.0});
    if (!approxEqVec(b.min, {-1.0, 1.5, -1.0}) || !approxEqVec(b.max, {3.0, 2.5, 7.0})) {
      std::cerr << "[test_geometry_aabb] fromCenterExtents expected c +/- abs(e)\n";
      ++fails;
    }
  }

  // contains + clampPoint
  {
    const math::Aabb3d b = math::Aabb3d::fromMinMax({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    if (!b.contains({0.0, 0.5, 1.0}) || b.contains({-1.0, 0.0, 0.0})) {
      std::cerr << "[test_geometry_aabb] contains() basic cases failed\n";
      ++fails;
    }

    const math::Vec3d c = b.clampPoint({2.0, -1.0, 0.25});
    if (!approxEqVec(c, {1.0, 0.0, 0.25})) {
      std::cerr << "[test_geometry_aabb] clampPoint expected per-axis clamp\n";
      ++fails;
    }
  }

  // distanceSqToPoint
  {
    const math::Aabb3d b = math::Aabb3d::fromMinMax({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    if (!approxEq(b.distanceSqToPoint({0.25, 0.25, 0.25}), 0.0)) {
      std::cerr << "[test_geometry_aabb] distanceSqToPoint inside should be 0\n";
      ++fails;
    }
    if (!approxEq(b.distanceSqToPoint({2.0, 0.5, 0.5}), 1.0)) {
      std::cerr << "[test_geometry_aabb] distanceSqToPoint x-outside expected 1\n";
      ++fails;
    }
    if (!approxEq(b.distanceSqToPoint({2.0, 2.0, 2.0}), 3.0)) {
      std::cerr << "[test_geometry_aabb] distanceSqToPoint corner expected 3\n";
      ++fails;
    }

    math::Aabb3d empty;
    if (!std::isinf(empty.distanceSqToPoint({0.0, 0.0, 0.0}))) {
      std::cerr << "[test_geometry_aabb] empty box distance should be inf\n";
      ++fails;
    }
  }

  // intersectsSphere
  {
    const math::Aabb3d b = math::Aabb3d::fromMinMax({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    if (!b.intersectsSphere({0.5, 0.5, 0.5}, 0.1)) {
      std::cerr << "[test_geometry_aabb] sphere inside should intersect\n";
      ++fails;
    }
    if (b.intersectsSphere({2.0, 0.5, 0.5}, 0.999)) {
      std::cerr << "[test_geometry_aabb] sphere just short should not intersect\n";
      ++fails;
    }
    if (!b.intersectsSphere({2.0, 0.5, 0.5}, 1.0)) {
      std::cerr << "[test_geometry_aabb] sphere touching should intersect\n";
      ++fails;
    }
  }

  // expand point + expand box
  {
    math::Aabb3d b;
    b.expand({1.0, 2.0, 3.0});
    b.expand({-1.0, 5.0, 0.0});

    if (!approxEqVec(b.min, {-1.0, 2.0, 0.0}) || !approxEqVec(b.max, {1.0, 5.0, 3.0})) {
      std::cerr << "[test_geometry_aabb] expand(point) expected updated bounds\n";
      ++fails;
    }

    b.expand(math::Aabb3d::fromMinMax({10.0, -1.0, 2.0}, {11.0, 0.0, 4.0}));
    if (!approxEqVec(b.min, {-1.0, -1.0, 0.0}) || !approxEqVec(b.max, {11.0, 5.0, 4.0})) {
      std::cerr << "[test_geometry_aabb] expand(box) expected merged bounds\n";
      ++fails;
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    b.expand({nan, 0.0, 0.0}); // ignored
    if (!b.isFinite()) {
      std::cerr << "[test_geometry_aabb] expand(NaN) should not poison the box\n";
      ++fails;
    }
  }

  return fails;
}
