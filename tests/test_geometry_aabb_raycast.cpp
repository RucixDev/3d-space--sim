#include "stellar/math/Geometry.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approxEq(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

int test_geometry_aabb_raycast() {
  int fails = 0;

  // Basic ray hit with normalized direction.
  {
    const math::Aabb3d box = math::Aabb3d::fromMinMax({0.0, 0.0, -6.0}, {2.0, 2.0, -4.0});
    double tEnter = 0.0, tExit = 0.0;
    if (!box.rayIntersectionT({1.0, 1.0, 0.0}, {0.0, 0.0, -1.0}, tEnter, tExit)) {
      std::cerr << "[test_geometry_aabb_raycast] expected ray hit\n";
      ++fails;
    } else {
      if (!approxEq(tEnter, 4.0) || !approxEq(tExit, 6.0)) {
        std::cerr << "[test_geometry_aabb_raycast] expected tEnter=4, tExit=6 (got "
                  << tEnter << ", " << tExit << ")\n";
        ++fails;
      }
    }
  }

  // Ray hit with non-normalized direction should report world distances.
  {
    const math::Aabb3d box = math::Aabb3d::fromMinMax({0.0, 0.0, -6.0}, {2.0, 2.0, -4.0});
    double tEnter = 0.0, tExit = 0.0;
    if (!box.rayIntersectionT({1.0, 1.0, 0.0}, {0.0, 0.0, -2.0}, tEnter, tExit)) {
      std::cerr << "[test_geometry_aabb_raycast] expected non-normalized ray hit\n";
      ++fails;
    } else {
      if (!approxEq(tEnter, 4.0) || !approxEq(tExit, 6.0)) {
        std::cerr << "[test_geometry_aabb_raycast] expected same world distances for scaled dir\n";
        ++fails;
      }
    }
  }

  // Ray starting inside the box clamps entry distance to 0.
  {
    const math::Aabb3d box = math::Aabb3d::fromMinMax({0.0, 0.0, -6.0}, {2.0, 2.0, -4.0});
    double tEnter = 123.0, tExit = 123.0;
    if (!box.rayIntersectionT({1.0, 1.0, -5.0}, {1.0, 0.0, 0.0}, tEnter, tExit)) {
      std::cerr << "[test_geometry_aabb_raycast] expected inside ray hit\n";
      ++fails;
    } else {
      if (!approxEq(tEnter, 0.0) || !approxEq(tExit, 1.0)) {
        std::cerr << "[test_geometry_aabb_raycast] expected inside ray tEnter=0, tExit=1 (got "
                  << tEnter << ", " << tExit << ")\n";
        ++fails;
      }
    }
  }

  // Parallel ray missing the slab should not hit.
  {
    const math::Aabb3d box = math::Aabb3d::fromMinMax({0.0, 0.0, -6.0}, {2.0, 2.0, -4.0});
    double tEnter = 0.0, tExit = 0.0;
    if (box.rayIntersectionT({3.0, 1.0, -5.0}, {0.0, 1.0, 0.0}, tEnter, tExit)) {
      std::cerr << "[test_geometry_aabb_raycast] expected parallel miss\n";
      ++fails;
    }
  }

  // Segment intersection returns normalized entry/exit parameters.
  {
    const math::Aabb3d box = math::Aabb3d::fromMinMax({0.0, 0.0, -6.0}, {2.0, 2.0, -4.0});
    double tEnter = 0.0, tExit = 0.0;
    if (!box.segmentIntersectionT({-1.0, 1.0, -5.0}, {3.0, 1.0, -5.0}, tEnter, tExit)) {
      std::cerr << "[test_geometry_aabb_raycast] expected segment hit\n";
      ++fails;
    } else {
      if (!approxEq(tEnter, 0.25) || !approxEq(tExit, 0.75)) {
        std::cerr << "[test_geometry_aabb_raycast] expected segment tEnter=0.25, tExit=0.75 (got "
                  << tEnter << ", " << tExit << ")\n";
        ++fails;
      }
    }
  }

  // Segment completely outside.
  {
    const math::Aabb3d box = math::Aabb3d::fromMinMax({0.0, 0.0, -6.0}, {2.0, 2.0, -4.0});
    double tEnter = 0.0, tExit = 0.0;
    if (box.segmentIntersectionT({-1.0, 5.0, -5.0}, {3.0, 5.0, -5.0}, tEnter, tExit)) {
      std::cerr << "[test_geometry_aabb_raycast] expected segment miss\n";
      ++fails;
    }
  }

  return fails;
}
