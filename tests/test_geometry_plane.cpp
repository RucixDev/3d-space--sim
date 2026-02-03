#include "stellar/math/Geometry.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_geometry_plane() {
  int fails = 0;

  const math::Vec3d planePoint{0, 0, 0};
  const math::Vec3d nUnit{0, 1, 0};

  // --- Basic crossing ---
  {
    const math::Vec3d a{0, -1, 0};
    const math::Vec3d b{0, 1, 0};

    double t = -1.0;
    math::Vec3d p{};
    const bool hit = math::segmentPlaneIntersection(a, b, planePoint, nUnit, t, p);
    if (!hit || !approx(t, 0.5, 1e-12) || !approx(p.y, 0.0, 1e-12)) {
      std::cerr << "[test_geometry_plane] crossing failed. hit=" << hit
                << " t=" << t << " p=(" << p.x << "," << p.y << "," << p.z << ")\n";
      ++fails;
    }
  }

  // --- Non-unit normal should behave identically ---
  {
    const math::Vec3d a{0, -2, 0};
    const math::Vec3d b{0, 2, 0};

    double t = -1.0;
    const bool hit = math::segmentPlaneIntersectionT(a, b, planePoint, {0, 2, 0}, t);
    if (!hit || !approx(t, 0.5, 1e-12)) {
      std::cerr << "[test_geometry_plane] non-unit normal failed. hit=" << hit
                << " t=" << t << " expected 0.5\n";
      ++fails;
    }
  }

  // --- Parallel segment not on the plane should miss ---
  {
    const math::Vec3d a{0, 1, 0};
    const math::Vec3d b{1, 1, 0};

    double t = 0.0;
    const bool hit = math::segmentPlaneIntersectionT(a, b, planePoint, nUnit, t);
    if (hit) {
      std::cerr << "[test_geometry_plane] parallel miss failed (unexpected hit). t=" << t << "\n";
      ++fails;
    }
  }

  // --- Segment starting on the plane is treated as a hit at t=0 ---
  {
    const math::Vec3d a{0, 0, 0};
    const math::Vec3d b{0, 1, 0};

    double t = -1.0;
    const bool hit = math::segmentPlaneIntersectionT(a, b, planePoint, nUnit, t);
    if (!hit || !approx(t, 0.0, 1e-12)) {
      std::cerr << "[test_geometry_plane] start-on-plane failed. hit=" << hit
                << " t=" << t << " expected 0.0\n";
      ++fails;
    }
  }

  // --- Segment entirely on the plane is also treated as a hit (t=0) ---
  {
    const math::Vec3d a{-2, 0, 0};
    const math::Vec3d b{2, 0, 0};

    double t = -1.0;
    const bool hit = math::segmentPlaneIntersectionT(a, b, planePoint, nUnit, t);
    if (!hit || !approx(t, 0.0, 1e-12)) {
      std::cerr << "[test_geometry_plane] on-plane segment failed. hit=" << hit
                << " t=" << t << " expected 0.0\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_geometry_plane] PASS\n";
  }

  return fails;
}
