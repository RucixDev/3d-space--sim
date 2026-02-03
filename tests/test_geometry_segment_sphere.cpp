#include "stellar/math/Geometry.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_geometry_segment_sphere() {
  int fails = 0;

  // --- Simple hit with known entry/exit params ---
  {
    const math::Vec3d a{0, 0, 0};
    const math::Vec3d b{0, 0, 10};
    const math::Vec3d c{0, 0, 5};

    double tEnter = -1.0;
    double tExit = -1.0;
    const bool hit = math::segmentSphereIntersectionT(a, b, c, 1.0, tEnter, tExit);

    if (!hit || !approx(tEnter, 0.4, 1e-12) || !approx(tExit, 0.6, 1e-12)) {
      std::cerr << "[test_geometry_segment_sphere] basic hit failed. hit=" << hit
                << " tEnter=" << tEnter << " tExit=" << tExit
                << " expected (0.4,0.6)\n";
      ++fails;
    }
  }

  // --- Starting inside should clamp entry to 0 ---
  {
    const math::Vec3d a{0, 0, 5};
    const math::Vec3d b{0, 0, 10};
    const math::Vec3d c{0, 0, 5};

    double tEnter = -1.0;
    double tExit = -1.0;
    const bool hit = math::segmentSphereIntersectionT(a, b, c, 1.0, tEnter, tExit);

    if (!hit || !approx(tEnter, 0.0, 1e-12) || !approx(tExit, 0.2, 1e-12)) {
      std::cerr << "[test_geometry_segment_sphere] start-inside failed. hit=" << hit
                << " tEnter=" << tEnter << " tExit=" << tExit
                << " expected (0.0,0.2)\n";
      ++fails;
    }
  }

  // --- Miss should return false ---
  {
    const math::Vec3d a{0, 0, 0};
    const math::Vec3d b{0, 0, 10};
    const math::Vec3d c{1, 0, 5};

    double tEnter = -1.0;
    double tExit = -1.0;
    const bool hit = math::segmentSphereIntersectionT(a, b, c, 0.5, tEnter, tExit);

    if (hit) {
      std::cerr << "[test_geometry_segment_sphere] miss failed (unexpected hit). tEnter=" << tEnter
                << " tExit=" << tExit << "\n";
      ++fails;
    }
  }

  // --- Degenerate segment (point) ---
  {
    const math::Vec3d p{0, 0, 0};
    const math::Vec3d c{0, 0, 0};

    double tEnter = -1.0;
    double tExit = -1.0;
    const bool hit = math::segmentSphereIntersectionT(p, p, c, 1.0, tEnter, tExit);
    if (!hit || !approx(tEnter, 0.0, 1e-12) || !approx(tExit, 0.0, 1e-12)) {
      std::cerr << "[test_geometry_segment_sphere] degenerate-inside failed. hit=" << hit
                << " tEnter=" << tEnter << " tExit=" << tExit << "\n";
      ++fails;
    }

    tEnter = -1.0;
    tExit = -1.0;
    const bool miss = math::segmentSphereIntersectionT(p, p, {5, 0, 0}, 1.0, tEnter, tExit);
    if (miss) {
      std::cerr << "[test_geometry_segment_sphere] degenerate-outside failed (unexpected hit).\n";
      ++fails;
    }
  }

  // --- Convenience entry-only overload ---
  {
    const math::Vec3d a{0, 0, 0};
    const math::Vec3d b{0, 0, 10};
    const math::Vec3d c{0, 0, 5};

    double tEnter = -1.0;
    const bool hit = math::segmentSphereIntersectionEnterT(a, b, c, 1.0, tEnter);
    if (!hit || !approx(tEnter, 0.4, 1e-12)) {
      std::cerr << "[test_geometry_segment_sphere] entry-only overload failed. hit=" << hit
                << " tEnter=" << tEnter << " expected 0.4\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_geometry_segment_sphere] PASS\n";
  }

  return fails;
}
