#include "stellar/math/Geometry.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_geometry_raycast() {
  int fails = 0;

  // --- Ray-sphere entry distance ---
  {
    const math::Vec3d o{0, 0, 0};
    const math::Vec3d d{0, 0, 1};
    const math::Vec3d c{0, 0, 10};

    double tEnter = 0.0;
    double tExit = 0.0;
    const bool hit = math::raySphereIntersect(o, d, c, 1.0, tEnter, tExit);
    if (!hit || !approx(tEnter, 9.0, 1e-9) || !approx(tExit, 11.0, 1e-9)) {
      std::cerr << "[test_geometry_raycast] raySphereIntersect wrong. hit=" << hit
                << " tEnter=" << tEnter << " tExit=" << tExit << " expected (9,11)\n";
      ++fails;
    }
  }

  // --- Ray-sphere should reject spheres behind the origin ---
  {
    const math::Vec3d o{0, 0, 0};
    const math::Vec3d d{0, 0, 1};
    const math::Vec3d c{0, 0, -10};
    double tEnter = 0.0;
    double tExit = 0.0;
    const bool hit = math::raySphereIntersect(o, d, c, 1.0, tEnter, tExit);
    if (hit) {
      std::cerr << "[test_geometry_raycast] raySphereIntersect should not hit behind. tEnter=" << tEnter
                << " tExit=" << tExit << "\n";
      ++fails;
    }
  }

  // --- Ray-sphere should treat dir as a direction (not necessarily normalized) ---
  {
    const math::Vec3d o{0, 0, 0};
    const math::Vec3d d{0, 0, 2};
    const math::Vec3d c{0, 0, 10};

    double tEnter = 0.0;
    double tExit = 0.0;
    const bool hit = math::raySphereIntersect(o, d, c, 1.0, tEnter, tExit);
    if (!hit || !approx(tEnter, 9.0, 1e-9) || !approx(tExit, 11.0, 1e-9)) {
      std::cerr << "[test_geometry_raycast] raySphereIntersect should be scale-invariant. hit=" << hit
                << " tEnter=" << tEnter << " tExit=" << tExit << " expected (9,11)\n";
      ++fails;
    }
  }

  // --- Ray-sphere origin inside should clamp entry to 0 ---
  {
    const math::Vec3d o{0, 0, 0};
    const math::Vec3d d{0, 0, 1};
    const math::Vec3d c{0, 0, 0};

    double tEnter = 123.0;
    double tExit = 0.0;
    const bool hit = math::raySphereIntersect(o, d, c, 1.0, tEnter, tExit);
    if (!hit || !approx(tEnter, 0.0, 1e-9) || !(tExit > 0.99 && tExit < 1.01)) {
      std::cerr << "[test_geometry_raycast] raySphereIntersect inside clamp wrong. hit=" << hit
                << " tEnter=" << tEnter << " tExit=" << tExit << " expected (0,1)\n";
      ++fails;
    }
  }

  // --- Segment-sphere hit test ---
  {
    const math::Vec3d a{0, 0, 0};
    const math::Vec3d b{0, 0, 10};
    const math::Vec3d c{0, 0, 5};

    if (!math::segmentHitsSphere(a, b, c, 0.5)) {
      std::cerr << "[test_geometry_raycast] segmentHitsSphere expected hit.\n";
      ++fails;
    }

    if (math::segmentHitsSphere(a, b, {10, 0, 5}, 0.5)) {
      std::cerr << "[test_geometry_raycast] segmentHitsSphere expected miss.\n";
      ++fails;
    }
  }

  // --- Segment closest t ---
  {
    const math::Vec3d a{0, 0, 0};
    const math::Vec3d b{0, 0, 10};

    const double t = math::segmentClosestT(a, b, {0, 0, 7});
    if (!approx(t, 0.7, 1e-12)) {
      std::cerr << "[test_geometry_raycast] segmentClosestT wrong. t=" << t << " expected 0.7\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_geometry_raycast] PASS\n";
  }
  return fails;
}
