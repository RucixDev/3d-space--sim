#include "stellar/sim/Ballistics.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-6) {
  return std::abs(a - b) <= eps;
}

int test_ballistics() {
  int fails = 0;

  // --- Stationary target, stationary shooter: time = dist / speed ---
  {
    const math::Vec3d sPos{0, 0, 0};
    const math::Vec3d sVel{0, 0, 0};
    const math::Vec3d tPos{1000, 0, 0};
    const math::Vec3d tVel{0, 0, 0};
    const double speed = 100.0;

    const auto t = sim::solveInterceptTimeSec(sPos, sVel, tPos, tVel, speed);
    if (!t || !approx(*t, 10.0, 1e-9)) {
      std::cerr << "[test_ballistics] expected t=10, got " << (t ? *t : -1.0) << "\n";
      ++fails;
    }

    const auto lead = sim::solveProjectileLead(sPos, sVel, tPos, tVel, speed);
    if (!lead) {
      std::cerr << "[test_ballistics] expected lead solution.\n";
      ++fails;
    } else {
      if (!approx(lead->tSec, 10.0, 1e-9)) {
        std::cerr << "[test_ballistics] lead t mismatch.\n";
        ++fails;
      }
      if ((lead->leadPointKm - tPos).length() > 1e-6) {
        std::cerr << "[test_ballistics] lead point mismatch.\n";
        ++fails;
      }
      if (lead->aimDirWorld.x < 0.99) {
        std::cerr << "[test_ballistics] aim dir not pointing toward +x.\n";
        ++fails;
      }
    }
  }

  // --- Moving target: verify the solution satisfies |r + v t| ~= s t ---
  {
    const math::Vec3d sPos{0, 0, 0};
    const math::Vec3d sVel{0, 0, 0};
    const math::Vec3d tPos{2000, 0, 0};
    const math::Vec3d tVel{0, 20, 0};
    const double speed = 200.0;

    const auto t = sim::solveInterceptTimeSec(sPos, sVel, tPos, tVel, speed);
    if (!t) {
      std::cerr << "[test_ballistics] expected moving-target solution.\n";
      ++fails;
    } else {
      const math::Vec3d r = tPos - sPos;
      const math::Vec3d v = tVel - sVel;
      const double lhs = (r + v * (*t)).length();
      const double rhs = speed * (*t);
      if (!approx(lhs, rhs, 1e-5)) {
        std::cerr << "[test_ballistics] intercept equation mismatch: lhs=" << lhs << " rhs=" << rhs << "\n";
        ++fails;
      }
    }
  }

  // --- Unsolvable case: target outruns projectile directly away. ---
  {
    const math::Vec3d sPos{0, 0, 0};
    const math::Vec3d sVel{0, 0, 0};
    const math::Vec3d tPos{1000, 0, 0};
    const math::Vec3d tVel{120, 0, 0};
    const double speed = 80.0;

    const auto t = sim::solveInterceptTimeSec(sPos, sVel, tPos, tVel, speed);
    if (t) {
      std::cerr << "[test_ballistics] expected no solution, got t=" << *t << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_ballistics] PASS\n";
  }
  return fails;
}
