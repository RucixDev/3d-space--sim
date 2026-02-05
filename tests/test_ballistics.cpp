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

  
  // --- Accelerating target: verify |r + v t + 0.5 a t^2| ~= s t ---
  {
    const math::Vec3d sPos{0, 0, 0};
    const math::Vec3d sVel{0, 0, 0};
    const math::Vec3d tPos{1500, 0, 0};
    const math::Vec3d tVel{0, 15, 0};
    const math::Vec3d tAcc{0, 0.02, 0}; // 0.02 km/s^2 = 20 m/s^2
    const double speed = 220.0;

    const auto sol = sim::solveProjectileLeadAccel(sPos, sVel, tPos, tVel, tAcc, speed, /*maxTimeSec*/200.0);
    if (!sol) {
      std::cerr << "[test_ballistics] expected accel lead solution.\n";
      ++fails;
    } else {
      const double tSec = sol->tSec;
      const math::Vec3d r = (tPos - sPos) + (tVel - sVel) * tSec + tAcc * (0.5 * tSec * tSec);
      const double lhs = r.length();
      const double rhs = speed * tSec;
      if (!approx(lhs, rhs, 2e-4)) {
        std::cerr << "[test_ballistics] accel intercept mismatch: lhs=" << lhs << " rhs=" << rhs << "\n";
        ++fails;
      }

      const math::Vec3d expectedLead = tPos + tVel * tSec + tAcc * (0.5 * tSec * tSec);
      if ((sol->leadPointKm - expectedLead).length() > 1e-6) {
        std::cerr << "[test_ballistics] accel lead point mismatch.\n";
        ++fails;
      }
    }
  }

  // --- Accel = 0: accel solver matches constant-velocity solver ---
  {
    const math::Vec3d sPos{0, 0, 0};
    const math::Vec3d sVel{0, 0, 0};
    const math::Vec3d tPos{1200, 300, 0};
    const math::Vec3d tVel{-5, 8, 0};
    const math::Vec3d tAcc{0, 0, 0};
    const double speed = 180.0;

    const auto sol0 = sim::solveProjectileLead(sPos, sVel, tPos, tVel, speed, /*maxTimeSec*/200.0);
    const auto sol1 = sim::solveProjectileLeadAccel(sPos, sVel, tPos, tVel, tAcc, speed, /*maxTimeSec*/200.0);
    if (!sol0 || !sol1) {
      std::cerr << "[test_ballistics] expected both CV and accel(0) solutions.\n";
      ++fails;
    } else {
      if (!approx(sol0->tSec, sol1->tSec, 1e-6)) {
        std::cerr << "[test_ballistics] accel(0) time mismatch.\n";
        ++fails;
      }
      if ((sol0->leadPointKm - sol1->leadPointKm).length() > 1e-6) {
        std::cerr << "[test_ballistics] accel(0) lead point mismatch.\n";
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
