#include "stellar/sim/LambertSolver.h"

#include "stellar/sim/Gravity.h"

#include "stellar/math/Math.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) { return std::abs(a - b) <= eps; }

static bool nearVec(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps) {
  return near(a.x, b.x, eps) && near(a.y, b.y, eps) && near(a.z, b.z, eps);
}

int test_lambert_solver_multirev() {
  int fails = 0;

  const double mu = stellar::sim::muFromEarthMass(1.0);

  const double rKm = 7000.0;
  const double vKmS = std::sqrt(mu / rKm);
  const double period = 2.0 * stellar::math::kPi * std::sqrt((rKm * rKm * rKm) / mu);

  const stellar::math::Vec3d r1{rKm, 0, 0};
  const stellar::math::Vec3d r2{0, rKm, 0};

  // Quarter-orbit geometry, but delayed by one full revolution.
  // A circular-orbit solution should exist with M=1.
  const double dt = period * 1.25;

  stellar::sim::LambertOptions opt;
  opt.longWay = false;
  opt.prograde = true;
  opt.refNormal = {0, 0, 1};
  opt.maxIterations = 160;
  opt.tolSec = 1e-4;

  const auto multi = stellar::sim::solveLambertUniversalMultiRev(r1, r2, dt, mu, 1, opt);
  if (!multi.ok || multi.solutions.empty()) {
    std::cerr << "[test_lambert_solver_multirev] solver returned no solutions\n";
    return ++fails;
  }

  bool foundM1 = false;
  for (const auto& sol : multi.solutions) {
    if (!sol.ok) continue;
    if (sol.revolutions != 1) continue;

    foundM1 = true;

    const stellar::math::Vec3d v1Exp{0, vKmS, 0};
    const stellar::math::Vec3d v2Exp{-vKmS, 0, 0};

    const double eps = 2.5e-3; // 2.5 m/s
    if (!nearVec(sol.v1KmS, v1Exp, eps) || !nearVec(sol.v2KmS, v2Exp, eps)) {
      std::cerr << "[test_lambert_solver_multirev] M=1: v mismatch\n"
                << "  got v1=(" << sol.v1KmS.x << "," << sol.v1KmS.y << "," << sol.v1KmS.z << ")\n"
                << "  exp v1=(" << v1Exp.x << "," << v1Exp.y << "," << v1Exp.z << ")\n"
                << "  got v2=(" << sol.v2KmS.x << "," << sol.v2KmS.y << "," << sol.v2KmS.z << ")\n"
                << "  exp v2=(" << v2Exp.x << "," << v2Exp.y << "," << v2Exp.z << ")\n";
      ++fails;
    }

    // sanity: sqrt(z) should be near (2π + π/2) for a circular M=1 quarter-phase
    const double s = std::sqrt(std::max(0.0, sol.z));
    const double sExp = 2.5 * stellar::math::kPi;
    if (!near(s, sExp, 1e-2)) {
      std::cerr << "[test_lambert_solver_multirev] M=1: expected sqrt(z) ~ " << sExp << ", got " << s << "\n";
      ++fails;
    }
  }

  if (!foundM1) {
    std::cerr << "[test_lambert_solver_multirev] expected an M=1 solution but none found\n";
    ++fails;
  }

  if (fails == 0) std::cout << "[test_lambert_solver_multirev] pass\n";
  return fails;
}
