#include "stellar/sim/LambertSolver.h"

#include "stellar/sim/Gravity.h"

#include "stellar/math/Math.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) { return std::abs(a - b) <= eps; }

static bool nearVec(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps) {
  return near(a.x, b.x, eps) && near(a.y, b.y, eps) && near(a.z, b.z, eps);
}

int test_lambert_solver() {
  int fails = 0;

  // Earth-like mu (derived from physical constants used by sim::Gravity).
  const double mu = stellar::sim::muFromEarthMass(1.0);

  // Simple circular orbit geometry.
  const double rKm = 7000.0;
  const double vKmS = std::sqrt(mu / rKm);
  const double period = 2.0 * stellar::math::kPi * std::sqrt((rKm * rKm * rKm) / mu);

  const stellar::math::Vec3d r1{rKm, 0, 0};
  const stellar::math::Vec3d r2{0, rKm, 0};

  // --- Short-way (90 deg), quarter period ---
  {
    const double dt = period * 0.25;

    stellar::sim::LambertOptions opt;
    opt.longWay = false;
    opt.prograde = true;
    opt.refNormal = {0, 0, 1};
    opt.maxIterations = 80;
    opt.tolSec = 1e-4;

    const auto sol = stellar::sim::solveLambertUniversal(r1, r2, dt, mu, opt);
    if (!sol.ok) {
      std::cerr << "[test_lambert_solver] short-way: solver failed\n";
      ++fails;
    } else {
      const stellar::math::Vec3d v1Exp{0, vKmS, 0};
      const stellar::math::Vec3d v2Exp{-vKmS, 0, 0};

      // The bisection-based solver should be very close for this canonical case.
      const double eps = 2.0e-3; // 2 m/s
      if (!nearVec(sol.v1KmS, v1Exp, eps) || !nearVec(sol.v2KmS, v2Exp, eps)) {
        std::cerr << "[test_lambert_solver] short-way: v mismatch\n"
                  << "  got v1=(" << sol.v1KmS.x << "," << sol.v1KmS.y << "," << sol.v1KmS.z << ")\n"
                  << "  exp v1=(" << v1Exp.x << "," << v1Exp.y << "," << v1Exp.z << ")\n"
                  << "  got v2=(" << sol.v2KmS.x << "," << sol.v2KmS.y << "," << sol.v2KmS.z << ")\n"
                  << "  exp v2=(" << v2Exp.x << "," << v2Exp.y << "," << v2Exp.z << ")\n";
        ++fails;
      }

      if (!near(sol.transferAngleRad, stellar::math::kPi * 0.5, 1e-3)) {
        std::cerr << "[test_lambert_solver] short-way: expected transfer angle ~90 deg, got "
                  << stellar::math::radToDeg(sol.transferAngleRad) << " deg\n";
        ++fails;
      }
    }
  }

  // --- Long-way (270 deg), three-quarter period ---
  {
    const double dt = period * 0.75;

    stellar::sim::LambertOptions opt;
    opt.longWay = true;
    opt.prograde = true;
    opt.refNormal = {0, 0, 1};
    opt.maxIterations = 120;
    opt.tolSec = 1e-4;

    const auto sol = stellar::sim::solveLambertUniversal(r1, r2, dt, mu, opt);
    if (!sol.ok) {
      std::cerr << "[test_lambert_solver] long-way: solver failed\n";
      ++fails;
    } else {
      // For this timing, the solution matches a retrograde circular orbit.
      const stellar::math::Vec3d v1Exp{0, -vKmS, 0};
      const stellar::math::Vec3d v2Exp{vKmS, 0, 0};

      const double eps = 2.0e-3; // 2 m/s
      if (!nearVec(sol.v1KmS, v1Exp, eps) || !nearVec(sol.v2KmS, v2Exp, eps)) {
        std::cerr << "[test_lambert_solver] long-way: v mismatch\n"
                  << "  got v1=(" << sol.v1KmS.x << "," << sol.v1KmS.y << "," << sol.v1KmS.z << ")\n"
                  << "  exp v1=(" << v1Exp.x << "," << v1Exp.y << "," << v1Exp.z << ")\n"
                  << "  got v2=(" << sol.v2KmS.x << "," << sol.v2KmS.y << "," << sol.v2KmS.z << ")\n"
                  << "  exp v2=(" << v2Exp.x << "," << v2Exp.y << "," << v2Exp.z << ")\n";
        ++fails;
      }

      if (!near(sol.transferAngleRad, stellar::math::kPi * 1.5, 1e-3)) {
        std::cerr << "[test_lambert_solver] long-way: expected transfer angle ~270 deg, got "
                  << stellar::math::radToDeg(sol.transferAngleRad) << " deg\n";
        ++fails;
      }
    }
  }

  if (fails == 0) std::cout << "[test_lambert_solver] pass\n";
  return fails;
}
