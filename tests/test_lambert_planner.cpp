#include "stellar/sim/LambertPlanner.h"

#include "stellar/sim/Gravity.h"
#include "stellar/math/Math.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) { return std::abs(a - b) <= eps; }

int test_lambert_planner() {
  int fails = 0;

  // Earth-like mu (derived from physical constants used by sim::Gravity).
  const double mu = stellar::sim::muFromEarthMass(1.0);

  // Canonical circular orbit case (same as test_lambert_solver, but exercised via planner).
  const double rKm = 7000.0;
  const double vKmS = std::sqrt(mu / rKm);
  const double period = 2.0 * stellar::math::kPi * std::sqrt((rKm * rKm * rKm) / mu);

  const stellar::math::Vec3d r1{rKm, 0, 0};
  const stellar::math::Vec3d v1{0, vKmS, 0};
  const stellar::math::Vec3d r2{0, rKm, 0};
  const stellar::math::Vec3d v2{-vKmS, 0, 0};

  stellar::sim::EphemerisFn dep = [&](double /*timeDays*/, stellar::math::Vec3d& posKm, stellar::math::Vec3d& velKmS) {
    posKm = r1;
    velKmS = v1;
  };
  stellar::sim::EphemerisFn arr = [&](double /*timeDays*/, stellar::math::Vec3d& posKm, stellar::math::Vec3d& velKmS) {
    posKm = r2;
    velKmS = v2;
  };

  stellar::sim::LambertPorkchopParams p;
  p.departMinSec = 0.0;
  p.departMaxSec = 0.0;
  p.departSteps = 1;

  // Choose a grid where 0.25*period is an exact sample:
  //  tof = 0.05P + k*0.05P, k=0..18 => includes k=4 -> 0.25P.
  p.tofMinSec = period * 0.05;
  p.tofMaxSec = period * 0.95;
  p.tofSteps = 19;

  p.scoreMode = stellar::sim::LambertScoreMode::MinTotalDv;
  p.topK = 3;
  p.storeGrid = true;
  p.lambertOpt.longWay = false;
  p.lambertOpt.prograde = true;
  p.lambertOpt.refNormal = {0, 0, 0}; // auto
  p.lambertOpt.maxIterations = 80;
  p.lambertOpt.tolSec = 1e-4;

  const auto res = stellar::sim::searchLambertPorkchop(0.0, dep, arr, mu, p);

  if (res.departSteps != 1 || res.tofSteps != 19) {
    std::cerr << "[test_lambert_planner] unexpected grid dims\n";
    ++fails;
  }
  if (res.grid.size() != 19) {
    std::cerr << "[test_lambert_planner] expected grid size 19, got " << res.grid.size() << "\n";
    ++fails;
  }
  if (res.best.empty()) {
    std::cerr << "[test_lambert_planner] no best candidates returned\n";
    ++fails;
  } else {
    const auto& best = res.best[0];
    if (!best.ok) {
      std::cerr << "[test_lambert_planner] best candidate not ok\n";
      ++fails;
    }

    const double expTof = period * 0.25;
    if (!near(best.tofSec, expTof, 1e-6)) {
      std::cerr << "[test_lambert_planner] expected TOF " << expTof << " got " << best.tofSec << "\n";
      ++fails;
    }

    // The Lambert solver is already tested; here we just sanity check the planner wiring.
    const double tolKmS = 3.0e-3; // 3 m/s
    if (best.dvDepartMagKmS > tolKmS || best.arriveRelSpeedKmS > tolKmS) {
      std::cerr << "[test_lambert_planner] expected near-zero burns, got dvDep=" << best.dvDepartMagKmS
                << " km/s arrRel=" << best.arriveRelSpeedKmS << " km/s\n";
      ++fails;
    }
  }

  // Determinism: same query -> same best result.
  const auto res2 = stellar::sim::searchLambertPorkchop(0.0, dep, arr, mu, p);
  if (!res.best.empty() && !res2.best.empty()) {
    if (!near(res.best[0].tofSec, res2.best[0].tofSec, 1e-12) ||
        !near(res.best[0].score, res2.best[0].score, 1e-12)) {
      std::cerr << "[test_lambert_planner] nondeterministic result\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_lambert_planner] pass\n";
  return fails;
}
