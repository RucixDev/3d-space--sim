#include "stellar/sim/LambertPlanner.h"

#include "stellar/sim/Gravity.h"
#include "stellar/math/Math.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) { return std::abs(a - b) <= eps; }

int test_lambert_planner_stepper() {
  int fails = 0;

  // Earth-like mu (derived from physical constants used by sim::Gravity).
  const double mu = stellar::sim::muFromEarthMass(1.0);

  // Canonical circular orbit case (same as test_lambert_planner).
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

  p.tofMinSec = period * 0.05;
  p.tofMaxSec = period * 0.95;
  p.tofSteps = 19;

  p.scoreMode = stellar::sim::LambertScoreMode::MinTotalDv;
  p.topK = 3;
  p.storeGrid = true;
  p.lambertOpt.longWay = false;
  p.lambertOpt.prograde = true;
  p.lambertOpt.refNormal = {0, 0, 0};
  p.lambertOpt.maxIterations = 80;
  p.lambertOpt.tolSec = 1e-4;

  // Reference: one-shot helper.
  const auto ref = stellar::sim::searchLambertPorkchop(0.0, dep, arr, mu, p);

  // Stepper should reproduce the same best result and grid size regardless of chunking.
  stellar::sim::LambertPorkchopStepper s;
  s.start(0.0, dep, arr, mu, p);

  int guard = 0;
  while (!s.done() && guard < 100000) {
    s.step(1); // worst-case chunking
    ++guard;
  }
  if (!s.done()) {
    std::cerr << "[test_lambert_planner_stepper] stepper did not finish\n";
    ++fails;
  }

  const auto res = s.result();

  if (res.grid.size() != ref.grid.size()) {
    std::cerr << "[test_lambert_planner_stepper] grid size mismatch\n";
    ++fails;
  }
  if (res.best.empty() || ref.best.empty()) {
    std::cerr << "[test_lambert_planner_stepper] missing best candidates\n";
    ++fails;
  } else {
    if (!near(res.best[0].tofSec, ref.best[0].tofSec, 1e-12) ||
        !near(res.best[0].score, ref.best[0].score, 1e-12)) {
      std::cerr << "[test_lambert_planner_stepper] best candidate mismatch\n";
      ++fails;
    }

    // Expected TOF from the canonical scenario.
    const double expTof = period * 0.25;
    if (!near(res.best[0].tofSec, expTof, 1e-6)) {
      std::cerr << "[test_lambert_planner_stepper] expected TOF " << expTof << " got " << res.best[0].tofSec << "\n";
      ++fails;
    }
  }

  // progress01() should be 1 when done.
  if (!near(s.progress01(), 1.0, 1e-12)) {
    std::cerr << "[test_lambert_planner_stepper] progress not 1 when done\n";
    ++fails;
  }

  if (fails == 0) std::cout << "[test_lambert_planner_stepper] pass\n";
  return fails;
}
