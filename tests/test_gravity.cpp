#include "stellar/sim/Gravity.h"
#include "stellar/sim/OrbitalMechanics.h"

#include "stellar/math/Math.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps = 1e-6) { return std::abs(a - b) <= eps; }

int test_gravity() {
  int fails = 0;

  // --- Single-body acceleration sanity ---
  {
    const double mu = stellar::sim::muFromEarthMass(1.0);
    const double rKm = 7000.0;
    const auto a = stellar::sim::gravityAccelFromBodyKmS2({0,0,0}, mu, {rKm,0,0}, /*minRadiusKm=*/0.0);
    const double expectedAx = -mu / (rKm * rKm);
    if (!near(a.x, expectedAx, 1e-12) || !near(a.y, 0.0, 1e-12) || !near(a.z, 0.0, 1e-12)) {
      std::cerr << "[test_gravity] accel mismatch: got (" << a.x << "," << a.y << "," << a.z
                << ") expected (" << expectedAx << ",0,0)\n";
      ++fails;
    }
  }

  // --- Circular orbit (e ~ 0, a ~ r) ---
  {
    const double mu = stellar::sim::muFromEarthMass(1.0);
    const double rKm = 7000.0;
    const double vKmS = std::sqrt(mu / rKm);
    const auto orb = stellar::sim::solveTwoBodyOrbit({rKm, 0, 0}, {0, vKmS, 0}, mu);

    if (orb.eccentricity > 1e-6) {
      std::cerr << "[test_gravity] expected near-circular e~0, got e=" << orb.eccentricity << "\n";
      ++fails;
    }
    if (!near(orb.semiMajorAxisKm, rKm, 1e-3)) {
      std::cerr << "[test_gravity] expected a~" << rKm << " km, got a=" << orb.semiMajorAxisKm << "\n";
      ++fails;
    }

    const double expectedPeriod = 2.0 * stellar::math::kPi * std::sqrt((rKm * rKm * rKm) / mu);
    if (!near(orb.periodSec, expectedPeriod, 1e-3)) {
      std::cerr << "[test_gravity] expected period~" << expectedPeriod << " s, got " << orb.periodSec << " s\n";
      ++fails;
    }
  }

  // --- Elliptical orbit (perigee 7000 km, apogee 14000 km) ---
  {
    const double mu = stellar::sim::muFromEarthMass(1.0);
    const double rp = 7000.0;
    const double ra = 14000.0;
    const double a = 0.5 * (rp + ra);
    const double e = (ra - rp) / (ra + rp);

    // Speed at periapsis: v = sqrt(mu*(2/r - 1/a))
    const double vp = std::sqrt(mu * (2.0 / rp - 1.0 / a));
    const auto orb = stellar::sim::solveTwoBodyOrbit({rp, 0, 0}, {0, vp, 0}, mu);

    if (!near(orb.semiMajorAxisKm, a, 1e-2)) {
      std::cerr << "[test_gravity] expected a~" << a << " km, got a=" << orb.semiMajorAxisKm << "\n";
      ++fails;
    }
    if (!near(orb.eccentricity, e, 2e-4)) {
      std::cerr << "[test_gravity] expected e~" << e << ", got e=" << orb.eccentricity << "\n";
      ++fails;
    }
    if (!near(orb.periapsisKm, rp, 1e-2) || !near(orb.apoapsisKm, ra, 1e-2)) {
      std::cerr << "[test_gravity] expected rp/ra (" << rp << "," << ra << ") got (" << orb.periapsisKm
                << "," << orb.apoapsisKm << ")\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_gravity] pass\n";
  return fails;
}
