#include "stellar/sim/Ship.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_ship() {
  int fails = 0;

  // --- Boost caps default to legacy multipliers (and track base when not customized). ---
  {
    sim::Ship s;

    s.setMaxLinearAccelKmS2(0.10);
    if (!approx(s.maxLinearAccelBoostKmS2(), 0.18, 1e-12)) {
      std::cerr << "[test_ship] default linear boost cap mismatch after base set. got="
                << s.maxLinearAccelBoostKmS2() << " expected=0.18\n";
      ++fails;
    }

    s.setMaxAngularAccelRadS2(2.0);
    if (!approx(s.maxAngularAccelBoostRadS2(), 2.8, 1e-12)) {
      std::cerr << "[test_ship] default angular boost cap mismatch after base set. got="
                << s.maxAngularAccelBoostRadS2() << " expected=2.8\n";
      ++fails;
    }
  }

  // --- Custom boost caps persist even if base caps change later. ---
  {
    sim::Ship s;
    s.setMaxLinearAccelBoostKmS2(0.50);
    s.setMaxLinearAccelKmS2(0.20);
    if (!approx(s.maxLinearAccelBoostKmS2(), 0.50, 1e-12)) {
      std::cerr << "[test_ship] custom linear boost cap should persist. got="
                << s.maxLinearAccelBoostKmS2() << " expected=0.50\n";
      ++fails;
    }
  }

  // --- Physics step uses the configured boost cap (not a baked multiplier). ---
  {
    sim::ShipInput in{};
    in.thrustLocal = {0, 0, 1};
    in.torqueLocal = {0, 0, 0};
    in.dampers = false;
    in.brake = false;

    // Base
    {
      sim::Ship s;
      s.setMaxLinearAccelKmS2(0.10);
      s.setMaxLinearAccelBoostKmS2(0.40);
      in.boost = false;
      s.step(1.0, in);
      if (!approx(s.velocityKmS().z, 0.10, 1e-6)) {
        std::cerr << "[test_ship] base accel step mismatch. got vz=" << s.velocityKmS().z
                  << " expected=0.10\n";
        ++fails;
      }
    }

    // Boost
    {
      sim::Ship s;
      s.setMaxLinearAccelKmS2(0.10);
      s.setMaxLinearAccelBoostKmS2(0.40);
      in.boost = true;
      s.step(1.0, in);
      if (!approx(s.velocityKmS().z, 0.40, 1e-6)) {
        std::cerr << "[test_ship] boost accel step mismatch. got vz=" << s.velocityKmS().z
                  << " expected=0.40\n";
        ++fails;
      }
    }
  }

  // --- Regression guard: defaults should preserve the legacy 1.8x boost multiplier. ---
  {
    sim::Ship s;
    sim::ShipInput in{};
    in.thrustLocal = {0, 0, 1};
    in.dampers = false;
    in.brake = false;
    in.boost = true;
    s.step(1.0, in);

    // Default base accel is 0.05 km/s^2, legacy boost multiplier is 1.8.
    const double expectedVz = 0.05 * 1.8;
    if (!approx(s.velocityKmS().z, expectedVz, 1e-6)) {
      std::cerr << "[test_ship] default boost behavior mismatch. got vz=" << s.velocityKmS().z
                << " expected=" << expectedVz << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_ship] PASS\n";
  }
  return fails;
}
