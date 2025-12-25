#include "stellar/sim/Orbit.h"

#include "stellar/math/Math.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps = 1e-4) { return std::abs(a - b) <= eps; }

static bool near3(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps = 1e-4) {
  return near(a.x, b.x, eps) && near(a.y, b.y, eps) && near(a.z, b.z, eps);
}

int test_orbit() {
  int fails = 0;

  stellar::sim::OrbitElements e{};
  e.semiMajorAxisAU = 1.0;
  e.eccentricity = 0.0;
  e.inclinationRad = 0.0;
  e.ascendingNodeRad = 0.0;
  e.argPeriapsisRad = 0.0;
  e.meanAnomalyAtEpochRad = 0.0;
  e.epochDays = 0.0;
  e.periodDays = 100.0;

  auto p0 = stellar::sim::orbitPosition3DAU(e, 0.0);
  if (!near(p0.x, 1.0) || !near(p0.y, 0.0)) {
    std::cerr << "[test_orbit] expected (1,0) at t=0, got (" << p0.x << "," << p0.y << ")\n";
    ++fails;
  }

  auto p25 = stellar::sim::orbitPosition3DAU(e, 25.0);
  if (!near(p25.x, 0.0, 5e-3) || !near(p25.y, 1.0, 5e-3)) {
    std::cerr << "[test_orbit] expected ~ (0,1) at quarter period, got (" << p25.x << "," << p25.y << ")\n";
    ++fails;
  }

  // Velocity sanity (AU/day)
  const double vExpected = (2.0 * stellar::math::kPi * e.semiMajorAxisAU) / e.periodDays;
  const auto v0 = stellar::sim::orbitVelocity3DAU(e, 0.0);
  if (!near3(v0, {0.0, vExpected, 0.0}, 2e-4)) {
    std::cerr << "[test_orbit] expected velocity ~(0," << vExpected << ",0) at t=0, got ("
              << v0.x << "," << v0.y << "," << v0.z << ")\n";
    ++fails;
  }

  const auto v25 = stellar::sim::orbitVelocity3DAU(e, 25.0);
  if (!near3(v25, {-vExpected, 0.0, 0.0}, 2e-3)) {
    std::cerr << "[test_orbit] expected velocity ~(-" << vExpected << ",0,0) at quarter period, got ("
              << v25.x << "," << v25.y << "," << v25.z << ")\n";
    ++fails;
  }

  // Eccentric orbit sanity: should be within [a(1-e), a(1+e)]
  e.eccentricity = 0.4;
  auto pe = stellar::sim::orbitPosition3DAU(e, 0.0);
  const double r = std::sqrt(pe.x*pe.x + pe.y*pe.y + pe.z*pe.z);
  const double rmin = e.semiMajorAxisAU * (1.0 - e.eccentricity) - 1e-6;
  const double rmax = e.semiMajorAxisAU * (1.0 + e.eccentricity) + 1e-6;
  if (!(r >= rmin && r <= rmax)) {
    std::cerr << "[test_orbit] eccentric radius out of bounds r=" << r << " expected [" << rmin << "," << rmax << "]\n";
    ++fails;
  }

  // Velocity matches numeric derivative (eccentric)
  {
    const double t = 12.345;
    const double dt = 1e-4; // days (~8.6s)
    const auto vA = stellar::sim::orbitVelocity3DAU(e, t);
    const auto pF = stellar::sim::orbitPosition3DAU(e, t + dt);
    const auto pB = stellar::sim::orbitPosition3DAU(e, t - dt);
    const auto vN = (pF - pB) / (2.0 * dt);
    const double err = (vA - vN).length();
    if (err > 2e-3) {
      std::cerr << "[test_orbit] velocity mismatch vs numeric derivative err=" << err << " AU/day\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_orbit] pass\n";
  return fails;
}
