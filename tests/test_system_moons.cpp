#include "stellar/proc/SystemGenerator.h"

#include "stellar/sim/Faction.h"
#include "stellar/sim/Units.h"

#include "test_harness.h"

#include <cmath>
#include <iostream>
#include <limits>

using namespace stellar;

namespace {

constexpr double kEarthMassSol = 3.0034146856628466e-6; // (M_earth / M_sun)
constexpr double kEarthRadiusAU = 4.258750455597227e-5; // (R_earth / AU)

double hillRadiusAU(const sim::Planet& p, const sim::Star& s) {
  const double mStar = std::max(0.08, s.massSol);
  const double mPlanet = std::max(0.0, p.massEarth * kEarthMassSol);
  if (mPlanet <= 0.0) {
    return 0.0;
  }
  return p.orbit.semiMajorAxisAU * std::cbrt(mPlanet / (3.0 * mStar));
}

double planetRadiusAU(const sim::Planet& p) {
  return std::max(0.0, p.radiusEarth) * kEarthRadiusAU;
}

bool approxEq(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

sim::SystemStub makeTestStub(core::u64 seed) {
  sim::SystemStub stub{};
  stub.id = seed;
  stub.seed = seed;
  stub.name = "Test";
  stub.posLy = {0.0, 0.0, 0.0};
  stub.primaryClass = sim::StarClass::G;
  stub.planetCount = 12;
  stub.stationCount = 2;
  stub.factionId = 0;
  return stub;
}

} // namespace

int test_system_moons() {
  int failures = 0;
  const core::u64 universeSeed = 0xC0FFEEULL;
  const auto factions = sim::generateFactions(universeSeed, 4);

  // Determinism sanity: same stub -> identical moons.
  {
    const sim::SystemStub stub = makeTestStub(1337ULL);
    const sim::StarSystem a = proc::generateSystem(universeSeed, stub, factions);
    const sim::StarSystem b = proc::generateSystem(universeSeed, stub, factions);

    CHECK(a.moons.size() == b.moons.size());
    if (a.moons.size() != b.moons.size()) {
      return 1;
    }

    for (std::size_t i = 0; i < a.moons.size(); ++i) {
      const auto& ma = a.moons[i];
      const auto& mb = b.moons[i];

      CHECK(ma.id == mb.id);
      CHECK(ma.parentPlanetIndex == mb.parentPlanetIndex);
      CHECK(ma.type == mb.type);

      CHECK(approxEq(ma.radiusEarth, mb.radiusEarth));
      CHECK(approxEq(ma.massEarth, mb.massEarth));

      CHECK(approxEq(ma.orbit.semiMajorAxisAU, mb.orbit.semiMajorAxisAU));
      CHECK(approxEq(ma.orbit.eccentricity, mb.orbit.eccentricity));
      CHECK(approxEq(ma.orbit.inclinationRad, mb.orbit.inclinationRad));
      CHECK(approxEq(ma.orbit.ascendingNodeRad, mb.orbit.ascendingNodeRad));
      CHECK(approxEq(ma.orbit.argPeriapsisRad, mb.orbit.argPeriapsisRad));
      CHECK(approxEq(ma.orbit.meanAnomalyAtEpochRad, mb.orbit.meanAnomalyAtEpochRad));
      CHECK(approxEq(ma.orbit.epochDays, mb.orbit.epochDays));
      CHECK(approxEq(ma.orbit.periodDays, mb.orbit.periodDays));
    }
  }

  // Invariants: find a seed with at least one moon, then validate stability constraints.
  sim::StarSystem sys{};
  bool found = false;
  for (core::u64 k = 1; k <= 512; ++k) {
    auto stub = makeTestStub(0xA11CE000ULL + k * 9973ULL);
    auto candidate = proc::generateSystem(universeSeed, stub, factions);
    if (!candidate.moons.empty()) {
      sys = std::move(candidate);
      found = true;
      break;
    }
  }

  if (!found) {
    CHECK(false && "couldn't find any moons in scan range");
    return 1;
  }

  for (const auto& m : sys.moons) {
    CHECK(m.parentPlanetIndex < sys.planets.size());
    if (m.parentPlanetIndex >= sys.planets.size()) {
      continue;
    }

    const auto& host = sys.planets[m.parentPlanetIndex];
    const double hillAU = hillRadiusAU(host, sys.star);
    const double aAU = m.orbit.semiMajorAxisAU;
    const double rMinAU = planetRadiusAU(host) * 4.0;

    CHECK(aAU > rMinAU);

    if (hillAU > 0.0) {
      CHECK(aAU < hillAU * 0.50);
    }

    const auto posKm = sim::moonPosKm(host, m, 0.0);
    CHECK(std::isfinite(posKm.x));
    CHECK(std::isfinite(posKm.y));
    CHECK(std::isfinite(posKm.z));
  }

  return failures ? 1 : 0;
}
