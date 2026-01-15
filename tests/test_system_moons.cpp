#include "stellar/proc/SystemGenerator.h"

#include "stellar/sim/Faction.h"
#include "stellar/sim/Units.h"

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

int main() {
  const core::u64 universeSeed = 0xC0FFEEULL;
  const auto factions = sim::generateFactions(universeSeed, 4);

  // Determinism sanity: same stub -> identical moons.
  {
    const sim::SystemStub stub = makeTestStub(1337ULL);
    const sim::StarSystem a = proc::generateSystem(universeSeed, stub, factions);
    const sim::StarSystem b = proc::generateSystem(universeSeed, stub, factions);

    if (a.moons.size() != b.moons.size()) {
      std::cerr << "FAIL: moon count differs across identical generation\n";
      return 1;
    }

    for (std::size_t i = 0; i < a.moons.size(); ++i) {
      const auto& ma = a.moons[i];
      const auto& mb = b.moons[i];

      if (ma.id != mb.id || ma.parentPlanetIndex != mb.parentPlanetIndex || ma.type != mb.type) {
        std::cerr << "FAIL: moon identity differs at index " << i << "\n";
        return 1;
      }

      if (!approxEq(ma.radiusEarth, mb.radiusEarth) || !approxEq(ma.massEarth, mb.massEarth)) {
        std::cerr << "FAIL: moon physicals differ at index " << i << "\n";
        return 1;
      }

      if (!approxEq(ma.orbit.semiMajorAxisAU, mb.orbit.semiMajorAxisAU) ||
          !approxEq(ma.orbit.eccentricity, mb.orbit.eccentricity) ||
          !approxEq(ma.orbit.inclinationRad, mb.orbit.inclinationRad) ||
          !approxEq(ma.orbit.ascendingNodeRad, mb.orbit.ascendingNodeRad) ||
          !approxEq(ma.orbit.argPeriapsisRad, mb.orbit.argPeriapsisRad) ||
          !approxEq(ma.orbit.meanAnomalyAtEpochRad, mb.orbit.meanAnomalyAtEpochRad) ||
          !approxEq(ma.orbit.epochDays, mb.orbit.epochDays) ||
          !approxEq(ma.orbit.periodDays, mb.orbit.periodDays)) {
        std::cerr << "FAIL: moon orbit differs at index " << i << "\n";
        return 1;
      }
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
    std::cerr << "FAIL: couldn't find any moons in scan range\n";
    return 1;
  }

  for (const auto& m : sys.moons) {
    if (m.parentPlanetIndex >= sys.planets.size()) {
      std::cerr << "FAIL: moon has invalid parent index\n";
      return 1;
    }

    const auto& host = sys.planets[m.parentPlanetIndex];
    const double hillAU = hillRadiusAU(host, sys.star);
    const double aAU = m.orbit.semiMajorAxisAU;
    const double rMinAU = planetRadiusAU(host) * 4.0;

    if (!(aAU > rMinAU)) {
      std::cerr << "FAIL: moon orbit is implausibly small (inside min safe radius)\n";
      return 1;
    }

    if (hillAU > 0.0 && !(aAU < hillAU * 0.50)) {
      std::cerr << "FAIL: moon orbit exceeds stability fraction of Hill radius\n";
      return 1;
    }

    const auto posKm = sim::moonPosKm(host, m, 0.0);
    if (!std::isfinite(posKm.x) || !std::isfinite(posKm.y) || !std::isfinite(posKm.z)) {
      std::cerr << "FAIL: moonPosKm returned non-finite values\n";
      return 1;
    }
  }

  std::cout << "OK: moons generated and invariants satisfied (moons=" << sys.moons.size() << ")\n";
  return 0;
}
