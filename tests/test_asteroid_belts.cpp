#include "test_harness.h"

#include "stellar/proc/AsteroidBeltGenerator.h"

#include "stellar/core/Hash.h"
#include "stellar/math/Math.h"

#include <cmath>
#include <cstdint>
#include <vector>

using namespace stellar;

static inline long long q(double v, double s) {
  return (long long)std::llround(v * s);
}

static core::u64 hashBeltPoints(const std::vector<proc::AsteroidBeltPoint>& pts) {
  core::u64 h = core::hashCombine(core::fnv1a64("belt_points"), (core::u64)pts.size());
  for (const auto& p : pts) {
    // Quantize in AU so the hash is robust to tiny floating differences.
    h = core::hashCombine(h, (core::u64)q(p.posAU.x, 16384.0));
    h = core::hashCombine(h, (core::u64)q(p.posAU.y, 16384.0));
    h = core::hashCombine(h, (core::u64)q(p.posAU.z, 16384.0));
    h = core::hashCombine(h, (core::u64)q(p.density01, 65535.0));
  }
  return h;
}

static double meanDensity(const proc::AsteroidBelt& b, double aAU) {
  const int n = 16;
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    const double th = (double)i / (double)n * (2.0 * math::kPi);
    sum += proc::asteroidBeltDensity01(b, aAU, th);
  }
  return sum / (double)n;
}

static sim::Planet makePlanet(const char* name, sim::PlanetType type, double aAU, double mEarth, double incDeg) {
  sim::Planet p;
  p.name = name;
  p.type = type;
  p.radiusEarth = (type == sim::PlanetType::GasGiant) ? 11.2 : 1.0;
  p.massEarth = mEarth;
  p.orbit.semiMajorAxisAU = aAU;
  p.orbit.eccentricity = 0.02;
  p.orbit.inclinationRad = math::degToRad(incDeg);
  p.orbit.ascendingNodeRad = math::degToRad(15.0);
  p.orbit.argPeriapsisRad = 0.0;
  p.orbit.meanAnomalyAtEpochRad = math::degToRad(30.0);
  p.orbit.epochDays = 0.0;
  p.orbit.periodDays = 365.25 * std::sqrt(aAU * aAU * aAU);
  return p;
}

int test_asteroid_belts() {
  int failures = 0;

  sim::StarSystem sys;
  sys.stub.id = 99;
  sys.stub.seed = 42;
  sys.stub.name = "TestSys";
  sys.star.cls = sim::StarClass::G;
  sys.star.massSol = 1.0;
  sys.star.luminositySol = 1.0;
  sys.star.radiusSol = 1.0;
  sys.star.temperatureK = 5778.0;

  // A simple Solar-ish layout: inner terrestrials + a Jupiter-like giant.
  sys.planets.push_back(makePlanet("A", sim::PlanetType::Rocky,   0.40,   0.06, 0.5));
  sys.planets.push_back(makePlanet("B", sim::PlanetType::Ocean,   1.00,   1.00, 0.8));
  sys.planets.push_back(makePlanet("C", sim::PlanetType::Desert,  1.60,   0.11, 1.2));
  sys.planets.push_back(makePlanet("D", sim::PlanetType::GasGiant, 5.20, 317.0, 1.0));
  sys.planets.push_back(makePlanet("E", sim::PlanetType::Ice,     10.0,  14.0, 1.6));

  const core::u64 universeSeed = 1337;

  const auto p0 = proc::generateAsteroidBelts(universeSeed, sys);
  const auto p1 = proc::generateAsteroidBelts(universeSeed, sys);

  CHECK(p0.belts.size() == p1.belts.size());
  CHECK(!p0.belts.empty());

  // Find the main belt and ensure it has at least one resonance feature.
  const proc::AsteroidBelt* mainBelt = nullptr;
  for (const auto& b : p0.belts) {
    if (b.kind == proc::AsteroidBeltKind::MainBelt) {
      mainBelt = &b;
      break;
    }
  }
  CHECK(mainBelt != nullptr);
  CHECK(mainBelt->aOuterAU > mainBelt->aInnerAU);
  CHECK(mainBelt->controllingPlanetIndex >= 0);
  CHECK(!mainBelt->resonances.empty());

  // The 3:1 resonance should carve a visible dip in mean density.
  const proc::BeltResonanceFeature* r31 = nullptr;
  for (const auto& r : mainBelt->resonances) {
    if (r.m == 3 && r.n == 1 && r.strength01 > 0.0) {
      r31 = &r;
      break;
    }
  }
  CHECK(r31 != nullptr);

  const double dAt = meanDensity(*mainBelt, r31->aAU);
  const double dOff = meanDensity(*mainBelt, r31->aAU + r31->halfWidthAU * 3.0);
  CHECK(dAt + 0.03 < dOff);

  // Deterministic point sampling.
  const auto pts0 = proc::sampleAsteroidBeltPoints(universeSeed, sys, *mainBelt, 512, 18);
  const auto pts1 = proc::sampleAsteroidBeltPoints(universeSeed, sys, *mainBelt, 512, 18);
  CHECK(pts0.size() == pts1.size());
  CHECK(hashBeltPoints(pts0) == hashBeltPoints(pts1));

  // Changing universe seed should change the sampled belt.
  const auto pts2 = proc::sampleAsteroidBeltPoints(universeSeed + 1u, sys, *mainBelt, 512, 18);
  CHECK(hashBeltPoints(pts2) != hashBeltPoints(pts0));

  return failures;
}
