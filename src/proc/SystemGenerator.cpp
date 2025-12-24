#include "stellar/proc/SystemGenerator.h"

#include "stellar/core/Random.h"
#include "stellar/econ/Economy.h"
#include "stellar/math/Math.h"
#include "stellar/proc/NameGenerator.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {

static sim::Star makeStar(sim::StarClass cls, core::SplitMix64& rng) {
  sim::Star s{};
  s.cls = cls;

  auto rr = [&](double a, double b) { return rng.range(a, b); };

  switch (cls) {
    case sim::StarClass::O:
      s.massSol = rr(16.0, 60.0);
      s.radiusSol = rr(6.0, 15.0);
      s.luminositySol = rr(30000.0, 500000.0);
      s.temperatureK = rr(30000.0, 50000.0);
      break;
    case sim::StarClass::B:
      s.massSol = rr(2.1, 16.0);
      s.radiusSol = rr(2.0, 6.0);
      s.luminositySol = rr(25.0, 30000.0);
      s.temperatureK = rr(10000.0, 30000.0);
      break;
    case sim::StarClass::A:
      s.massSol = rr(1.4, 2.1);
      s.radiusSol = rr(1.4, 2.5);
      s.luminositySol = rr(5.0, 25.0);
      s.temperatureK = rr(7500.0, 10000.0);
      break;
    case sim::StarClass::F:
      s.massSol = rr(1.04, 1.4);
      s.radiusSol = rr(1.15, 1.6);
      s.luminositySol = rr(1.5, 5.0);
      s.temperatureK = rr(6000.0, 7500.0);
      break;
    case sim::StarClass::G:
      s.massSol = rr(0.8, 1.04);
      s.radiusSol = rr(0.9, 1.2);
      s.luminositySol = rr(0.6, 1.6);
      s.temperatureK = rr(5200.0, 6000.0);
      break;
    case sim::StarClass::K:
      s.massSol = rr(0.45, 0.8);
      s.radiusSol = rr(0.7, 0.95);
      s.luminositySol = rr(0.08, 0.6);
      s.temperatureK = rr(3700.0, 5200.0);
      break;
    case sim::StarClass::M:
    default:
      s.massSol = rr(0.08, 0.45);
      s.radiusSol = rr(0.1, 0.7);
      s.luminositySol = rr(0.0001, 0.08);
      s.temperatureK = rr(2400.0, 3700.0);
      break;
  }
  return s;
}

static sim::PlanetType pickPlanetType(double aAU, core::SplitMix64& rng) {
  // crude zone heuristic + randomness
  const double r = rng.nextDouble();
  if (aAU < 0.6) {
    return (r < 0.7) ? sim::PlanetType::Rocky : sim::PlanetType::Desert;
  }
  if (aAU < 2.0) {
    if (r < 0.5) return sim::PlanetType::Rocky;
    if (r < 0.8) return sim::PlanetType::Ocean;
    return sim::PlanetType::Ice;
  }
  if (aAU < 6.0) {
    return (r < 0.3) ? sim::PlanetType::GasGiant : sim::PlanetType::Ice;
  }
  return (r < 0.6) ? sim::PlanetType::GasGiant : sim::PlanetType::Ice;
}

static void setPlanetMassRadius(sim::Planet& p, core::SplitMix64& rng) {
  auto rr = [&](double a, double b) { return rng.range(a, b); };
  switch (p.type) {
    case sim::PlanetType::GasGiant:
      p.radiusEarth = rr(3.0, 11.0);
      p.massEarth = rr(20.0, 320.0);
      break;
    case sim::PlanetType::Ocean:
      p.radiusEarth = rr(0.8, 2.2);
      p.massEarth = rr(0.6, 8.0);
      break;
    case sim::PlanetType::Ice:
      p.radiusEarth = rr(0.4, 1.6);
      p.massEarth = rr(0.2, 5.0);
      break;
    case sim::PlanetType::Desert:
      p.radiusEarth = rr(0.6, 2.0);
      p.massEarth = rr(0.4, 7.0);
      break;
    case sim::PlanetType::Rocky:
    default:
      p.radiusEarth = rr(0.3, 1.8);
      p.massEarth = rr(0.1, 6.0);
      break;
  }
}

static econ::StationType pickStationType(core::SplitMix64& rng, double bias) {
  // bias -1..+1 (agri <-> industrial)
  const double r = rng.nextDouble();
  if (r < 0.15) return econ::StationType::Outpost;
  if (r < 0.35) return econ::StationType::Mining;
  if (r < 0.55) return econ::StationType::Refinery;
  if (r < 0.75) {
    return (bias < 0.0) ? econ::StationType::Agricultural : econ::StationType::Industrial;
  }
  if (r < 0.88) return econ::StationType::TradeHub;
  if (r < 0.96) return econ::StationType::Research;
  return econ::StationType::Shipyard;
}

static const sim::Faction* findFaction(core::u32 id, const std::vector<sim::Faction>& factions) {
  for (const auto& f : factions) if (f.id == id) return &f;
  return nullptr;
}

sim::StarSystem generateSystem(const sim::SystemStub& stub, const std::vector<sim::Faction>& factions) {
  core::SplitMix64 rng(stub.seed);

  sim::StarSystem sys{};
  sys.stub = stub;

  sys.star = makeStar(stub.primaryClass, rng);

  NameGenerator ng(stub.seed);

  // Planets
  double a = rng.range(0.25, 0.6); // start AU
  const int nPlanets = std::max(0, stub.planetCount);

  sys.planets.reserve(static_cast<std::size_t>(nPlanets));

  for (int i = 0; i < nPlanets; ++i) {
    sim::Planet p{};
    p.name = ng.planetName(stub.name, i);
    a *= rng.range(1.35, 1.9);
    a += rng.range(0.05, 0.25);

    p.orbit.semiMajorAxisAU = a;
    p.orbit.eccentricity = rng.range(0.0, 0.18);
    p.orbit.inclinationRad = rng.range(0.0, stellar::math::degToRad(6.0));
    p.orbit.ascendingNodeRad = rng.range(0.0, 2.0*stellar::math::kPi);
    p.orbit.argPeriapsisRad = rng.range(0.0, 2.0*stellar::math::kPi);
    p.orbit.meanAnomalyAtEpochRad = rng.range(0.0, 2.0*stellar::math::kPi);
    p.orbit.epochDays = 0.0;

    // Kepler-ish: P(years)^2 = a(AU)^3 / M(star)
    const double years = std::sqrt((a*a*a) / std::max(0.08, sys.star.massSol));
    p.orbit.periodDays = years * 365.25;

    p.type = pickPlanetType(a, rng);
    setPlanetMassRadius(p, rng);

    sys.planets.push_back(std::move(p));
  }

  // Stations
  const int nStations = std::max(0, stub.stationCount);
  sys.stations.reserve(static_cast<std::size_t>(nStations));

  const sim::Faction* fac = findFaction(stub.factionId, factions);
  const double fee = fac ? fac->taxRate : 0.02;
  const double bias = fac ? fac->industryBias : 0.0;

  for (int i = 0; i < nStations; ++i) {
    sim::Station st{};
    st.id = core::hashCombine(static_cast<core::u64>(stub.id), static_cast<core::u64>(i + 1));
    st.name = ng.stationName(stub.name, i);
    st.factionId = stub.factionId;
    st.feeRate = fee;
    st.type = pickStationType(rng, bias);
    st.economyModel = econ::makeEconomyModel(st.type, bias);

    // ----------------------
    // Physical placement
    // ----------------------
    // For now stations are placed on simple heliocentric orbits, biased toward
    // existing planet semi-major axes so they feel "in-system".
    double aAU = 0.0;
    if (!sys.planets.empty()) {
      const auto& p = sys.planets[(std::size_t)i % sys.planets.size()];
      const double jitter = rng.range(0.01, 0.08) * (rng.nextDouble() < 0.5 ? -1.0 : 1.0);
      aAU = std::max(0.15, p.orbit.semiMajorAxisAU + jitter);
    } else {
      aAU = rng.range(0.35, 4.0);
    }

    st.orbit.semiMajorAxisAU = aAU;
    st.orbit.eccentricity = rng.range(0.0, 0.03);
    st.orbit.inclinationRad = rng.range(0.0, stellar::math::degToRad(2.5));
    st.orbit.ascendingNodeRad = rng.range(0.0, 2.0 * stellar::math::kPi);
    st.orbit.argPeriapsisRad = rng.range(0.0, 2.0 * stellar::math::kPi);
    st.orbit.meanAnomalyAtEpochRad = rng.range(0.0, 2.0 * stellar::math::kPi);
    st.orbit.epochDays = 0.0;
    {
      const double years = std::sqrt((aAU * aAU * aAU) / std::max(0.08, sys.star.massSol));
      st.orbit.periodDays = years * 365.25;
    }

    // ----------------------
    // Docking parameters
    // ----------------------
    // Keep this lightweight and tunable. Values are chosen so that docking is
    // achievable with early-game thruster settings.
    double radiusKm = 12.0;
    double speedLimitKmS = 0.10;
    double corridorLenKm = 3000.0;
    double halfAngleRad = stellar::math::degToRad(14.0);

    switch (st.type) {
      case econ::StationType::Outpost:
        radiusKm = 8.0;
        speedLimitKmS = 0.15;
        corridorLenKm = 2000.0;
        halfAngleRad = stellar::math::degToRad(16.0);
        break;
      case econ::StationType::Agricultural:
      case econ::StationType::Mining:
      case econ::StationType::Refinery:
      case econ::StationType::Research:
        radiusKm = 12.0;
        speedLimitKmS = 0.12;
        corridorLenKm = 2800.0;
        halfAngleRad = stellar::math::degToRad(14.0);
        break;
      case econ::StationType::Industrial:
      case econ::StationType::TradeHub:
        radiusKm = 18.0;
        speedLimitKmS = 0.10;
        corridorLenKm = 3500.0;
        halfAngleRad = stellar::math::degToRad(13.0);
        break;
      case econ::StationType::Shipyard:
        radiusKm = 26.0;
        speedLimitKmS = 0.08;
        corridorLenKm = 4200.0;
        halfAngleRad = stellar::math::degToRad(12.0);
        break;
      default:
        break;
    }

    st.radiusKm = radiusKm;
    st.dockingSpeedLimitKmS = speedLimitKmS;
    st.dockingCorridorLengthKm = corridorLenKm;
    st.dockingCorridorHalfAngleRad = halfAngleRad;
    st.dockingRangeKm = std::max(60.0, radiusKm * 4.0);

    sys.stations.push_back(std::move(st));
  }

  return sys;
}

} // namespace stellar::proc
