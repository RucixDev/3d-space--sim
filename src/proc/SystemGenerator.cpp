#include "stellar/proc/SystemGenerator.h"

#include "stellar/core/Random.h"
#include "stellar/econ/Economy.h"
#include "stellar/math/Math.h"
#include "stellar/proc/NameGenerator.h"

#include <algorithm>
#include <cmath>

// NOTE:
// This file intentionally stays deterministic from (SystemStub.seed + galaxy factions).
// Be careful when introducing new RNG draws; it will perturb downstream procedural
// content and may require updating signature-based regression tests.

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

static double stationControlWeight(econ::StationType t) {
  // Mirrors the lightweight weighting in sim::computeSystemControl().
  // Keep this local to proc so we don't add a proc->sim dependency.
  using econ::StationType;
  switch (t) {
    case StationType::TradeHub: return 1.35;
    case StationType::Shipyard: return 1.30;
    case StationType::Research: return 1.20;
    case StationType::Industrial: return 1.15;
    case StationType::Refinery: return 1.10;
    case StationType::Mining: return 1.05;
    case StationType::Agricultural: return 1.05;
    case StationType::Outpost: return 1.0;
    default: return 1.0;
  }
}

static void applyFactionToStation(sim::Station& st,
                                  core::u32 factionId,
                                  const std::vector<sim::Faction>& factions) {
  st.factionId = factionId;
  const sim::Faction* fac = findFaction(factionId, factions);
  const double fee = fac ? fac->taxRate : 0.02;
  const double bias = fac ? fac->industryBias : 0.0;
  st.feeRate = fee;
  st.economyModel = econ::makeEconomyModel(st.type, bias);
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

  // ---------------------------------------------------------------------------
  // Stations
  // ---------------------------------------------------------------------------
  const int nStations = std::max(0, stub.stationCount);
  sys.stations.reserve(static_cast<std::size_t>(nStations));

  const sim::Faction* fac = findFaction(stub.factionId, factions);
  const double fee = fac ? fac->taxRate : 0.02;
  const double bias = fac ? fac->industryBias : 0.0;

  // Helper: physical/docking parameters.
  auto setStationPhysicals = [&](sim::Station& st) {
    auto rr = [&](double x, double y) { return rng.range(x, y); };

    double baseRadius = 6000.0;
    double speed = 0.20;

    switch (st.type) {
      case econ::StationType::Outpost:
        baseRadius = 4500.0;
        speed = 0.18;
        break;
      case econ::StationType::Mining:
        baseRadius = 6500.0;
        speed = 0.20;
        break;
      case econ::StationType::Refinery:
        baseRadius = 7000.0;
        speed = 0.20;
        break;
      case econ::StationType::Agricultural:
        baseRadius = 6500.0;
        speed = 0.22;
        break;
      case econ::StationType::Industrial:
        baseRadius = 8000.0;
        speed = 0.22;
        break;
      case econ::StationType::Research:
        baseRadius = 6000.0;
        speed = 0.20;
        break;
      case econ::StationType::TradeHub:
        baseRadius = 11000.0;
        speed = 0.25;
        break;
      case econ::StationType::Shipyard:
        baseRadius = 13000.0;
        speed = 0.28;
        break;
      default:
        break;
    }

    st.radiusKm = baseRadius * rr(0.85, 1.20);

    // Slot scaled off radius.
    st.slotWidthKm  = st.radiusKm * rr(0.75, 0.95);
    st.slotHeightKm = st.radiusKm * rr(0.30, 0.45);
    st.slotDepthKm  = st.radiusKm * rr(0.90, 1.35);

    // Approach corridor
    st.approachLengthKm = st.radiusKm * rr(8.0, 14.0);
    st.approachRadiusKm = st.radiusKm * rr(1.2, 2.2);
    st.maxApproachSpeedKmS = speed * rr(0.85, 1.15);

    // Comms range for clearance
    st.commsRangeKm = st.radiusKm * rr(10.0, 16.0);
  };

  for (int i = 0; i < nStations; ++i) {
    sim::Station st{};
    st.id = core::hashCombine(static_cast<core::u64>(stub.id), static_cast<core::u64>(i + 1));
    st.name = ng.stationName(stub.name, i);
    st.factionId = stub.factionId;
    st.feeRate = fee;
    st.type = pickStationType(rng, bias);
    st.economyModel = econ::makeEconomyModel(st.type, bias);

    // Place stations on simple orbits around the primary.
    // Prefer to place near an existing planet orbit (feels like inhabited space), else pick a random AU.
    double aAU = rng.range(0.8, 5.5);
    if (!sys.planets.empty()) {
      const int idx = rng.range<int>(0, (int)sys.planets.size() - 1);
      aAU = sys.planets[(std::size_t)idx].orbit.semiMajorAxisAU * rng.range(0.92, 1.08);
    }

    st.orbit.semiMajorAxisAU = std::max(0.35, aAU);
    st.orbit.eccentricity = rng.range(0.0, 0.06);
    st.orbit.inclinationRad = rng.range(0.0, stellar::math::degToRad(3.0));
    st.orbit.ascendingNodeRad = rng.range(0.0, 2.0 * stellar::math::kPi);
    st.orbit.argPeriapsisRad = rng.range(0.0, 2.0 * stellar::math::kPi);
    st.orbit.meanAnomalyAtEpochRad = rng.range(0.0, 2.0 * stellar::math::kPi);
    st.orbit.epochDays = 0.0;

    const double years = std::sqrt((st.orbit.semiMajorAxisAU * st.orbit.semiMajorAxisAU * st.orbit.semiMajorAxisAU) /
                                   std::max(0.08, sys.star.massSol));
    st.orbit.periodDays = years * 365.25;

    setStationPhysicals(st);

    sys.stations.push_back(std::move(st));
  }

  // Border outposts:
  // If a system sits inside the overlapping influence radii of multiple factions
  // (and it has enough stations), spawn a minority-owned outpost for the strongest
  // nearby rival.
  //
  // This makes contested systems naturally appear at faction borders, which then
  // feeds into the SecurityModel (contest01) and any systems that consume it
  // (encounters, mission rewards, etc.).
  //
  // IMPORTANT: This logic does NOT consume RNG; it should only change station
  // ownership/economy parameters, not orbital layouts.
  if (stub.factionId != 0 && nStations >= 2 && sys.stations.size() >= 2) {
    const sim::Faction* ctrlFac = findFaction(stub.factionId, factions);
    if (ctrlFac && ctrlFac->influenceRadiusLy > 1e-6) {
      const double dCtrl = (stub.posLy - ctrlFac->homePosLy).length();
      const double wCtrl = std::max(0.0, 1.0 - dCtrl / ctrlFac->influenceRadiusLy);

      core::u32 bestOtherId = 0;
      double bestOtherW = 0.0;
      double bestOtherD = 1e30;

      for (const auto& f : factions) {
        if (f.id == 0) continue;
        if (f.id == stub.factionId) continue;
        if (f.influenceRadiusLy <= 1e-6) continue;
        const double d = (stub.posLy - f.homePosLy).length();
        if (d >= f.influenceRadiusLy) continue;
        const double w = std::max(0.0, 1.0 - d / f.influenceRadiusLy);
        // Prefer the strongest overlapping influence; tie-break by closer distance then faction id.
        if (w > bestOtherW + 1e-12 ||
            (std::abs(w - bestOtherW) <= 1e-12 && (d < bestOtherD - 1e-9 || (std::abs(d - bestOtherD) <= 1e-9 && f.id < bestOtherId)))) {
          bestOtherW = w;
          bestOtherD = d;
          bestOtherId = f.id;
        }
      }

      if (bestOtherId != 0 && bestOtherW > 0.18) {
        const double border01 = bestOtherW / std::max(1e-9, wCtrl + bestOtherW);

        // Threshold chosen so deep-core systems remain uniform, while overlap regions
        // gain contested presence.
        if (border01 > 0.28) {
          // Keep the "primary" station (Alpha, index 0) under the controlling faction.
          //
          // For 2-station systems we can only do this safely if Alpha is strictly heavier
          // than the outpost we would flip; otherwise the computed system-control could
          // end up ambiguous (tie) or dominated by the minority.
          if (sys.stations.size() == 2) {
            const double w0 = stationControlWeight(sys.stations[0].type);
            const double w1 = stationControlWeight(sys.stations[1].type);
            if (w0 > w1 + 1e-12) {
              applyFactionToStation(sys.stations[1], bestOtherId, factions);
            }
          } else {
            // Choose a "minor" station (lowest weight) among indices [1..].
            std::size_t pickIdx = 1;
            double bestW = 1e30;
            for (std::size_t i = 1; i < sys.stations.size(); ++i) {
              const double sw = stationControlWeight(sys.stations[i].type);
              if (sw < bestW - 1e-12 || (std::abs(sw - bestW) <= 1e-12 && i > pickIdx)) {
                bestW = sw;
                pickIdx = i;
              }
            }

            applyFactionToStation(sys.stations[pickIdx], bestOtherId, factions);
          }
        }
      }
    }
  }

  return sys;
}

} // namespace stellar::proc
