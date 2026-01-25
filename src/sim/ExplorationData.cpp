#include "stellar/sim/ExplorationData.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Scan keys
// -----------------------------------------------------------------------------

// We use stable salted hashCombine chains rather than packing bits so we can
// freely change the ID space of underlying objects without worrying about
// overflow, and so we can keep each key type collision-resistant.
//
// IMPORTANT: Salts are versioned. Changing a salt will invalidate existing
// SaveGame::scannedKeys, so only do that when bumping save versions.
static constexpr core::u64 kSaltStar = 0x3e8d2d1b2c8d4a11ULL;   // arbitrary constant
static constexpr core::u64 kSaltPlanet = 0x7a18b9e5d1c2f3a4ULL;
static constexpr core::u64 kSaltStation = 0x9c5c61f29a7e4d12ULL;
static constexpr core::u64 kSaltSignal = 0x6f0c2e3a9d1b7c55ULL;
static constexpr core::u64 kSaltAsteroid = 0x4b19d6aa31e0c113ULL;
static constexpr core::u64 kSaltSystemComplete = 0x12c8e4a9f0d31b77ULL;

// A helper to mix a human-readable tag into the salt in a deterministic way.
// We still keep explicit constants above so the IDs are fully stable even if
// the hash implementation changes.
static core::u64 tagSalt(core::u64 base, const char* tag) {
  return core::hashCombine(base, core::fnv1a64(tag));
}

core::u64 scanKeyStar(SystemId systemId) {
  return core::hashCombine(tagSalt(kSaltStar, "scan_star_v1"), systemId);
}

core::u64 scanKeyPlanet(SystemId systemId, std::size_t planetIndex) {
  core::u64 h = tagSalt(kSaltPlanet, "scan_planet_v1");
  h = core::hashCombine(h, systemId);
  h = core::hashCombine(h, static_cast<core::u64>(planetIndex));
  return h;
}

core::u64 scanKeyStation(StationId stationId) {
  return core::hashCombine(tagSalt(kSaltStation, "scan_station_v1"), stationId);
}

core::u64 scanKeySignal(core::u64 signalId) {
  return core::hashCombine(tagSalt(kSaltSignal, "scan_signal_v1"), signalId);
}

core::u64 scanKeyAsteroid(core::u64 asteroidId) {
  return core::hashCombine(tagSalt(kSaltAsteroid, "scan_asteroid_v1"), asteroidId);
}

core::u64 scanKeySystemComplete(SystemId systemId) {
  return core::hashCombine(tagSalt(kSaltSystemComplete, "scan_system_complete_v1"), systemId);
}

// -----------------------------------------------------------------------------
// Value formulas
// -----------------------------------------------------------------------------

static double clampFinite(double v, double lo, double hi) {
  if (!std::isfinite(v)) return lo;
  return std::clamp(v, lo, hi);
}

static double starClassRarityMul(StarClass cls) {
  // Rough rarity ladder: hotter, more massive stars are rarer and thus worth more.
  // These are intentionally tame so exploration doesn't dominate the economy.
  switch (cls) {
    case StarClass::O: return 8.0;
    case StarClass::B: return 5.5;
    case StarClass::A: return 3.8;
    case StarClass::F: return 2.4;
    case StarClass::G: return 1.8;
    case StarClass::K: return 1.4;
    case StarClass::M: return 1.0;
    default: return 1.0;
  }
}

double scanValueStar(const Star& star) {
  const double base = 45.0;
  const double rarity = starClassRarityMul(star.cls);

  const double tempK = clampFinite(star.temperatureK, 2000.0, 50000.0);
  const double lum = clampFinite(star.luminositySol, 0.01, 500.0);
  const double mass = clampFinite(star.massSol, 0.05, 80.0);

  // Temperature contributes gently (log scale) and luminosity/mass add a mild
  // bonus. We intentionally avoid huge swings to keep payouts readable.
  const double tempBonus = 0.75 + 0.25 * std::log10(tempK / 2500.0 + 1.0);
  const double lumBonus = 0.85 + 0.15 * std::log10(lum + 1.0);
  const double massBonus = 0.85 + 0.15 * std::log10(mass + 1.0);

  double v = base * rarity * tempBonus * lumBonus * massBonus;
  v = clampFinite(v, 10.0, 1200.0);
  return v;
}

static double planetTypeBase(PlanetType t) {
  switch (t) {
    case PlanetType::Ocean: return 220.0;
    case PlanetType::Ice: return 150.0;
    case PlanetType::Desert: return 135.0;
    case PlanetType::GasGiant: return 170.0;
    case PlanetType::Rocky: return 105.0;
    default: return 105.0;
  }
}

double scanValuePlanet(const Planet& planet) {
  const double base = planetTypeBase(planet.type);
  const double r = clampFinite(planet.radiusEarth, 0.15, 25.0);
  const double m = clampFinite(planet.massEarth, 0.01, 400.0);

  // Size drives value more than mass (players intuitively see radius differences).
  // Use a sublinear curve so gas giants don't explode payouts.
  const double sizeMul = std::pow(r, 0.42);
  const double massMul = 0.85 + 0.15 * std::log10(m + 1.0);

  double v = base * sizeMul * massMul;

  // Small extra for potentially habitable-ish categories.
  if (planet.type == PlanetType::Ocean) v *= 1.15;
  if (planet.type == PlanetType::Desert) v *= 1.05;

  v = clampFinite(v, 15.0, 2500.0);
  return v;
}

static double stationTypeBase(econ::StationType t) {
  switch (t) {
    case econ::StationType::Research: return 140.0;
    case econ::StationType::Shipyard: return 130.0;
    case econ::StationType::TradeHub: return 115.0;
    case econ::StationType::Industrial: return 105.0;
    case econ::StationType::Refinery: return 95.0;
    case econ::StationType::Mining: return 90.0;
    case econ::StationType::Agricultural: return 85.0;
    case econ::StationType::Outpost: return 75.0;
    default: return 75.0;
  }
}

double scanValueStation(const Station& station) {
  const double base = stationTypeBase(station.type);
  // Large stations are a bit more valuable to catalogue.
  const double r = clampFinite(station.radiusKm, 500.0, 30000.0);
  const double sizeMul = 0.85 + 0.15 * std::log10(r / 500.0 + 1.0);
  double v = base * sizeMul;
  v = clampFinite(v, 10.0, 600.0);
  return v;
}

double scanValueSignal(SignalKind kind) {
  double v = 0.0;
  switch (kind) {
    case SignalKind::ResourceField: v = 65.0; break;
    case SignalKind::Derelict: v = 120.0; break;
    case SignalKind::Distress: v = 105.0; break;
    case SignalKind::MissionSalvage: v = 75.0; break;
    case SignalKind::TrafficConvoy: v = 55.0; break;
    default: v = 60.0; break;
  }
  return clampFinite(v, 5.0, 500.0);
}

double scanValueAsteroidProspect() {
  // Prospecting is fast and repeatable across many asteroids, so keep it small.
  return 35.0;
}

double scanValueSystemSurveyBonus(int planetCount, int stationCount) {
  const int pc = std::max(0, planetCount);
  const int sc = std::max(0, stationCount);

  // Completion bonus scales gently with total objects.
  double v = 120.0 + 45.0 * std::sqrt((double)pc) + 35.0 * std::sqrt((double)sc);

  // Tiny nudge for "busy" systems.
  v += 8.0 * std::log1p((double)pc + 0.5 * (double)sc);

  return clampFinite(v, 40.0, 2000.0);
}

// -----------------------------------------------------------------------------
// Broker multipliers
// -----------------------------------------------------------------------------

static double distanceLy(const SystemStub& a, const SystemStub& b) {
  const math::Vec3d d = a.posLy - b.posLy;
  const double dist = d.length();
  if (!std::isfinite(dist)) return 0.0;
  return std::max(0.0, dist);
}

double explorationDataStationDemandMultiplier(econ::StationType stationType,
                                             LogbookEntryKind kind,
                                             core::u8 /*subKind*/) {
  // Base multipliers per station archetype.
  //
  // Philosophy:
  //  - Research stations value scientific scans (stars/planets/derelicts/distress)
  //  - Mining/refinery value resource fields + prospecting
  //  - Trade hubs value traffic signals + station scans
  //  - Shipyards value station scans + system survey completion (nav charts)
  //  - Agricultural value habitable-ish planet types slightly more
  //
  // Sub-kinds are currently unused here (kept for future fine-tuning).

  double m = 1.0;
  switch (stationType) {
    case econ::StationType::Research: {
      if (kind == LogbookEntryKind::StarScan || kind == LogbookEntryKind::PlanetScan) m = 1.20;
      else if (kind == LogbookEntryKind::SignalScan) m = 1.12;
      else if (kind == LogbookEntryKind::SystemSurveyBonus) m = 1.15;
      else m = 0.98;
    } break;
    case econ::StationType::Mining:
    case econ::StationType::Refinery: {
      if (kind == LogbookEntryKind::AsteroidProspect) m = 1.22;
      else if (kind == LogbookEntryKind::SignalScan) m = 1.10;
      else if (kind == LogbookEntryKind::PlanetScan) m = 0.95;
      else m = 0.98;
    } break;
    case econ::StationType::TradeHub: {
      if (kind == LogbookEntryKind::StationScan) m = 1.12;
      else if (kind == LogbookEntryKind::SignalScan) m = 1.08;
      else if (kind == LogbookEntryKind::SystemSurveyBonus) m = 1.10;
      else m = 1.00;
    } break;
    case econ::StationType::Shipyard: {
      if (kind == LogbookEntryKind::StationScan) m = 1.20;
      else if (kind == LogbookEntryKind::SystemSurveyBonus) m = 1.18;
      else if (kind == LogbookEntryKind::StarScan) m = 1.05;
      else m = 0.98;
    } break;
    case econ::StationType::Industrial: {
      if (kind == LogbookEntryKind::StationScan) m = 1.08;
      else if (kind == LogbookEntryKind::SignalScan) m = 1.05;
      else m = 1.00;
    } break;
    case econ::StationType::Agricultural: {
      if (kind == LogbookEntryKind::PlanetScan) m = 1.08;
      else if (kind == LogbookEntryKind::SystemSurveyBonus) m = 1.05;
      else m = 0.98;
    } break;
    case econ::StationType::Outpost:
    default:
      m = 1.0;
      break;
  }

  if (!std::isfinite(m) || m <= 0.0) m = 1.0;
  return m;
}

double explorationDataBrokerMultiplier(const ExplorationDataBrokerParams& params,
                                       const SystemStub& saleSystem,
                                       const Station& saleStation,
                                       const SystemStub& scanSystem,
                                       LogbookEntryKind kind,
                                       core::u8 subKind) {
  double m = 1.0;

  // Distance premium.
  if (params.enableDistancePremium) {
    const double dist = distanceLy(saleSystem, scanSystem);
    const double scale = std::max(1e-6, params.distanceScaleLy);
    const double t = std::clamp(dist / scale, 0.0, 1.0);
    m *= (1.0 + std::clamp(params.maxDistancePremium, 0.0, 2.0) * t);
  }

  // Faction alignment. We treat scanSystem.factionId as the "local" authority.
  if (saleStation.factionId != 0 && scanSystem.factionId != 0) {
    if (saleStation.factionId == scanSystem.factionId) {
      m *= (1.0 + std::clamp(params.sameFactionBonus, -1.0, 2.0));
    } else {
      m *= (1.0 - std::clamp(params.otherFactionPenalty, 0.0, 1.0));
    }
  }

  // Station demand shaping.
  if (params.enableStationDemand && params.demandStrength > 1e-6) {
    const double raw = explorationDataStationDemandMultiplier(saleStation.type, kind, subKind);
    const double s = std::clamp(params.demandStrength, 0.0, 1.0);
    // Blend toward the demand multiplier.
    const double blended = (1.0 - s) * 1.0 + s * raw;
    m *= blended;
  }

  m = clampFinite(m, params.minMultiplier, params.maxMultiplier);
  return m;
}

} // namespace stellar::sim
