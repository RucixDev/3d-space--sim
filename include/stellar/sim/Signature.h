#pragma once

#include "stellar/core/StableHash.h"
#include "stellar/sim/Faction.h"
#include "stellar/sim/System.h"

#include <cstddef>
#include <vector>

namespace stellar::sim {

// Build stable, portable 64-bit signatures for deterministic/procedural content.
//
// These signatures are intended for:
//  - tooling output (easy regression checking)
//  - unit tests (detect accidental procedural drift)
//
// The functions are header-only and deliberately avoid pulling in heavy dependencies.

inline void signatureOrbit(core::StableHash64& h, const OrbitElements& o) {
  h.addDoubleQ(o.semiMajorAxisAU);
  h.addDoubleQ(o.eccentricity);
  h.addDoubleQ(o.inclinationRad);
  h.addDoubleQ(o.ascendingNodeRad);
  h.addDoubleQ(o.argPeriapsisRad);
  h.addDoubleQ(o.meanAnomalyAtEpochRad);
  h.addDoubleQ(o.epochDays);
  h.addDoubleQ(o.periodDays);
}

inline core::u64 signatureSystemStub(const SystemStub& s) {
  core::StableHash64 h;
  h.addU64(s.id);
  h.addU64(s.seed);
  h.addString(s.name);
  h.addDoubleQ(s.posLy.x);
  h.addDoubleQ(s.posLy.y);
  h.addDoubleQ(s.posLy.z);
  h.addU8(static_cast<core::u8>(s.primaryClass));
  h.addInt(s.planetCount);
  h.addInt(s.stationCount);
  h.addU32(s.factionId);
  return h.value();
}

// Deterministic signature for a list of SystemStubs.
//
// Notes:
//  - This function does NOT sort the input; callers should feed a stable order.
inline core::u64 signatureSystemStubList(const std::vector<SystemStub>& stubs) {
  core::StableHash64 h;
  h.addU64(static_cast<core::u64>(stubs.size()));
  for (const auto& s : stubs) {
    h.addU64(signatureSystemStub(s));
  }
  return h.value();
}

inline core::u64 signatureFaction(const Faction& f) {
  core::StableHash64 h;
  h.addU32(f.id);
  h.addString(f.name);
  h.addDoubleQ(f.homePosLy.x);
  h.addDoubleQ(f.homePosLy.y);
  h.addDoubleQ(f.homePosLy.z);
  h.addDoubleQ(f.color.x);
  h.addDoubleQ(f.color.y);
  h.addDoubleQ(f.color.z);
  h.addDoubleQ(f.influenceRadiusLy);
  h.addDoubleQ(f.taxRate);
  h.addDoubleQ(f.industryBias);
  return h.value();
}

inline core::u64 signatureFactionList(const std::vector<Faction>& factions) {
  core::StableHash64 h;
  h.addU64(static_cast<core::u64>(factions.size()));
  for (const auto& f : factions) {
    h.addU64(signatureFaction(f));
  }
  return h.value();
}

inline void signatureStar(core::StableHash64& h, const Star& s) {
  h.addU8(static_cast<core::u8>(s.cls));
  h.addDoubleQ(s.massSol);
  h.addDoubleQ(s.radiusSol);
  h.addDoubleQ(s.luminositySol);
  h.addDoubleQ(s.temperatureK);
}

inline void signaturePlanet(core::StableHash64& h, const Planet& p) {
  h.addString(p.name);
  h.addU8(static_cast<core::u8>(p.type));
  h.addDoubleQ(p.radiusEarth);
  h.addDoubleQ(p.massEarth);
  signatureOrbit(h, p.orbit);
}

inline void signatureEconomyModel(core::StableHash64& h, const econ::StationEconomyModel& m) {
  h.addU8(static_cast<core::u8>(m.type));
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) h.addDoubleQ(m.productionPerDay[i]);
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) h.addDoubleQ(m.consumptionPerDay[i]);
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) h.addDoubleQ(m.desiredStock[i]);
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) h.addDoubleQ(m.capacity[i]);
  h.addDoubleQ(m.priceVolatility);
  h.addDoubleQ(m.shockVolatility);
}

inline void signatureStation(core::StableHash64& h, const Station& s) {
  h.addU64(s.id);
  h.addString(s.name);
  h.addU8(static_cast<core::u8>(s.type));
  h.addU32(s.factionId);
  h.addDoubleQ(s.feeRate);
  signatureEconomyModel(h, s.economyModel);
  signatureOrbit(h, s.orbit);

  h.addDoubleQ(s.radiusKm);
  h.addDoubleQ(s.commsRangeKm);
  h.addDoubleQ(s.approachLengthKm);
  h.addDoubleQ(s.approachRadiusKm);
  h.addDoubleQ(s.maxApproachSpeedKmS);

  h.addDoubleQ(s.slotWidthKm);
  h.addDoubleQ(s.slotHeightKm);
  h.addDoubleQ(s.slotDepthKm);
}

inline core::u64 signatureStarSystem(const StarSystem& sys) {
  core::StableHash64 h;
  h.addU64(signatureSystemStub(sys.stub));
  signatureStar(h, sys.star);

  h.addU64(static_cast<core::u64>(sys.planets.size()));
  for (const auto& p : sys.planets) signaturePlanet(h, p);

  h.addU64(static_cast<core::u64>(sys.stations.size()));
  for (const auto& s : sys.stations) signatureStation(h, s);
  return h.value();
}

// Deterministic signature for a full StarSystem plus the faction set used to
// interpret it.
inline core::u64 signatureStarSystem(const StarSystem& sys, const std::vector<Faction>& factions) {
  core::StableHash64 h;
  h.addU64(signatureStarSystem(sys));
  h.addU64(signatureFactionList(factions));
  return h.value();
}

} // namespace stellar::sim
