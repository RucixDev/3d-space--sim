#pragma once

#include "stellar/core/Random.h"
#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/econ/Economy.h"
#include "stellar/sim/Celestial.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>

namespace stellar::sim {

// Bitmask helpers for simple contraband/legality rules.
// Commodity count is intentionally small in the prototype, so a 32-bit mask is sufficient.
inline core::u32 commodityBit(econ::CommodityId cid) {
  return (core::u32)1u << (core::u32)cid;
}

// Deterministic per-faction illegality mask.
//
// This is used by both the SDL prototype and the headless modules (missions / tooling) so that
// "illegal here" is consistent across UI, tests, and any CLI utilities.
//
// Design: each non-zero faction bans 1 "vice" good and, with some probability,
// bans 1 "controlled tech" good.
inline core::u32 illegalCommodityMask(core::u64 universeSeed, core::u32 factionId) {
  if (factionId == 0) return 0u;

  // Large odd constant to improve bit diffusion. Keep stable for save/behavior continuity.
  core::SplitMix64 r(universeSeed ^ (core::u64)factionId * 0x9E3779B97F4A7C15ull);

  core::u32 mask = 0u;

  // Keep these lists small and "gamey" for now: it creates consistent but varied
  // contraband rules without needing a full law simulation.
  const std::array<econ::CommodityId, 3> vice{
      econ::CommodityId::Luxury,
      econ::CommodityId::Medicine,
      econ::CommodityId::Stimulants,
  };
  const std::array<econ::CommodityId, 3> tech{
      econ::CommodityId::Electronics,
      econ::CommodityId::Machinery,
      econ::CommodityId::Weapons,
  };

  mask |= commodityBit(vice[(std::size_t)(r.nextU32() % (core::u32)vice.size())]);

  if (r.nextUnit() < 0.70) {
    mask |= commodityBit(tech[(std::size_t)(r.nextU32() % (core::u32)tech.size())]);
  }

  return mask;
}

inline bool isIllegalCommodity(core::u64 universeSeed, core::u32 factionId, econ::CommodityId cid) {
  if (factionId == 0) return false;
  return (illegalCommodityMask(universeSeed, factionId) & commodityBit(cid)) != 0u;
}

inline std::string illegalCommodityListString(core::u64 universeSeed, core::u32 factionId) {
  const core::u32 mask = illegalCommodityMask(universeSeed, factionId);
  if (mask == 0u) return "None";

  std::string out;
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    if ((mask & ((core::u32)1u << (core::u32)i)) == 0u) continue;
    if (!out.empty()) out += ", ";
    out += std::string(econ::commodityName((econ::CommodityId)i));
  }
  return out.empty() ? "None" : out;
}

inline bool hasIllegalCargo(core::u64 universeSeed,
                            core::u32 factionId,
                            const std::array<double, econ::kCommodityCount>& cargo) {
  const core::u32 mask = illegalCommodityMask(universeSeed, factionId);
  if (mask == 0u) return false;

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    if ((mask & ((core::u32)1u << (core::u32)i)) != 0u) {
      if (cargo[i] > 1e-6) return true;
    }
  }
  return false;
}

// --- Station-specific extensions --------------------------------------------
//
// Factions define a baseline contraband set, but individual stations may also
// apply additional restrictions (or special enforcement) depending on their
// role in the economy (type). This enables a richer smuggling loop:
//  - A commodity can be legal in one port but illegal in another under the same faction.
//
// IMPORTANT:
//  - The station mask is a *superset* of the faction mask in this prototype.
//    (Stations add bans but do not remove faction bans.)
//  - The station layer is deterministic and version-tagged so behavior remains stable.
inline core::u32 illegalCommodityMaskForStation(core::u64 universeSeed,
                                               core::u32 factionId,
                                               StationId stationId,
                                               econ::StationType stationType) {
  if (factionId == 0) return 0u;

  core::u32 mask = illegalCommodityMask(universeSeed, factionId);

  // Station overlay seed per (universe, faction, station).
  // Keep stable across versions via a text tag.
  core::u64 s = universeSeed;
  s ^= (core::u64)factionId * 0xD1B54A32D192ED03ull;
  s ^= (core::u64)stationId * 0x9E3779B97F4A7C15ull;
  s ^= core::seedFromText("station_contraband_v1");
  core::SplitMix64 r(s);

  // Rough security bias by station type.
  double bias = 0.30;
  switch (stationType) {
    case econ::StationType::Outpost:      bias = 0.20; break;
    case econ::StationType::Agricultural: bias = 0.28; break;
    case econ::StationType::Mining:       bias = 0.26; break;
    case econ::StationType::Refinery:     bias = 0.30; break;
    case econ::StationType::Industrial:   bias = 0.34; break;
    case econ::StationType::Research:     bias = 0.44; break;
    case econ::StationType::TradeHub:     bias = 0.38; break;
    case econ::StationType::Shipyard:     bias = 0.50; break;
    default:                              bias = 0.30; break;
  }

  // Chance to add *one* extra restricted commodity.
  const double extraChance = 0.22 + 0.30 * std::clamp(bias, 0.0, 1.0); // ~[0.22..0.37]
  if (r.nextUnit() >= extraChance) {
    return mask;
  }

  auto addIfNew = [&](econ::CommodityId cid) {
    const core::u32 bit = commodityBit(cid);
    if ((mask & bit) != 0u) return false;
    mask |= bit;
    return true;
  };

  // Candidate sets tuned per station type.
  // NOTE: keep candidates within the existing "vice/tech" set so the loop stays readable.
  switch (stationType) {
    case econ::StationType::Outpost: {
      const std::array<econ::CommodityId, 3> c{
          econ::CommodityId::Stimulants,
          econ::CommodityId::Weapons,
          econ::CommodityId::Luxury,
      };
      for (int tries = 0; tries < 6; ++tries) {
        if (addIfNew(c[(std::size_t)(r.nextU32() % (core::u32)c.size())])) break;
      }
      break;
    }
    case econ::StationType::Agricultural: {
      const std::array<econ::CommodityId, 3> c{
          econ::CommodityId::Stimulants,
          econ::CommodityId::Weapons,
          econ::CommodityId::Luxury,
      };
      for (int tries = 0; tries < 6; ++tries) {
        if (addIfNew(c[(std::size_t)(r.nextU32() % (core::u32)c.size())])) break;
      }
      break;
    }
    case econ::StationType::Mining: {
      const std::array<econ::CommodityId, 3> c{
          econ::CommodityId::Weapons,
          econ::CommodityId::Stimulants,
          econ::CommodityId::Luxury,
      };
      for (int tries = 0; tries < 6; ++tries) {
        if (addIfNew(c[(std::size_t)(r.nextU32() % (core::u32)c.size())])) break;
      }
      break;
    }
    case econ::StationType::Refinery: {
      const std::array<econ::CommodityId, 4> c{
          econ::CommodityId::Weapons,
          econ::CommodityId::Electronics,
          econ::CommodityId::Stimulants,
          econ::CommodityId::Medicine,
      };
      for (int tries = 0; tries < 6; ++tries) {
        if (addIfNew(c[(std::size_t)(r.nextU32() % (core::u32)c.size())])) break;
      }
      break;
    }
    case econ::StationType::Industrial: {
      const std::array<econ::CommodityId, 4> c{
          econ::CommodityId::Weapons,
          econ::CommodityId::Electronics,
          econ::CommodityId::Machinery,
          econ::CommodityId::Stimulants,
      };
      for (int tries = 0; tries < 6; ++tries) {
        if (addIfNew(c[(std::size_t)(r.nextU32() % (core::u32)c.size())])) break;
      }
      break;
    }
    case econ::StationType::Research: {
      const std::array<econ::CommodityId, 4> c{
          econ::CommodityId::Weapons,
          econ::CommodityId::Electronics,
          econ::CommodityId::Stimulants,
          econ::CommodityId::Medicine,
      };
      for (int tries = 0; tries < 6; ++tries) {
        if (addIfNew(c[(std::size_t)(r.nextU32() % (core::u32)c.size())])) break;
      }
      break;
    }
    case econ::StationType::TradeHub: {
      const std::array<econ::CommodityId, 4> c{
          econ::CommodityId::Stimulants,
          econ::CommodityId::Weapons,
          econ::CommodityId::Luxury,
          econ::CommodityId::Electronics,
      };
      for (int tries = 0; tries < 6; ++tries) {
        if (addIfNew(c[(std::size_t)(r.nextU32() % (core::u32)c.size())])) break;
      }
      break;
    }
    case econ::StationType::Shipyard: {
      const std::array<econ::CommodityId, 4> c{
          econ::CommodityId::Weapons,
          econ::CommodityId::Electronics,
          econ::CommodityId::Machinery,
          econ::CommodityId::Stimulants,
      };
      for (int tries = 0; tries < 6; ++tries) {
        if (addIfNew(c[(std::size_t)(r.nextU32() % (core::u32)c.size())])) break;
      }
      break;
    }
    default: {
      // Reasonable fallback.
      const std::array<econ::CommodityId, 3> c{
          econ::CommodityId::Weapons,
          econ::CommodityId::Stimulants,
          econ::CommodityId::Electronics,
      };
      for (int tries = 0; tries < 6; ++tries) {
        if (addIfNew(c[(std::size_t)(r.nextU32() % (core::u32)c.size())])) break;
      }
      break;
    }
  }

  return mask;
}

inline bool isIllegalCommodityAtStation(core::u64 universeSeed,
                                       core::u32 factionId,
                                       StationId stationId,
                                       econ::StationType stationType,
                                       econ::CommodityId cid) {
  if (factionId == 0) return false;
  const core::u32 mask = illegalCommodityMaskForStation(universeSeed, factionId, stationId, stationType);
  return (mask & commodityBit(cid)) != 0u;
}

inline std::string illegalCommodityListStringForStation(core::u64 universeSeed,
                                                       core::u32 factionId,
                                                       StationId stationId,
                                                       econ::StationType stationType) {
  const core::u32 mask = illegalCommodityMaskForStation(universeSeed, factionId, stationId, stationType);
  if (mask == 0u) return "None";

  std::string out;
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    if ((mask & ((core::u32)1u << (core::u32)i)) == 0u) continue;
    if (!out.empty()) out += ", ";
    out += std::string(econ::commodityName((econ::CommodityId)i));
  }
  return out.empty() ? "None" : out;
}

inline bool hasIllegalCargoForStation(core::u64 universeSeed,
                                     core::u32 factionId,
                                     StationId stationId,
                                     econ::StationType stationType,
                                     const std::array<double, econ::kCommodityCount>& cargo) {
  const core::u32 mask = illegalCommodityMaskForStation(universeSeed, factionId, stationId, stationType);
  if (mask == 0u) return false;

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    if ((mask & ((core::u32)1u << (core::u32)i)) != 0u) {
      if (cargo[i] > 1e-6) return true;
    }
  }
  return false;
}

} // namespace stellar::sim
