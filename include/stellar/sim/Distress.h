#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/Celestial.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Distress encounters
// -----------------------------------------------------------------------------
// A deterministic generator for lightweight "signal source" distress encounters.
//
// This module intentionally keeps the data model small so it can be used by both:
//  - the SDL prototype (stellar_game)
//  - headless tools / unit tests
//
// The generator does *not* spawn ships or manage runtime state. It only produces
// a stable plan (victim request + optional ambush parameters) that the caller can
// interpret.

enum class DistressScenario : core::u8 {
  Supplies = 0,   // Food / Water
  Fuel = 1,       // Fuel
  Medical = 2,    // Medicine
  Mechanical = 3, // Machinery / Metals

  // Pure hostile trap (no legitimate victim request).
  Ambush = 4,
};

inline const char* distressScenarioName(DistressScenario s) {
  switch (s) {
    case DistressScenario::Supplies: return "Supplies";
    case DistressScenario::Fuel: return "Fuel";
    case DistressScenario::Medical: return "Medical";
    case DistressScenario::Mechanical: return "Mechanical";
    case DistressScenario::Ambush: return "Ambush";
    default: return "Distress";
  }
}

struct DistressPlan {
  DistressScenario scenario{DistressScenario::Supplies};

  // When true, the caller should spawn a disabled/stranded ship that can be rescued.
  bool hasVictim{true};

  // The requested cargo commodity and quantity (units) to complete the rescue.
  // Only meaningful when hasVictim == true.
  econ::CommodityId needCommodity{econ::CommodityId::Food};
  double needUnits{0.0};

  // Reward schedule (credits + reputation with the payer faction).
  double rewardCr{0.0};
  double repReward{0.0};
  core::u32 payerFactionId{0}; // 0 means "no faction / independent".

  // Optional hostile response.
  bool ambush{false};
  int pirateCount{0};

  // 0..1-ish risk score used only for tuning/reward shaping.
  double risk{0.0};
};

// Generate a deterministic distress plan.
//
// Inputs:
//  - universeSeed: global seed
//  - systemId: used for locality (different systems feel different)
//  - signalId: unique per signal instance (provides per-signal variation)
//  - timeDays: only the integer day component is mixed in (stabilizes within a day)
//  - localFactionId: who is assumed to pay/award rep (0 allowed)
DistressPlan planDistressEncounter(core::u64 universeSeed,
                                   SystemId systemId,
                                   core::u64 signalId,
                                   double timeDays,
                                   core::u32 localFactionId = 0);

} // namespace stellar::sim
