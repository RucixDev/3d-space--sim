#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Celestial.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Traffic Escort Contracts (deterministic, headless)
// -----------------------------------------------------------------------------
//
// The prototype spawns "traffic convoys" as physical ships that represent a
// background trade-lane shipment. This module provides a deterministic way to
// plan a lightweight escort contract around such a convoy (time window + reward).
//
// Like the other encounter planning helpers (Distress, Derelict), this module
// intentionally does *not* spawn ships or manage runtime state; it only produces
// a stable plan that the caller can interpret.

struct TrafficEscortPlan {
  // When false, the caller should not offer a contract.
  bool offer{true};

  // 0..1-ish risk score used for tuning/reward shaping.
  double risk01{0.0};

  // Player must remain within this range of the convoy during the escort.
  double maxRangeKm{160000.0};

  // Duration the escort must be maintained (sim time, in days).
  double durationDays{0.0};

  // Base payout for completing the escort.
  double rewardCr{0.0};

  // Optional bonus per pirate destroyed while defending the convoy.
  double bonusPerPirateCr{0.0};

  // Reputation awarded with the payer faction on success.
  double repReward{0.0};
};

// Plan an escort contract for a traffic convoy.
//
// Inputs:
//  - universeSeed: global seed
//  - systemId: current system (locality)
//  - convoyId: unique id for the convoy instance
//  - timeDays: current time (only integer day mixed in, stabilizes within a day)
//  - timeToArriveDays: scheduled remaining time until the convoy reaches its destination
//  - cargoValueCr: approximate cargo value of the shipment
//  - piracy01/security01/contest01: system-level security knobs (0..1)
//  - piratesPresent: if pirates are already on-grid attacking the convoy
TrafficEscortPlan planTrafficEscortContract(core::u64 universeSeed,
                                            SystemId systemId,
                                            core::u64 convoyId,
                                            double timeDays,
                                            double timeToArriveDays,
                                            double cargoValueCr,
                                            double piracy01,
                                            double security01,
                                            double contest01,
                                            bool piratesPresent = false);

} // namespace stellar::sim
