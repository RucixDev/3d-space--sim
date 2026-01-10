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
// intentionally does *not* spawn ships or manage runtime behaviour; it only
// produces a stable plan that the caller can interpret.
//
// The game does keep a small *serializable* runtime state for an accepted
// escort contract (progress, bonuses, etc.). The structs below are plain data
// containers so SaveGame + UI can persist/display contracts without coupling
// to any rendering/engine types.

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

// -----------------------------------------------------------------------------
// Serializable runtime state (no simulation logic)
// -----------------------------------------------------------------------------

// Active escort contract the player has accepted for a given traffic convoy.
//
// Notes:
//  - This is a lightweight objective: stay within range for a short window.
//  - The contract may be offered by the convoy's faction (or the local faction
//    if the convoy is unaffiliated).
struct TrafficEscortContractState {
  bool active{false};
  core::u64 convoyId{0};
  core::u32 payerFactionId{0};
  StationId toStationId{0};

  double startDays{0.0};
  double untilDays{0.0};
  double maxRangeKm{160000.0};
  double tooFarSec{0.0};

  double rewardCr{0.0};
  double bonusPerPirateCr{0.0};
  double repReward{0.0};

  int piratesKilled{0};
  bool piratesPresentAtStart{false};
};

// Anti-farm record: convoys for which the player has already been paid.
// The game may discard old entries after a retention window.
struct TrafficEscortSettlementState {
  core::u64 convoyId{0};
  double settledDay{0.0};
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
