#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/Law.h"

#include <array>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Police scans / contraband helpers (headless)
// -----------------------------------------------------------------------------
//
// The SDL prototype historically implemented contraband scans, bribes, and fines
// directly inside apps/stellar_game/main.cpp.
//
// This module extracts the deterministic *math* (scan rates, bribe chance, fine
// schedule application, and cargo confiscation) into the core library so:
//  - tests can validate behavior deterministically
//  - other apps (stellar_sandbox, future server sims) can reuse the same rules
//  - main.cpp stays focused on orchestration/UI
//
// Design goals:
//  - no renderer/UI dependencies
//  - deterministic given the inputs
//  - stable, bounded outputs

// Helper: scan frequency multiplier from smuggle compartments.
// Matches the prototype tuning in apps/stellar_game.
//
//  mk=0 => 1.00 (no reduction)
//  mk=1 => 0.72
//  mk=2 => 0.50
//  mk=3 => 0.35
// Values outside [0,3] are clamped.
double smuggleHoldScanRateMult(int smuggleHoldMk);

// Helper: scan duration multiplier from smuggle compartments.
//
//  mk=0 => 1.00
//  mk=1 => 1.15
//  mk=2 => 1.35
//  mk=3 => 1.60
// Values outside [0,3] are clamped.
double smuggleHoldScanDurationMult(int smuggleHoldMk);

// Compute the scan start rate (per real second).
//
// NOTE: The caller is expected to convert this rate into a per-frame probability
// using a hazard-rate style transform, e.g.
//   p = 1 - exp(-ratePerSec * dtSeconds)
// so that results are independent of frame rate.
double cargoScanStartRatePerSec(bool hasContraband, const LawProfile& law, int smuggleHoldMk);

// Compute scan durations (seconds) used by the prototype for station/police scans.
double cargoScanDurationSecStation(bool hasContraband, int smuggleHoldMk);
double cargoScanDurationSecPolice(int smuggleHoldMk);

// Result of scanning the player's cargo for contraband.
struct IllegalCargoScanResult {
  double illegalValueCr{0.0};
  std::array<double, econ::kCommodityCount> scannedIllegalUnits{}; // what the scanner "saw"
};

// Scan the given cargo array and compute the total illegal cargo value under a precomputed
// illegality mask.
//
// This is useful when "illegal here" depends on station-specific rules rather than faction-wide rules.
IllegalCargoScanResult scanIllegalCargoMask(core::u32 illegalMask,
                                           const std::array<double, econ::kCommodityCount>& cargoUnits,
                                           const std::array<double, econ::kCommodityCount>* midPriceOverrideCr = nullptr);

// Scan the given cargo array and compute the total illegal cargo value under the
// faction's deterministic contraband rules.
//
// If midPriceOverrideCr is provided, it supplies a per-commodity mid price used
// to value the goods (useful for sizing fines in local market terms). If null,
// commodity base prices are used.
IllegalCargoScanResult scanIllegalCargo(core::u64 universeSeed,
                                       core::u32 factionId,
                                       const std::array<double, econ::kCommodityCount>& cargoUnits,
                                       const std::array<double, econ::kCommodityCount>* midPriceOverrideCr = nullptr);

// Compute the probability that a police scan offers a bribe instead of immediate
// confiscation + fine.
//
// Inputs follow the prototype semantics:
//  - playerRep is the player's reputation with the scanning faction.
//  - playerHeat is the ship heat gauge (0..~100).
double bribeOfferChance(const LawProfile& law, double playerRep, double playerHeat, double illegalValueCr);

// Bribe amount for keeping the cargo (cr). The caller may round for UI.
double bribeAmountCr(const LawProfile& law, double illegalValueCr);

// Convenience: compute a bribe amount and round it to a pleasant UI increment.
// roundToCr <= 0 disables rounding.
double bribeAmountCrRounded(const LawProfile& law, double illegalValueCr, double roundToCr = 10.0);

// Result of enforcing contraband (confiscation + fine).
//
// This function is pure: it does not mutate inputs. The caller applies the
// resulting credits/cargo and any rep/bounty side effects.
struct ContrabandEnforcementResult {
  double creditsAfter{0.0};
  std::array<double, econ::kCommodityCount> cargoAfter{};

  // Confiscated units per commodity (<= scannedIllegalUnits).
  std::array<double, econ::kCommodityCount> confiscatedUnits{};

  // Fine schedule.
  double fineCr{0.0};
  double paidCr{0.0};
  double unpaidCr{0.0};

  // Reputation penalty to apply (negative).
  double repPenalty{0.0};

  // Suggested bounty to add.
  //
  // For contraband compliance, this is intentionally *not* set to the unpaid fine;
  // callers can track unpaidCr as a fine debt (payable later) and only mint a
  // bounty/warrant when a fine becomes overdue or when the player evades.
  double bountyAddedCr{0.0};

  // Suggested local security escalation knobs (prototype tuning).
  double policeHeatDelta{0.0};
  double policeAlertSeconds{0.0};
  double nextPoliceSpawnDelaySeconds{0.0};
};

ContrabandEnforcementResult enforceContraband(const LawProfile& law,
                                             double credits,
                                             const std::array<double, econ::kCommodityCount>& cargoUnits,
                                             const std::array<double, econ::kCommodityCount>& scannedIllegalUnits,
                                             double illegalValueCr);

// Result of evading/refusing a contraband compliance window.
struct ContrabandEvadeResult {
  double repPenalty{0.0};
  double bountyAddedCr{0.0};

  double policeHeatDelta{0.0};
  double policeAlertSeconds{0.0};
  double nextPoliceSpawnDelaySeconds{0.0};

  // Suggested ship heat increase for the prototype (adds urgency).
  double shipHeatDelta{0.0};
};

ContrabandEvadeResult evadeContraband(const LawProfile& law,
                                     double fineCr,
                                     double illegalValueCr,
                                     bool escalateLocal);

} // namespace stellar::sim
