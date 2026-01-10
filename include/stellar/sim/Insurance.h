#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/ShipLoadout.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Insurance / rebuy (core gameplay)
// -----------------------------------------------------------------------------
//
// The prototype originally applied a flat penalty on ship destruction.
// This module formalizes a replacement-cost model so:
//   - death has meaningful risk that scales with ship progression
//   - the economy loop (shipyard upgrades) ties into combat/survival
//   - the logic is deterministic + unit testable

struct InsurancePolicy {
  // Rebuy is computed as max(minRebuyCr, shipValueCr * rebuyRate).
  double rebuyRate{0.05};
  double minRebuyCr{200.0};

  // Maximum amount of cumulative debt the insurer will extend.
  // If the player cannot cover the rebuy (credits + remaining loan headroom),
  // they are issued a basic loaner ship.
  double loanMaxCr{5000.0};

  // When bankrupt, how many credits (if any) the player keeps.
  // Keeping this at 0.0 is a clean "fresh start".
  double bankruptcyKeepCreditsCr{0.0};
};

// Snapshot of the player's ship progression relevant to replacement cost.
// (Kept separate from SaveGame so it can be used by apps and tests.)
struct PlayerShipEconomyState {
  ShipHullClass hull{ShipHullClass::Scout};
  int thrusterMk{1};
  int shieldMk{1};
  int distributorMk{1};
  WeaponType weaponPrimary{WeaponType::BeamLaser};
  WeaponType weaponSecondary{WeaponType::Cannon};

  // Upgrades purchased at shipyards.
  double cargoCapacityKg{420.0};
  double fuelMax{45.0};
  int passengerSeats{2};
  double fsdRangeLy{18.0};
  int smuggleHoldMk{0};
  int fuelScoopMk{0};
};

struct ShipValueBreakdown {
  double hullCr{0.0};
  double coreModulesCr{0.0};
  double weaponsCr{0.0};
  double upgradesCr{0.0};
  double totalCr{0.0};
};

ShipValueBreakdown computeShipValue(const PlayerShipEconomyState& s);

double computeRebuyCost(const InsurancePolicy& policy, double shipValueCr);

enum class RebuyResult : core::u8 {
  Paid = 0,
  Loan = 1,
  BankruptReset = 2,
};

struct RebuyOutcome {
  RebuyResult result{RebuyResult::Paid};
  double shipValueCr{0.0};
  double rebuyCostCr{0.0};

  double paidFromCreditsCr{0.0};
  double loanTakenCr{0.0};
  double debtAfterCr{0.0};

  bool shipReset{false};
};

// Applies insurance rebuy to credits/debt. May reset ship progression to a basic
// loaner if the player cannot cover the rebuy even with the loan headroom.
RebuyOutcome applyRebuyOnDeath(const InsurancePolicy& policy,
                              PlayerShipEconomyState& ship,
                              double& credits,
                              double& debtCr);

} // namespace stellar::sim
