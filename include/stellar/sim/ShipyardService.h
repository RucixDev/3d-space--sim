#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Economy.h"
#include "stellar/sim/FactionProfile.h"
#include "stellar/sim/ShipLoadout.h"
#include "stellar/sim/System.h"

#include <cstddef>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// ShipyardService (headless)
// -----------------------------------------------------------------------------
//
// The SDL/OpenGL prototype historically kept shipyard pricing and upgrade
// application logic in apps/stellar_game/main.cpp.
//
// This module extracts the deterministic *business logic* into the core library
// so:
//  - shipyard pricing is consistent across UI, tools, and tests
//  - save/load semantics remain stable as ship progression is tuned
//  - game UI code focuses on presentation

struct ShipyardPriceModel {
  // Effective station fee rate (0..1). The game typically computes this as
  // applyRepToFee(st.feeRate, rep).
  double feeRateEff{0.0};

  // Player reputation with the station faction. Conventionally ~[-100,+100].
  // Used for trade-in value and (optionally) pricing knobs.
  double rep{0.0};
};

// Coarse tier for price scaling. Higher tiers are more sensitive to faction tech.
enum class ShipyardItemTier : core::u8 {
  Basic = 0,
  Advanced = 1,
  Exotic = 2,
};

// Simple tier classifier based on the item's base price.
// (Callers can override if they want more explicit control.)
ShipyardItemTier shipyardTierForPrice(double basePriceCr);

// Clamp rep into a reasonable stable range.
double clampRep(double rep);

// Trade-in multiplier for old gear/hulls.
// Returns a value in ~[0.50, 0.80] for rep in [-100, +100].
double tradeInMultiplier(double rep);

// Price a shipyard item with a given base cost, scaled by:
//  - station fees (feeRateEff)
//  - faction procedural profile (wealth/tech)
//  - the item's coarse tier
//
// This keeps shipyard prices deterministic but system-dependent.
double shipyardPrice(core::u64 universeSeed,
                     core::u32 factionId,
                     econ::StationType stationType,
                     ShipyardItemTier tier,
                     double basePriceCr,
                     const ShipyardPriceModel& model);

// -----------------------------------------------------------------------------
// Incremental upgrades
// -----------------------------------------------------------------------------

struct UpgradeQuote {
  bool ok{false};
  const char* reason{nullptr};

  double costCr{0.0};
  double delta{0.0};
  int nextMk{0};
};

UpgradeQuote quoteCargoRackUpgrade(core::u64 universeSeed,
                                  const Station& station,
                                  const ShipyardPriceModel& model);
bool applyCargoRackUpgrade(double& credits, double& cargoCapacityKg, const UpgradeQuote& q);

UpgradeQuote quotePassengerCabinUpgrade(core::u64 universeSeed,
                                       const Station& station,
                                       const ShipyardPriceModel& model);
bool applyPassengerCabinUpgrade(double& credits, int& passengerSeats, const UpgradeQuote& q);

UpgradeQuote quoteFuelTankUpgrade(core::u64 universeSeed,
                                 const Station& station,
                                 const ShipyardPriceModel& model);
bool applyFuelTankUpgrade(double& credits, double& fuelMax, double& fuel, const UpgradeQuote& q);

UpgradeQuote quoteFsdTuningUpgrade(core::u64 universeSeed,
                                  const Station& station,
                                  const ShipyardPriceModel& model);
bool applyFsdTuningUpgrade(double& credits, double& fsdRangeLy, const UpgradeQuote& q);

UpgradeQuote quoteSmuggleHoldUpgrade(core::u64 universeSeed,
                                    const Station& station,
                                    const ShipyardPriceModel& model,
                                    int currentMk);
bool applySmuggleHoldUpgrade(double& credits, int& smuggleHoldMk, const UpgradeQuote& q);

// -----------------------------------------------------------------------------
// Trade-in purchases (hulls/modules/weapons)
// -----------------------------------------------------------------------------

struct PurchaseQuote {
  bool ok{false};
  const char* reason{nullptr};

  double buyCostCr{0.0};
  double tradeInCr{0.0};
  double finalCostCr{0.0};
};

// Hull purchase quote.
// The shipyard uses the player's current cargo/fuel capacities as "installed"
// equipment and rescales them with the hull's cargo/fuel multipliers.
struct HullPurchaseQuote : PurchaseQuote {
  double newCargoCapacityKg{0.0};
  double newFuelMax{0.0};
  bool cargoFits{true};
};

HullPurchaseQuote quoteHullPurchase(core::u64 universeSeed,
                                   const Station& station,
                                   const ShipyardPriceModel& model,
                                   ShipHullClass currentHull,
                                   ShipHullClass newHull,
                                   double cargoCapacityKg,
                                   double fuelMax,
                                   double cargoMassKg);

// Generic quote for a purchase that trades in old equipment for some fraction
// of its local shipyard value.
PurchaseQuote quoteTradeInPurchase(core::u64 universeSeed,
                                  const Station& station,
                                  const ShipyardPriceModel& model,
                                  ShipyardItemTier oldTier,
                                  double oldBasePriceCr,
                                  ShipyardItemTier newTier,
                                  double newBasePriceCr);

// Weapon purchase quote (uses weaponDef prices; trade-in is based on current weapon).
PurchaseQuote quoteWeaponPurchase(core::u64 universeSeed,
                                 const Station& station,
                                 const ShipyardPriceModel& model,
                                 WeaponType currentWeapon,
                                 WeaponType newWeapon);

// Helper: apply a successful purchase quote to credits.
// Returns true on success (sufficient credits).
bool applyPurchase(double& credits, const PurchaseQuote& q);

} // namespace stellar::sim
