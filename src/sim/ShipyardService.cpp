#include "stellar/sim/ShipyardService.h"

#include "stellar/core/Clamp.h"

#include <algorithm>

namespace stellar::sim {

static double clampFactor(double x) {
  return std::clamp(x, 0.55, 1.75);
}

double clampRep(double rep) {
  // The rest of the codebase tends to treat reputation as roughly [-100,+100].
  // Clamp to keep trade-in/pricing stable even if callers temporarily exceed.
  return std::clamp(rep, -100.0, 100.0);
}

double tradeInMultiplier(double rep) {
  const double r = clampRep(rep);
  const double t = (r + 100.0) / 200.0; // 0..1
  // Low rep: the station offers poor trade-in.
  // High rep: the station is more generous.
  return 0.50 + 0.30 * t; // [0.50..0.80]
}

ShipyardItemTier shipyardTierForPrice(double basePriceCr) {
  if (basePriceCr <= 1e-6) return ShipyardItemTier::Basic;
  if (basePriceCr < 9000.0) return ShipyardItemTier::Advanced;
  return ShipyardItemTier::Exotic;
}

static double techFactorForTier(const FactionProfile& fp, ShipyardItemTier tier) {
  const double tech = std::clamp(fp.tech, 0.0, 1.0);
  switch (tier) {
    case ShipyardItemTier::Basic:
      // Basic items are everywhere; only a gentle tech effect.
      return 1.00 - 0.06 * tech;
    case ShipyardItemTier::Advanced:
      return 1.18 - 0.38 * tech;
    case ShipyardItemTier::Exotic:
    default:
      return 1.30 - 0.58 * tech;
  }
}

static double wealthFactor(const FactionProfile& fp) {
  const double w = std::clamp(fp.wealth, 0.0, 1.0);
  // Wealthy factions tend to have better supply chains.
  return 1.12 - 0.27 * w; // ~[0.85..1.12]
}

static double stationTypeFactor(econ::StationType t) {
  // Shipyards are assumed to be "hardware-centric" and slightly cheaper.
  // Other station types (if ever used) can bias prices up/down.
  switch (t) {
    case econ::StationType::Shipyard:     return 0.97;
    case econ::StationType::Research:     return 1.08;
    case econ::StationType::TradeHub:     return 1.02;
    case econ::StationType::Industrial:   return 0.99;
    case econ::StationType::Refinery:     return 1.00;
    case econ::StationType::Mining:       return 1.01;
    case econ::StationType::Agricultural: return 1.03;
    case econ::StationType::Outpost:      return 1.05;
    default:                              return 1.02;
  }
}

double shipyardPrice(core::u64 universeSeed,
                     core::u32 factionId,
                     econ::StationType stationType,
                     ShipyardItemTier tier,
                     double basePriceCr,
                     const ShipyardPriceModel& model) {
  if (basePriceCr <= 0.0) return 0.0;

  const FactionProfile fp = factionProfile(universeSeed, factionId);

  const double f = clampFactor(wealthFactor(fp) * techFactorForTier(fp, tier) * stationTypeFactor(stationType));
  const double feeMult = 1.0 + std::max(0.0, model.feeRateEff);
  return basePriceCr * f * feeMult;
}

static bool requireShipyard(const Station& station, UpgradeQuote& q) {
  if (station.type != econ::StationType::Shipyard) {
    q.ok = false;
    q.reason = "no_service";
    return false;
  }
  return true;
}

UpgradeQuote quoteCargoRackUpgrade(core::u64 universeSeed,
                                  const Station& station,
                                  const ShipyardPriceModel& model) {
  UpgradeQuote q{};
  if (!requireShipyard(station, q)) return q;
  q.ok = true;
  q.delta = 200.0; // kg
  q.costCr = shipyardPrice(universeSeed, station.factionId, station.type, ShipyardItemTier::Basic, 2000.0, model);
  return q;
}

bool applyCargoRackUpgrade(double& credits, double& cargoCapacityKg, const UpgradeQuote& q) {
  if (!q.ok || q.costCr <= 0.0) return false;
  if (credits + 1e-6 < q.costCr) return false;
  credits -= q.costCr;
  cargoCapacityKg += q.delta;
  return true;
}

UpgradeQuote quotePassengerCabinUpgrade(core::u64 universeSeed,
                                       const Station& station,
                                       const ShipyardPriceModel& model) {
  UpgradeQuote q{};
  if (!requireShipyard(station, q)) return q;
  q.ok = true;
  q.delta = 2.0; // seats
  q.costCr = shipyardPrice(universeSeed, station.factionId, station.type, ShipyardItemTier::Basic, 3000.0, model);
  return q;
}

bool applyPassengerCabinUpgrade(double& credits, int& passengerSeats, const UpgradeQuote& q) {
  if (!q.ok || q.costCr <= 0.0) return false;
  if (credits + 1e-6 < q.costCr) return false;
  credits -= q.costCr;
  passengerSeats += (int)std::llround(q.delta);
  return true;
}

UpgradeQuote quoteFuelTankUpgrade(core::u64 universeSeed,
                                 const Station& station,
                                 const ShipyardPriceModel& model) {
  UpgradeQuote q{};
  if (!requireShipyard(station, q)) return q;
  q.ok = true;
  q.delta = 10.0;
  q.costCr = shipyardPrice(universeSeed, station.factionId, station.type, ShipyardItemTier::Basic, 2500.0, model);
  return q;
}

bool applyFuelTankUpgrade(double& credits, double& fuelMax, double& fuel, const UpgradeQuote& q) {
  if (!q.ok || q.costCr <= 0.0) return false;
  if (credits + 1e-6 < q.costCr) return false;
  credits -= q.costCr;
  fuelMax += q.delta;
  fuel = std::min(fuel, fuelMax);
  return true;
}

UpgradeQuote quoteFsdTuningUpgrade(core::u64 universeSeed,
                                  const Station& station,
                                  const ShipyardPriceModel& model) {
  UpgradeQuote q{};
  if (!requireShipyard(station, q)) return q;
  q.ok = true;
  q.delta = 2.0;
  q.costCr = shipyardPrice(universeSeed, station.factionId, station.type, ShipyardItemTier::Exotic, 9000.0, model);
  return q;
}

bool applyFsdTuningUpgrade(double& credits, double& fsdRangeLy, const UpgradeQuote& q) {
  if (!q.ok || q.costCr <= 0.0) return false;
  if (credits + 1e-6 < q.costCr) return false;
  credits -= q.costCr;
  fsdRangeLy += q.delta;
  return true;
}

UpgradeQuote quoteSmuggleHoldUpgrade(core::u64 universeSeed,
                                    const Station& station,
                                    const ShipyardPriceModel& model,
                                    int currentMk) {
  UpgradeQuote q{};
  if (!requireShipyard(station, q)) return q;

  const int cur = std::clamp(currentMk, 0, 3);
  if (cur >= 3) {
    q.ok = false;
    q.reason = "max";
    return q;
  }

  static const double kBaseCost[4] = {0.0, 4500.0, 8500.0, 14000.0};
  const int next = cur + 1;
  q.ok = true;
  q.nextMk = next;

  // Smuggling gear is (in-fiction) black-market adjacent, so corruption makes it cheaper.
  const FactionProfile fp = factionProfile(universeSeed, station.factionId);
  const double corruption = std::clamp(fp.corruption, 0.0, 1.0);
  const double corruptFactor = std::clamp(1.18 - 0.40 * corruption, 0.70, 1.25);

  q.costCr = shipyardPrice(universeSeed,
                           station.factionId,
                           station.type,
                           ShipyardItemTier::Advanced,
                           kBaseCost[next],
                           model) * corruptFactor;
  return q;
}

bool applySmuggleHoldUpgrade(double& credits, int& smuggleHoldMk, const UpgradeQuote& q) {
  if (!q.ok || q.costCr <= 0.0) return false;
  if (credits + 1e-6 < q.costCr) return false;
  if (q.nextMk <= 0) return false;
  credits -= q.costCr;
  smuggleHoldMk = std::clamp(q.nextMk, 0, 3);
  return true;
}

HullPurchaseQuote quoteHullPurchase(core::u64 universeSeed,
                                   const Station& station,
                                   const ShipyardPriceModel& model,
                                   ShipHullClass currentHull,
                                   ShipHullClass newHull,
                                   double cargoCapacityKg,
                                   double fuelMax,
                                   double cargoMassKg) {
  HullPurchaseQuote q{};
  if (station.type != econ::StationType::Shipyard) {
    q.ok = false;
    q.reason = "no_service";
    return q;
  }

  if (currentHull == newHull) {
    q.ok = false;
    q.reason = "same";
    return q;
  }

  const HullDef& cur = hullDef(currentHull);
  const HullDef& dst = hullDef(newHull);

  const double newCargo = cargoCapacityKg * (dst.cargoMult / std::max(1e-6, cur.cargoMult));
  const double newFuel = fuelMax * (dst.fuelMult / std::max(1e-6, cur.fuelMult));
  q.newCargoCapacityKg = newCargo;
  q.newFuelMax = newFuel;
  q.cargoFits = (cargoMassKg <= newCargo + 1e-6);

  const ShipyardItemTier tierCur = shipyardTierForPrice(cur.priceCr);
  const ShipyardItemTier tierNew = shipyardTierForPrice(dst.priceCr);

  const PurchaseQuote base = quoteTradeInPurchase(universeSeed,
                                                  station,
                                                  model,
                                                  tierCur,
                                                  cur.priceCr,
                                                  tierNew,
                                                  dst.priceCr);
  q.ok = base.ok;
  q.reason = base.reason;
  q.buyCostCr = base.buyCostCr;
  q.tradeInCr = base.tradeInCr;
  q.finalCostCr = base.finalCostCr;
  return q;
}

PurchaseQuote quoteTradeInPurchase(core::u64 universeSeed,
                                  const Station& station,
                                  const ShipyardPriceModel& model,
                                  ShipyardItemTier oldTier,
                                  double oldBasePriceCr,
                                  ShipyardItemTier newTier,
                                  double newBasePriceCr) {
  PurchaseQuote q{};
  if (station.type != econ::StationType::Shipyard) {
    q.ok = false;
    q.reason = "no_service";
    return q;
  }

  if (newBasePriceCr < 0.0) {
    q.ok = false;
    q.reason = "invalid";
    return q;
  }

  q.ok = true;

  q.buyCostCr = shipyardPrice(universeSeed,
                              station.factionId,
                              station.type,
                              newTier,
                              newBasePriceCr,
                              model);

  // Trade-in is based on the *local* value of the old item.
  const double oldValue = shipyardPrice(universeSeed,
                                        station.factionId,
                                        station.type,
                                        oldTier,
                                        std::max(0.0, oldBasePriceCr),
                                        model);

  q.tradeInCr = oldValue * tradeInMultiplier(model.rep);

  // Never pay the player for "upgrades"; clamp at 0.
  q.finalCostCr = std::max(0.0, q.buyCostCr - q.tradeInCr);
  return q;
}

PurchaseQuote quoteWeaponPurchase(core::u64 universeSeed,
                                 const Station& station,
                                 const ShipyardPriceModel& model,
                                 WeaponType currentWeapon,
                                 WeaponType newWeapon) {
  PurchaseQuote q{};
  if (station.type != econ::StationType::Shipyard) {
    q.ok = false;
    q.reason = "no_service";
    return q;
  }

  if (currentWeapon == newWeapon) {
    q.ok = false;
    q.reason = "same";
    return q;
  }

  const double oldPrice = weaponDef(currentWeapon).priceCr;
  const double newPrice = weaponDef(newWeapon).priceCr;
  if (newPrice < 0.0) {
    q.ok = false;
    q.reason = "invalid";
    return q;
  }

  return quoteTradeInPurchase(universeSeed,
                              station,
                              model,
                              shipyardTierForPrice(oldPrice),
                              oldPrice,
                              shipyardTierForPrice(newPrice),
                              newPrice);
}

bool applyPurchase(double& credits, const PurchaseQuote& q) {
  if (!q.ok) return false;
  if (q.finalCostCr <= 0.0) return true;
  if (credits + 1e-6 < q.finalCostCr) return false;
  credits -= q.finalCostCr;
  return true;
}

} // namespace stellar::sim
