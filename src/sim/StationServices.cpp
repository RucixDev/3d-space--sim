#include "stellar/sim/StationServices.h"

#include "stellar/econ/Commodity.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::sim {

namespace {

static constexpr std::size_t idx(econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

static double clamp01Safe(double v) {
  if (!std::isfinite(v)) return 0.0;
  return std::clamp(v, 0.0, 1.0);
}

static double budgetOrInf(double creditsBudgetCr) {
  if (!std::isfinite(creditsBudgetCr) || creditsBudgetCr < 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  return std::max(0.0, creditsBudgetCr);
}

static econ::MarketQuote q(const econ::StationEconomyState& stEcon,
                           const econ::StationEconomyModel& model,
                           econ::CommodityId id,
                           const StationServicePriceModel& pm) {
  return econ::quote(stEcon, model, id, clamp01Safe(pm.bidAskSpread));
}

struct CmRecipe {
  double weapons{0.0};
  double electronics{0.0};
  double metals{0.0};
  double machinery{0.0};
  double utility{0.0};
};

static constexpr CmRecipe kFlareRecipe{0.50, 0.10, 0.0, 0.0, 1.0};
static constexpr CmRecipe kChaffRecipe{0.55, 0.15, 0.0, 0.0, 1.0};
// Heat sinks are rarer and expensive; recipe is intentionally heavier.
static constexpr CmRecipe kHeatSinkRecipe{0.60, 2.00, 2.00, 1.00, 3.0};

struct MissileRecipe {
  double weapons{0.0};
  double electronics{0.0};
  double utility{0.0};
};

static MissileRecipe missileRecipeForWeapon(WeaponType w) {
  if (!weaponUsesAmmo(w)) return {};
  // Utility is proportional to nominal damage; this is a cheap heuristic for "best use" of scarce ammo.
  const double dmg = std::max(0.0, weaponDef(w).dmg);

  switch (w) {
    case WeaponType::HomingMissile:
      return {1.00, 0.00, dmg};
    case WeaponType::RadarMissile:
      return {0.85, 0.55, dmg};
    default:
      return {};
  }
}

static double costForCommodities(double weaponsUnits,
                                double electronicsUnits,
                                double metalsUnits,
                                double machineryUnits,
                                double weaponsAsk,
                                double electronicsAsk,
                                double metalsAsk,
                                double machineryAsk,
                                double feeRateEff) {
  const double fee = clamp01Safe(feeRateEff);
  const double base =
      std::max(0.0, weaponsUnits) * std::max(0.0, weaponsAsk) +
      std::max(0.0, electronicsUnits) * std::max(0.0, electronicsAsk) +
      std::max(0.0, metalsUnits) * std::max(0.0, metalsAsk) +
      std::max(0.0, machineryUnits) * std::max(0.0, machineryAsk);
  return base * (1.0 + fee);
}

static int clampInt(int v, int lo, int hi) {
  return std::max(lo, std::min(v, hi));
}

} // namespace

HullRepairQuote quoteHullRepairToFull(const econ::StationEconomyState& stEcon,
                                     const econ::StationEconomyModel& model,
                                     double hullCurrent,
                                     double hullMax,
                                     const StationServicePriceModel& pm,
                                     double creditsBudgetCr) {
  HullRepairQuote out{};

  if (!std::isfinite(hullCurrent) || !std::isfinite(hullMax) || hullMax <= 0.0) {
    out.ok = false;
    out.reason = "invalid";
    return out;
  }

  const double missing = std::max(0.0, hullMax - std::max(0.0, hullCurrent));
  out.hullMissing = missing;

  if (missing <= 1e-6) {
    out.ok = false;
    out.reason = "no_need";
    return out;
  }

  if (kRepairMetalsPerHull <= 0.0 || kRepairMachineryPerHull <= 0.0) {
    out.ok = false;
    out.reason = "bad_recipe";
    return out;
  }

  const auto qMet = q(stEcon, model, econ::CommodityId::Metals, pm);
  const auto qMach = q(stEcon, model, econ::CommodityId::Machinery, pm);

  out.metalsAsk = qMet.ask;
  out.machineryAsk = qMach.ask;

  const double maxByMetals = (qMet.inventory <= 0.0) ? 0.0 : (qMet.inventory / kRepairMetalsPerHull);
  const double maxByMach = (qMach.inventory <= 0.0) ? 0.0 : (qMach.inventory / kRepairMachineryPerHull);

  double hullByStock = std::min({missing, maxByMetals, maxByMach});
  if (!std::isfinite(hullByStock) || hullByStock <= 1e-6) {
    out.ok = false;
    out.reason = "out_of_stock";
    out.limitedByStock = true;
    return out;
  }

  out.limitedByStock = (hullByStock + 1e-6 < missing);

  const double fee = clamp01Safe(pm.feeRateEff);
  const double costPerHull =
      (kRepairMetalsPerHull * qMet.ask + kRepairMachineryPerHull * qMach.ask) * (1.0 + fee);

  double hull = hullByStock;
  const double budget = budgetOrInf(creditsBudgetCr);
  out.limitedByCredits = false;
  if (std::isfinite(budget) && budget < std::numeric_limits<double>::infinity()) {
    if (costPerHull > 1e-9) {
      const double hullByCredits = budget / costPerHull;
      if (hullByCredits + 1e-9 < hull) {
        hull = std::max(0.0, hullByCredits);
        out.limitedByCredits = true;
      }
    } else {
      // Degenerate pricing.
      out.limitedByCredits = false;
    }
  }

  hull = std::max(0.0, std::min(hull, missing));

  if (hull <= 1e-6) {
    out.ok = false;
    out.reason = out.limitedByCredits ? "need_credits" : "out_of_stock";
    return out;
  }

  out.hullToRepair = hull;
  out.metalsUnits = hull * kRepairMetalsPerHull;
  out.machineryUnits = hull * kRepairMachineryPerHull;
  out.costCr = hull * costPerHull;

  out.ok = true;
  return out;
}

HullRepairResult applyHullRepairToFull(econ::StationEconomyState& stEcon,
                                      const econ::StationEconomyModel& model,
                                      double& ioCredits,
                                      double& ioHullCurrent,
                                      double hullMax,
                                      const StationServicePriceModel& pm) {
  HullRepairResult r{};

  const auto qRep = quoteHullRepairToFull(stEcon, model, ioHullCurrent, hullMax, pm, ioCredits);
  if (!qRep.ok) {
    r.ok = false;
    r.reason = qRep.reason;
    return r;
  }

  const double fee = clamp01Safe(pm.feeRateEff);
  const double costPerHull =
      (kRepairMetalsPerHull * qRep.metalsAsk + kRepairMachineryPerHull * qRep.machineryAsk) * (1.0 + fee);

  // Consume commodities.
  const double metalsTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Metals, qRep.metalsUnits);
  const double machTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Machinery, qRep.machineryUnits);

  // Compute actual repaired amount (defensive against any unexpected inventory clamping).
  double hullRepaired = qRep.hullToRepair;
  if (kRepairMetalsPerHull > 0.0) {
    hullRepaired = std::min(hullRepaired, metalsTaken / kRepairMetalsPerHull);
  }
  if (kRepairMachineryPerHull > 0.0) {
    hullRepaired = std::min(hullRepaired, machTaken / kRepairMachineryPerHull);
  }

  if (!std::isfinite(hullRepaired) || hullRepaired <= 1e-6) {
    r.ok = false;
    r.reason = "out_of_stock";
    return r;
  }

  hullRepaired = std::max(0.0, std::min(hullRepaired, qRep.hullToRepair));
  const double cost = hullRepaired * costPerHull;

  if (ioCredits + 1e-9 < cost) {
    // Should not happen (qRep already considered credit budget), but clamp defensively.
    const double maxHullByCredits = (costPerHull > 1e-9) ? std::max(0.0, ioCredits / costPerHull) : 0.0;
    hullRepaired = std::min(hullRepaired, maxHullByCredits);
  }

  const double costPaid = std::min(std::max(0.0, ioCredits), hullRepaired * costPerHull);
  ioCredits -= costPaid;

  ioHullCurrent = std::min(hullMax, std::max(0.0, ioHullCurrent) + hullRepaired);

  r.ok = true;
  r.hullRepaired = hullRepaired;
  r.creditsPaid = costPaid;
  r.metalsTaken = metalsTaken;
  r.machineryTaken = machTaken;
  return r;
}

RefuelQuote quoteRefuelToFull(const econ::StationEconomyState& stEcon,
                             const econ::StationEconomyModel& model,
                             double fuelCurrent,
                             double fuelMax,
                             const StationServicePriceModel& pm,
                             double creditsBudgetCr) {
  RefuelQuote out{};

  if (!std::isfinite(fuelCurrent) || !std::isfinite(fuelMax) || fuelMax < 0.0) {
    out.ok = false;
    out.reason = "invalid";
    return out;
  }

  const double need = std::max(0.0, fuelMax - std::max(0.0, fuelCurrent));
  out.fuelMissing = need;

  if (need <= 1e-6) {
    out.ok = false;
    out.reason = "no_need";
    return out;
  }

  const auto qFuel = q(stEcon, model, econ::CommodityId::Fuel, pm);
  out.fuelAsk = qFuel.ask;
  out.fuelInventory = qFuel.inventory;

  double byStock = std::min(need, std::max(0.0, qFuel.inventory));
  if (byStock <= 1e-6) {
    out.ok = false;
    out.reason = "out_of_stock";
    out.limitedByStock = true;
    return out;
  }

  out.limitedByStock = (byStock + 1e-6 < need);

  const double fee = clamp01Safe(pm.feeRateEff);
  const double costPerUnit = qFuel.ask * (1.0 + fee);

  double buy = byStock;
  const double budget = budgetOrInf(creditsBudgetCr);
  out.limitedByCredits = false;
  if (std::isfinite(budget) && budget < std::numeric_limits<double>::infinity()) {
    if (costPerUnit > 1e-9) {
      const double byCredits = budget / costPerUnit;
      if (byCredits + 1e-9 < buy) {
        buy = std::max(0.0, byCredits);
        out.limitedByCredits = true;
      }
    }
  }

  buy = std::max(0.0, std::min(buy, need));

  if (buy <= 1e-6) {
    out.ok = false;
    out.reason = out.limitedByCredits ? "need_credits" : "out_of_stock";
    return out;
  }

  out.fuelToBuy = buy;
  out.costCr = buy * costPerUnit;

  out.ok = true;
  return out;
}

RefuelResult applyRefuelToFull(econ::StationEconomyState& stEcon,
                              const econ::StationEconomyModel& model,
                              double& ioCredits,
                              double& ioFuelCurrent,
                              double fuelMax,
                              const StationServicePriceModel& pm) {
  RefuelResult r{};

  const auto qRef = quoteRefuelToFull(stEcon, model, ioFuelCurrent, fuelMax, pm, ioCredits);
  if (!qRef.ok) {
    r.ok = false;
    r.reason = qRef.reason;
    return r;
  }

  const double fee = clamp01Safe(pm.feeRateEff);
  const double costPerUnit = qRef.fuelAsk * (1.0 + fee);

  const double taken = econ::takeInventory(stEcon, model, econ::CommodityId::Fuel, qRef.fuelToBuy);
  if (taken <= 1e-6) {
    r.ok = false;
    r.reason = "out_of_stock";
    return r;
  }

  const double cost = taken * costPerUnit;
  if (ioCredits + 1e-9 < cost) {
    // Defensive: should not happen if quote used current credits.
    r.ok = false;
    r.reason = "need_credits";
    return r;
  }

  ioCredits -= cost;
  ioFuelCurrent = std::min(fuelMax, std::max(0.0, ioFuelCurrent) + taken);

  r.ok = true;
  r.fuelBought = taken;
  r.creditsPaid = cost;
  return r;
}

namespace {

struct CmPlan {
  int buyF{0};
  int buyC{0};
  int buyH{0};
  double useW{0.0};
  double useE{0.0};
  double useMet{0.0};
  double useMach{0.0};
  double cost{0.0};
  double util{0.0};
};

static CmPlan bestCountermeasurePlan(const econ::StationEconomyState& stEcon,
                                    const econ::StationEconomyModel& model,
                                    int needF,
                                    int needC,
                                    int needH,
                                    const StationServicePriceModel& pm,
                                    double creditsBudgetCr) {
  CmPlan best{};

  needF = std::max(0, needF);
  needC = std::max(0, needC);
  needH = std::max(0, needH);

  const auto qW = q(stEcon, model, econ::CommodityId::Weapons, pm);
  const auto qE = q(stEcon, model, econ::CommodityId::Electronics, pm);
  const auto qMet = q(stEcon, model, econ::CommodityId::Metals, pm);
  const auto qMach = q(stEcon, model, econ::CommodityId::Machinery, pm);

  const double fee = clamp01Safe(pm.feeRateEff);
  const double budget = budgetOrInf(creditsBudgetCr);

  const double invW = std::max(0.0, qW.inventory);
  const double invE = std::max(0.0, qE.inventory);
  const double invMet = std::max(0.0, qMet.inventory);
  const double invMach = std::max(0.0, qMach.inventory);

  auto consider = [&](int f, int c, int h) {
    if (f < 0 || c < 0 || h < 0) return;
    if (f > needF || c > needC || h > needH) return;

    const double useW = f * kFlareRecipe.weapons + c * kChaffRecipe.weapons + h * kHeatSinkRecipe.weapons;
    const double useE = f * kFlareRecipe.electronics + c * kChaffRecipe.electronics + h * kHeatSinkRecipe.electronics;
    const double useMet = h * kHeatSinkRecipe.metals;
    const double useMach = h * kHeatSinkRecipe.machinery;

    if (useW > invW + 1e-9) return;
    if (useE > invE + 1e-9) return;
    if (useMet > invMet + 1e-9) return;
    if (useMach > invMach + 1e-9) return;

    const double cost = costForCommodities(useW, useE, useMet, useMach,
                                          qW.ask, qE.ask, qMet.ask, qMach.ask, fee);
    if (cost > budget + 1e-9) return;

    const double util = f * kFlareRecipe.utility + c * kChaffRecipe.utility + h * kHeatSinkRecipe.utility;
    const int total = f + c + h;
    const int bestTotal = best.buyF + best.buyC + best.buyH;

    if (util > best.util + 1e-9 ||
        (std::abs(util - best.util) <= 1e-9 && total > bestTotal) ||
        (std::abs(util - best.util) <= 1e-9 && total == bestTotal && cost + 1e-9 < best.cost) ||
        (std::abs(util - best.util) <= 1e-9 && total == bestTotal && std::abs(cost - best.cost) <= 1e-9 && h > best.buyH)) {
      best.buyF = f;
      best.buyC = c;
      best.buyH = h;
      best.useW = useW;
      best.useE = useE;
      best.useMet = useMet;
      best.useMach = useMach;
      best.cost = cost;
      best.util = util;
    }
  };

  // Brute force is tiny (caps are ~ <= 10).
  for (int h = 0; h <= needH; ++h) {
    for (int c = 0; c <= needC; ++c) {
      for (int f = 0; f <= needF; ++f) {
        consider(f, c, h);
      }
    }
  }

  return best;
}

static CountermeasureRestockQuote buildCountermeasureQuote(const econ::StationEconomyState& stEcon,
                                                          const econ::StationEconomyModel& model,
                                                          ShipHullClass hullClass,
                                                          int haveF,
                                                          int haveC,
                                                          int haveH,
                                                          int needF,
                                                          int needC,
                                                          int needH,
                                                          const StationServicePriceModel& pm,
                                                          double creditsBudgetCr) {
  CountermeasureRestockQuote out{};
  out.caps = countermeasureCapsForHull(hullClass);

  // Defensive clamp.
  haveF = clampInt(haveF, 0, out.caps.flares);
  haveC = clampInt(haveC, 0, out.caps.chaff);
  haveH = clampInt(haveH, 0, out.caps.heatSinks);

  out.haveFlares = haveF;
  out.haveChaff = haveC;
  out.haveHeatSinks = haveH;

  out.needFlares = std::max(0, needF);
  out.needChaff = std::max(0, needC);
  out.needHeatSinks = std::max(0, needH);

  if (out.needFlares + out.needChaff + out.needHeatSinks <= 0) {
    out.ok = false;
    out.reason = "no_need";
    return out;
  }

  const CmPlan plan = bestCountermeasurePlan(stEcon, model,
                                            out.needFlares, out.needChaff, out.needHeatSinks,
                                            pm, creditsBudgetCr);

  out.buyFlares = plan.buyF;
  out.buyChaff = plan.buyC;
  out.buyHeatSinks = plan.buyH;

  out.useWeapons = plan.useW;
  out.useElectronics = plan.useE;
  out.useMetals = plan.useMet;
  out.useMachinery = plan.useMach;
  out.costCr = plan.cost;

  if (out.buyFlares + out.buyChaff + out.buyHeatSinks <= 0) {
    // Distinguish stock vs credit failure (best-effort).
    const CmPlan planStockOnly = bestCountermeasurePlan(stEcon, model,
                                                       out.needFlares, out.needChaff, out.needHeatSinks,
                                                       pm, -1.0);
    if (planStockOnly.buyF + planStockOnly.buyC + planStockOnly.buyH <= 0) {
      out.ok = false;
      out.reason = "out_of_stock";
      out.limitedByStock = true;
      return out;
    }

    out.ok = false;
    out.reason = "need_credits";
    out.limitedByCredits = true;
    return out;
  }

  // Determine limiting factors.
  const CmPlan planStockOnly = bestCountermeasurePlan(stEcon, model,
                                                     out.needFlares, out.needChaff, out.needHeatSinks,
                                                     pm, -1.0);

  out.limitedByStock = (planStockOnly.buyF < out.needFlares) ||
                       (planStockOnly.buyC < out.needChaff) ||
                       (planStockOnly.buyH < out.needHeatSinks);

  const double budget = budgetOrInf(creditsBudgetCr);
  if (std::isfinite(budget) && budget < std::numeric_limits<double>::infinity()) {
    out.limitedByCredits = (plan.buyF < planStockOnly.buyF) ||
                           (plan.buyC < planStockOnly.buyC) ||
                           (plan.buyH < planStockOnly.buyH);
  } else {
    out.limitedByCredits = false;
  }

  out.ok = true;
  return out;
}

} // namespace

CountermeasureRestockQuote quoteCountermeasureRestockAll(const econ::StationEconomyState& stEcon,
                                                        const econ::StationEconomyModel& model,
                                                        ShipHullClass hullClass,
                                                        int haveFlares,
                                                        int haveChaff,
                                                        int haveHeatSinks,
                                                        const StationServicePriceModel& pm,
                                                        double creditsBudgetCr) {
  const auto cap = countermeasureCapsForHull(hullClass);
  haveFlares = clampInt(haveFlares, 0, cap.flares);
  haveChaff = clampInt(haveChaff, 0, cap.chaff);
  haveHeatSinks = clampInt(haveHeatSinks, 0, cap.heatSinks);

  const int needF = cap.flares - haveFlares;
  const int needC = cap.chaff - haveChaff;
  const int needH = cap.heatSinks - haveHeatSinks;

  return buildCountermeasureQuote(stEcon, model, hullClass,
                                 haveFlares, haveChaff, haveHeatSinks,
                                 needF, needC, needH,
                                 pm, creditsBudgetCr);
}

CountermeasureRestockQuote quoteCountermeasureRestockFlaresToCap(const econ::StationEconomyState& stEcon,
                                                                 const econ::StationEconomyModel& model,
                                                                 ShipHullClass hullClass,
                                                                 int haveFlares,
                                                                 const StationServicePriceModel& pm,
                                                                 double creditsBudgetCr) {
  const auto cap = countermeasureCapsForHull(hullClass);
  haveFlares = clampInt(haveFlares, 0, cap.flares);
  const int needF = cap.flares - haveFlares;

  // Keep other caps/fields present but not targeted.
  return buildCountermeasureQuote(stEcon, model, hullClass,
                                 haveFlares, 0, 0,
                                 needF, 0, 0,
                                 pm, creditsBudgetCr);
}

CountermeasureRestockQuote quoteCountermeasureRestockChaffToCap(const econ::StationEconomyState& stEcon,
                                                                const econ::StationEconomyModel& model,
                                                                ShipHullClass hullClass,
                                                                int haveChaff,
                                                                const StationServicePriceModel& pm,
                                                                double creditsBudgetCr) {
  const auto cap = countermeasureCapsForHull(hullClass);
  haveChaff = clampInt(haveChaff, 0, cap.chaff);
  const int needC = cap.chaff - haveChaff;

  return buildCountermeasureQuote(stEcon, model, hullClass,
                                 0, haveChaff, 0,
                                 0, needC, 0,
                                 pm, creditsBudgetCr);
}

CountermeasureRestockQuote quoteCountermeasureRestockHeatSinksToCap(const econ::StationEconomyState& stEcon,
                                                                    const econ::StationEconomyModel& model,
                                                                    ShipHullClass hullClass,
                                                                    int haveHeatSinks,
                                                                    const StationServicePriceModel& pm,
                                                                    double creditsBudgetCr) {
  const auto cap = countermeasureCapsForHull(hullClass);
  haveHeatSinks = clampInt(haveHeatSinks, 0, cap.heatSinks);
  const int needH = cap.heatSinks - haveHeatSinks;

  return buildCountermeasureQuote(stEcon, model, hullClass,
                                 0, 0, haveHeatSinks,
                                 0, 0, needH,
                                 pm, creditsBudgetCr);
}

CountermeasureRestockResult applyCountermeasureRestockAll(econ::StationEconomyState& stEcon,
                                                         const econ::StationEconomyModel& model,
                                                         double& ioCredits,
                                                         ShipHullClass hullClass,
                                                         int& ioFlares,
                                                         int& ioChaff,
                                                         int& ioHeatSinks,
                                                         const StationServicePriceModel& pm) {
  CountermeasureRestockResult r{};

  const auto qRestock = quoteCountermeasureRestockAll(stEcon, model, hullClass,
                                                     ioFlares, ioChaff, ioHeatSinks,
                                                     pm, ioCredits);
  if (!qRestock.ok) {
    r.ok = false;
    r.reason = qRestock.reason;
    return r;
  }

  // Consume commodities.
  const double wTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Weapons, qRestock.useWeapons);
  const double eTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Electronics, qRestock.useElectronics);
  const double metTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Metals, qRestock.useMetals);
  const double machTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Machinery, qRestock.useMachinery);

  if (ioCredits + 1e-9 < qRestock.costCr) {
    // Should not happen (quote already considered credits).
    r.ok = false;
    r.reason = "need_credits";
    return r;
  }

  ioCredits -= qRestock.costCr;

  ioFlares += qRestock.buyFlares;
  ioChaff += qRestock.buyChaff;
  ioHeatSinks += qRestock.buyHeatSinks;

  clampCountermeasureAmmo(ioFlares, ioChaff, ioHeatSinks, hullClass);

  r.ok = true;
  r.flaresBought = qRestock.buyFlares;
  r.chaffBought = qRestock.buyChaff;
  r.heatSinksBought = qRestock.buyHeatSinks;
  r.creditsPaid = qRestock.costCr;
  r.weaponsTaken = wTaken;
  r.electronicsTaken = eTaken;
  r.metalsTaken = metTaken;
  r.machineryTaken = machTaken;
  return r;
}

CountermeasureRestockResult applyCountermeasureRestockFlaresToCap(econ::StationEconomyState& stEcon,
                                                                  const econ::StationEconomyModel& model,
                                                                  double& ioCredits,
                                                                  ShipHullClass hullClass,
                                                                  int& ioFlares,
                                                                  const StationServicePriceModel& pm) {
  int dummyC = 0;
  int dummyH = 0;
  auto q = quoteCountermeasureRestockFlaresToCap(stEcon, model, hullClass, ioFlares, pm, ioCredits);
  if (!q.ok) return {false, q.reason, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // Consume.
  const double wTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Weapons, q.useWeapons);
  const double eTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Electronics, q.useElectronics);
  const double metTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Metals, q.useMetals);
  const double machTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Machinery, q.useMachinery);

  if (ioCredits + 1e-9 < q.costCr) return {false, "need_credits"};
  ioCredits -= q.costCr;
  ioFlares += q.buyFlares;
  clampCountermeasureAmmo(ioFlares, dummyC, dummyH, hullClass);

  CountermeasureRestockResult r{};
  r.ok = true;
  r.flaresBought = q.buyFlares;
  r.creditsPaid = q.costCr;
  r.weaponsTaken = wTaken;
  r.electronicsTaken = eTaken;
  r.metalsTaken = metTaken;
  r.machineryTaken = machTaken;
  return r;
}

CountermeasureRestockResult applyCountermeasureRestockChaffToCap(econ::StationEconomyState& stEcon,
                                                                 const econ::StationEconomyModel& model,
                                                                 double& ioCredits,
                                                                 ShipHullClass hullClass,
                                                                 int& ioChaff,
                                                                 const StationServicePriceModel& pm) {
  int dummyF = 0;
  int dummyH = 0;
  auto q = quoteCountermeasureRestockChaffToCap(stEcon, model, hullClass, ioChaff, pm, ioCredits);
  if (!q.ok) return {false, q.reason, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0};

  const double wTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Weapons, q.useWeapons);
  const double eTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Electronics, q.useElectronics);
  const double metTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Metals, q.useMetals);
  const double machTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Machinery, q.useMachinery);

  if (ioCredits + 1e-9 < q.costCr) return {false, "need_credits"};
  ioCredits -= q.costCr;
  ioChaff += q.buyChaff;
  clampCountermeasureAmmo(dummyF, ioChaff, dummyH, hullClass);

  CountermeasureRestockResult r{};
  r.ok = true;
  r.chaffBought = q.buyChaff;
  r.creditsPaid = q.costCr;
  r.weaponsTaken = wTaken;
  r.electronicsTaken = eTaken;
  r.metalsTaken = metTaken;
  r.machineryTaken = machTaken;
  return r;
}

CountermeasureRestockResult applyCountermeasureRestockHeatSinksToCap(econ::StationEconomyState& stEcon,
                                                                     const econ::StationEconomyModel& model,
                                                                     double& ioCredits,
                                                                     ShipHullClass hullClass,
                                                                     int& ioHeatSinks,
                                                                     const StationServicePriceModel& pm) {
  int dummyF = 0;
  int dummyC = 0;
  auto q = quoteCountermeasureRestockHeatSinksToCap(stEcon, model, hullClass, ioHeatSinks, pm, ioCredits);
  if (!q.ok) return {false, q.reason, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0};

  const double wTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Weapons, q.useWeapons);
  const double eTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Electronics, q.useElectronics);
  const double metTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Metals, q.useMetals);
  const double machTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Machinery, q.useMachinery);

  if (ioCredits + 1e-9 < q.costCr) return {false, "need_credits"};
  ioCredits -= q.costCr;
  ioHeatSinks += q.buyHeatSinks;
  clampCountermeasureAmmo(dummyF, dummyC, ioHeatSinks, hullClass);

  CountermeasureRestockResult r{};
  r.ok = true;
  r.heatSinksBought = q.buyHeatSinks;
  r.creditsPaid = q.costCr;
  r.weaponsTaken = wTaken;
  r.electronicsTaken = eTaken;
  r.metalsTaken = metTaken;
  r.machineryTaken = machTaken;
  return r;
}

OrdnanceRearmQuote quoteOrdnanceRearmToFull(const econ::StationEconomyState& stEcon,
                                           const econ::StationEconomyModel& model,
                                           ShipHullClass hullClass,
                                           WeaponType weapon,
                                           core::u8 ammoCurrent,
                                           const StationServicePriceModel& pm,
                                           double creditsBudgetCr) {
  OrdnanceRearmQuote out{};
  out.weapon = weapon;

  const int ammoMax = weaponAmmoMax(weapon, hullClass);
  out.ammoMax = ammoMax;

  if (ammoMax <= 0) {
    out.ok = false;
    out.reason = "no_ammo_weapon";
    return out;
  }

  const int have = clampInt((int)ammoCurrent, 0, ammoMax);
  const int need = ammoMax - have;

  out.haveAmmo = have;
  out.needAmmo = need;

  if (need <= 0) {
    out.ok = false;
    out.reason = "no_need";
    return out;
  }

  const MissileRecipe rec = missileRecipeForWeapon(weapon);
  if (rec.weapons <= 0.0 && rec.electronics <= 0.0) {
    out.ok = false;
    out.reason = "bad_recipe";
    return out;
  }

  const auto qW = q(stEcon, model, econ::CommodityId::Weapons, pm);
  const auto qE = q(stEcon, model, econ::CommodityId::Electronics, pm);

  const double invW = std::max(0.0, qW.inventory);
  const double invE = std::max(0.0, qE.inventory);

  double maxByStock = std::numeric_limits<double>::infinity();
  if (rec.weapons > 0.0) maxByStock = std::min(maxByStock, invW / rec.weapons);
  if (rec.electronics > 0.0) maxByStock = std::min(maxByStock, invE / rec.electronics);

  if (!std::isfinite(maxByStock)) maxByStock = (double)need;

  const int maxByStockInt = std::max(0, (int)std::floor(maxByStock + 1e-9));

  const double fee = clamp01Safe(pm.feeRateEff);
  const double costPer = costForCommodities(rec.weapons, rec.electronics, 0.0, 0.0,
                                            qW.ask, qE.ask, 0.0, 0.0, fee);

  const double budget = budgetOrInf(creditsBudgetCr);
  int maxByCreditsInt = need;
  if (std::isfinite(budget) && budget < std::numeric_limits<double>::infinity()) {
    if (costPer > 1e-9) {
      maxByCreditsInt = std::max(0, (int)std::floor(budget / costPer + 1e-9));
    } else {
      maxByCreditsInt = need;
    }
  }

  const int buy = std::min({need, maxByStockInt, maxByCreditsInt});

  out.buyAmmo = buy;
  out.useWeapons = buy * rec.weapons;
  out.useElectronics = buy * rec.electronics;
  out.costCr = buy * costPer;

  out.limitedByStock = (maxByStockInt < need);
  out.limitedByCredits = (std::isfinite(budget) && budget < std::numeric_limits<double>::infinity()) && (buy < std::min(need, maxByStockInt));

  if (buy <= 0) {
    // Try to classify.
    if (maxByStockInt <= 0) {
      out.ok = false;
      out.reason = "out_of_stock";
    } else {
      out.ok = false;
      out.reason = "need_credits";
    }
    return out;
  }

  out.ok = true;
  return out;
}

OrdnanceRearmResult applyOrdnanceRearmToFull(econ::StationEconomyState& stEcon,
                                            const econ::StationEconomyModel& model,
                                            double& ioCredits,
                                            ShipHullClass hullClass,
                                            WeaponType weapon,
                                            core::u8& ioAmmo,
                                            const StationServicePriceModel& pm) {
  OrdnanceRearmResult r{};

  const auto qRear = quoteOrdnanceRearmToFull(stEcon, model, hullClass, weapon, ioAmmo, pm, ioCredits);
  if (!qRear.ok) {
    r.ok = false;
    r.reason = qRear.reason;
    return r;
  }

  if (ioCredits + 1e-9 < qRear.costCr) {
    r.ok = false;
    r.reason = "need_credits";
    return r;
  }

  const double wTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Weapons, qRear.useWeapons);
  const double eTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Electronics, qRear.useElectronics);

  ioCredits -= qRear.costCr;

  const int maxAmmo = weaponAmmoMax(weapon, hullClass);
  const int have = clampInt((int)ioAmmo, 0, maxAmmo);
  const int next = clampInt(have + qRear.buyAmmo, 0, maxAmmo);
  ioAmmo = (core::u8)next;

  r.ok = true;
  r.ammoBought = qRear.buyAmmo;
  r.creditsPaid = qRear.costCr;
  r.weaponsTaken = wTaken;
  r.electronicsTaken = eTaken;
  return r;
}

namespace {

struct RearmPlan {
  int buyP{0};
  int buyS{0};
  double useW{0.0};
  double useE{0.0};
  double cost{0.0};
  double util{0.0};
};

static RearmPlan bestRearmPlan(const econ::StationEconomyState& stEcon,
                              const econ::StationEconomyModel& model,
                              ShipHullClass hullClass,
                              WeaponType wP,
                              int needP,
                              WeaponType wS,
                              int needS,
                              const StationServicePriceModel& pm,
                              double creditsBudgetCr) {
  RearmPlan best{};

  needP = std::max(0, needP);
  needS = std::max(0, needS);

  const MissileRecipe recP = missileRecipeForWeapon(wP);
  const MissileRecipe recS = missileRecipeForWeapon(wS);

  const auto qW = q(stEcon, model, econ::CommodityId::Weapons, pm);
  const auto qE = q(stEcon, model, econ::CommodityId::Electronics, pm);

  const double invW = std::max(0.0, qW.inventory);
  const double invE = std::max(0.0, qE.inventory);

  const double fee = clamp01Safe(pm.feeRateEff);
  const double budget = budgetOrInf(creditsBudgetCr);

  auto consider = [&](int bp, int bs) {
    if (bp < 0 || bs < 0) return;
    if (bp > needP || bs > needS) return;

    const double useW = bp * recP.weapons + bs * recS.weapons;
    const double useE = bp * recP.electronics + bs * recS.electronics;

    if (useW > invW + 1e-9) return;
    if (useE > invE + 1e-9) return;

    const double cost = costForCommodities(useW, useE, 0.0, 0.0, qW.ask, qE.ask, 0.0, 0.0, fee);
    if (cost > budget + 1e-9) return;

    const double util = bp * recP.utility + bs * recS.utility;
    const int total = bp + bs;
    const int bestTotal = best.buyP + best.buyS;

    if (util > best.util + 1e-9 ||
        (std::abs(util - best.util) <= 1e-9 && total > bestTotal) ||
        (std::abs(util - best.util) <= 1e-9 && total == bestTotal && cost + 1e-9 < best.cost) ||
        (std::abs(util - best.util) <= 1e-9 && total == bestTotal && std::abs(cost - best.cost) <= 1e-9 && bp > best.buyP)) {
      best.buyP = bp;
      best.buyS = bs;
      best.useW = useW;
      best.useE = useE;
      best.cost = cost;
      best.util = util;
    }
  };

  for (int bs = 0; bs <= needS; ++bs) {
    for (int bp = 0; bp <= needP; ++bp) {
      consider(bp, bs);
    }
  }

  return best;
}

} // namespace

OrdnanceRearmAllQuote quoteOrdnanceRearmAllToFull(const econ::StationEconomyState& stEcon,
                                                 const econ::StationEconomyModel& model,
                                                 ShipHullClass hullClass,
                                                 WeaponType weaponPrimary,
                                                 core::u8 ammoPrimary,
                                                 WeaponType weaponSecondary,
                                                 core::u8 ammoSecondary,
                                                 const StationServicePriceModel& pm,
                                                 double creditsBudgetCr) {
  OrdnanceRearmAllQuote out{};
  out.weaponPrimary = weaponPrimary;
  out.weaponSecondary = weaponSecondary;

  const int maxP = weaponAmmoMax(weaponPrimary, hullClass);
  const int maxS = weaponAmmoMax(weaponSecondary, hullClass);

  const int haveP = (maxP > 0) ? clampInt((int)ammoPrimary, 0, maxP) : 0;
  const int haveS = (maxS > 0) ? clampInt((int)ammoSecondary, 0, maxS) : 0;

  const int needP = (maxP > 0) ? (maxP - haveP) : 0;
  const int needS = (maxS > 0) ? (maxS - haveS) : 0;

  if (needP + needS <= 0) {
    out.ok = false;
    out.reason = "no_need";
    return out;
  }

  const RearmPlan plan = bestRearmPlan(stEcon, model, hullClass,
                                      weaponPrimary, needP,
                                      weaponSecondary, needS,
                                      pm, creditsBudgetCr);

  out.buyPrimary = plan.buyP;
  out.buySecondary = plan.buyS;
  out.useWeapons = plan.useW;
  out.useElectronics = plan.useE;
  out.costCr = plan.cost;

  // Limit flags: compare with stock-only plan.
  const RearmPlan stockOnly = bestRearmPlan(stEcon, model, hullClass,
                                           weaponPrimary, needP,
                                           weaponSecondary, needS,
                                           pm, -1.0);

  out.limitedByStock = (stockOnly.buyP < needP) || (stockOnly.buyS < needS);

  const double budget = budgetOrInf(creditsBudgetCr);
  if (std::isfinite(budget) && budget < std::numeric_limits<double>::infinity()) {
    out.limitedByCredits = (plan.buyP < stockOnly.buyP) || (plan.buyS < stockOnly.buyS);
  } else {
    out.limitedByCredits = false;
  }

  if (out.buyPrimary + out.buySecondary <= 0) {
    if (stockOnly.buyP + stockOnly.buyS <= 0) {
      out.ok = false;
      out.reason = "out_of_stock";
    } else {
      out.ok = false;
      out.reason = "need_credits";
    }
    return out;
  }

  out.ok = true;
  return out;
}

OrdnanceRearmAllResult applyOrdnanceRearmAllToFull(econ::StationEconomyState& stEcon,
                                                  const econ::StationEconomyModel& model,
                                                  double& ioCredits,
                                                  ShipHullClass hullClass,
                                                  WeaponType weaponPrimary,
                                                  core::u8& ioAmmoPrimary,
                                                  WeaponType weaponSecondary,
                                                  core::u8& ioAmmoSecondary,
                                                  const StationServicePriceModel& pm) {
  OrdnanceRearmAllResult r{};

  const auto qAll = quoteOrdnanceRearmAllToFull(stEcon, model, hullClass,
                                               weaponPrimary, ioAmmoPrimary,
                                               weaponSecondary, ioAmmoSecondary,
                                               pm, ioCredits);
  if (!qAll.ok) {
    r.ok = false;
    r.reason = qAll.reason;
    return r;
  }

  if (ioCredits + 1e-9 < qAll.costCr) {
    r.ok = false;
    r.reason = "need_credits";
    return r;
  }

  const double wTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Weapons, qAll.useWeapons);
  const double eTaken = econ::takeInventory(stEcon, model, econ::CommodityId::Electronics, qAll.useElectronics);

  ioCredits -= qAll.costCr;

  const int maxP = weaponAmmoMax(weaponPrimary, hullClass);
  const int maxS = weaponAmmoMax(weaponSecondary, hullClass);

  if (maxP > 0) {
    const int haveP = clampInt((int)ioAmmoPrimary, 0, maxP);
    ioAmmoPrimary = (core::u8)clampInt(haveP + qAll.buyPrimary, 0, maxP);
  }

  if (maxS > 0) {
    const int haveS = clampInt((int)ioAmmoSecondary, 0, maxS);
    ioAmmoSecondary = (core::u8)clampInt(haveS + qAll.buySecondary, 0, maxS);
  }

  r.ok = true;
  r.primaryBought = qAll.buyPrimary;
  r.secondaryBought = qAll.buySecondary;
  r.creditsPaid = qAll.costCr;
  r.weaponsTaken = wTaken;
  r.electronicsTaken = eTaken;
  return r;
}

} // namespace stellar::sim
