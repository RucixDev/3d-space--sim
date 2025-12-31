#include "stellar/econ/RoutePlanner.h"

#include "stellar/econ/Commodity.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace stellar::econ {

static constexpr std::size_t idx(CommodityId id) { return static_cast<std::size_t>(id); }

std::vector<RouteOpportunity> bestRoutes(const StationEconomyState& fromState,
                                        const StationEconomyModel& fromModel,
                                        const StationEconomyState& toState,
                                        const StationEconomyModel& toModel,
                                        double bidAskSpread,
                                        std::size_t maxResults) {
  std::vector<RouteOpportunity> out;
  out.reserve(kCommodityCount);

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    const CommodityId cid = static_cast<CommodityId>(i);
    const auto qFrom = quote(fromState, fromModel, cid, bidAskSpread);
    const auto qTo   = quote(toState, toModel, cid, bidAskSpread);

    const double profit = qTo.bid - qFrom.ask;

    // Feasibility: don't recommend routes where you cannot actually buy/sell.
    const double unitsFrom = std::max(0.0, qFrom.inventory);
    const double unitsToSpace = std::max(0.0, toModel.capacity[i] - toState.inventory[i]);
    const double unitsPossible = std::min(unitsFrom, unitsToSpace);

    if (profit > 0.0 && unitsPossible > 1e-9) {
      const double massKg = std::max(1e-6, commodityDef(cid).massKg);

      RouteOpportunity r{};
      r.commodity = cid;
      r.profitPerUnit = profit;
      r.buyPrice = qFrom.ask;
      r.sellPrice = qTo.bid;

      r.unitsFrom = unitsFrom;
      r.unitsToSpace = unitsToSpace;
      r.unitsPossible = unitsPossible;

      r.unitMassKg = massKg;
      r.profitPerKg = profit / massKg;
      r.profitTotal = profit * unitsPossible;

      // Defaults: no fees applied.
      r.netProfitPerUnit = r.profitPerUnit;
      r.netProfitTotal = r.profitTotal;

      out.push_back(r);
    }
  }

  std::sort(out.begin(), out.end(), [](const RouteOpportunity& a, const RouteOpportunity& b) {
    return a.profitPerUnit > b.profitPerUnit;
  });

  if (out.size() > maxResults) out.resize(maxResults);
  return out;
}

std::vector<RouteOpportunity> bestRoutesForCargo(const StationEconomyState& fromState,
                                                const StationEconomyModel& fromModel,
                                                const StationEconomyState& toState,
                                                const StationEconomyModel& toModel,
                                                double cargoCapacityKg,
                                                double fromFeeRate,
                                                double toFeeRate,
                                                double bidAskSpread,
                                                std::size_t maxResults) {
  std::vector<RouteOpportunity> out;
  out.reserve(kCommodityCount);

  cargoCapacityKg = std::max(0.0, cargoCapacityKg);
  fromFeeRate = std::clamp(fromFeeRate, 0.0, 1.0);
  toFeeRate = std::clamp(toFeeRate, 0.0, 1.0);

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    const CommodityId cid = static_cast<CommodityId>(i);
    const auto qFrom = quote(fromState, fromModel, cid, bidAskSpread);
    const auto qTo   = quote(toState, toModel, cid, bidAskSpread);

    const double rawProfit = qTo.bid - qFrom.ask;
    if (rawProfit <= 0.0) continue;

    // Station feasibility limits.
    const double unitsFrom = std::max(0.0, qFrom.inventory);
    const double unitsToSpace = std::max(0.0, toModel.capacity[i] - toState.inventory[i]);
    double unitsPossible = std::min(unitsFrom, unitsToSpace);
    if (unitsPossible <= 1e-9) continue;

    const double massKg = std::max(1e-6, commodityDef(cid).massKg);

    // Cargo feasibility limit.
    if (cargoCapacityKg > 0.0) {
      const double cargoUnits = std::floor(cargoCapacityKg / massKg + 1e-9);
      unitsPossible = std::min(unitsPossible, std::max(0.0, cargoUnits));
      if (unitsPossible <= 1e-9) continue;
    }

    // Net-of-fees profit.
    const double netBuy = qFrom.ask * (1.0 + fromFeeRate);
    const double netSell = qTo.bid * (1.0 - toFeeRate);
    const double netProfit = netSell - netBuy;
    if (netProfit <= 0.0) continue;

    RouteOpportunity r{};
    r.commodity = cid;
    r.profitPerUnit = rawProfit;
    r.buyPrice = qFrom.ask;
    r.sellPrice = qTo.bid;
    r.unitsFrom = unitsFrom;
    r.unitsToSpace = unitsToSpace;
    r.unitsPossible = unitsPossible;
    r.unitMassKg = massKg;
    r.profitPerKg = rawProfit / massKg;
    r.profitTotal = rawProfit * unitsPossible;

    r.feeFrom = fromFeeRate;
    r.feeTo = toFeeRate;
    r.netProfitPerUnit = netProfit;
    r.netProfitTotal = netProfit * unitsPossible;

    out.push_back(r);
  }

  std::sort(out.begin(), out.end(), [](const RouteOpportunity& a, const RouteOpportunity& b) {
    // Primary: trip profit.
    if (a.netProfitTotal != b.netProfitTotal) return a.netProfitTotal > b.netProfitTotal;
    // Secondary: per-unit profit.
    return a.netProfitPerUnit > b.netProfitPerUnit;
  });

  if (out.size() > maxResults) out.resize(maxResults);
  return out;
}



static double clamp01(double x) {
  if (!std::isfinite(x)) return 0.0;
  return std::clamp(x, 0.0, 1.0);
}

static double safePos(double x) {
  if (!std::isfinite(x)) return 0.0;
  return std::max(0.0, x);
}

static void clampEconomyState(StationEconomyState& s, const StationEconomyModel& m) {
  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    double cap = m.capacity[i];
    if (!std::isfinite(cap)) cap = 0.0;
    cap = std::max(0.0, cap);

    double inv = s.inventory[i];
    if (!std::isfinite(inv)) inv = 0.0;
    s.inventory[i] = std::clamp(inv, 0.0, cap);
  }
}

CargoManifestPlan bestManifestForCargo(const StationEconomyState& fromState,
                                       const StationEconomyModel& fromModel,
                                       const StationEconomyState& toState,
                                       const StationEconomyModel& toModel,
                                       const CargoManifestParams& params) {
  CargoManifestPlan plan{};
  plan.cargoCapacityKg = safePos(params.cargoCapacityKg);
  if (plan.cargoCapacityKg <= 1e-9) return plan;

  CargoManifestParams p = params;
  p.fromFeeRate = clamp01(p.fromFeeRate);
  p.toFeeRate = clamp01(p.toFeeRate);

  if (!std::isfinite(p.bidAskSpread)) p.bidAskSpread = 0.10;
  p.bidAskSpread = std::clamp(p.bidAskSpread, 0.0, 1.0);

  double stepKg = p.stepKg;
  if (!std::isfinite(stepKg)) stepKg = 1.0;
  // Be conservative: overly small steps can be expensive in tools.
  stepKg = std::clamp(stepKg, 0.05, std::max(0.05, plan.cargoCapacityKg));

  double maxBuy = p.maxBuyCreditsCr;
  if (!std::isfinite(maxBuy)) maxBuy = 0.0;
  maxBuy = std::max(0.0, maxBuy);
  const bool useCredits = maxBuy > 0.0;

  // Mutable copies are used for:
  //  - enforcing per-commodity availability (inventory/capacity)
  //  - optionally simulating price impact (when p.simulatePriceImpact is true)
  StationEconomyState from = fromState;
  StationEconomyState to = toState;
  clampEconomyState(from, fromModel);
  clampEconomyState(to, toModel);

  struct Accum {
    double units{0.0};
    double massKg{0.0};
    double buyCr{0.0};
    double sellCr{0.0};
  };
  std::array<Accum, kCommodityCount> acc{};

  double filledKg = 0.0;
  double totalBuy = 0.0;
  double totalSell = 0.0;

  // Hard guard: prevents accidental infinite loops when stepKg is tiny.
  const std::size_t maxIter = static_cast<std::size_t>(std::ceil(plan.cargoCapacityKg / stepKg + 1e-6)) + 64;

  for (std::size_t iter = 0; iter < maxIter; ++iter) {
    const double remainingKg = std::max(0.0, plan.cargoCapacityKg - filledKg);
    if (remainingKg <= 1e-9) break;
    if (useCredits && totalBuy + 1e-9 >= maxBuy) break;

    const double thisStepKg = std::min(stepKg, remainingKg);

    int bestIdx = -1;
    double bestProfitPerKg = 0.0;
    double bestUnits = 0.0;
    double bestNetBuy = 0.0;
    double bestNetSell = 0.0;

    const StationEconomyState& priceFrom = p.simulatePriceImpact ? from : fromState;
    const StationEconomyState& priceTo = p.simulatePriceImpact ? to : toState;

    for (std::size_t i = 0; i < kCommodityCount; ++i) {
      const CommodityId cid = static_cast<CommodityId>(i);

      const double unitMass = std::max(1e-6, commodityDef(cid).massKg);

      // Availability constraints (always tracked using the mutable copies).
      const double unitsFrom = safePos(from.inventory[i]);
      const double capTo = safePos(toModel.capacity[i]);
      const double unitsToSpace = std::max(0.0, capTo - safePos(to.inventory[i]));
      const double unitsStation = std::min(unitsFrom, unitsToSpace);
      if (unitsStation <= 1e-9) continue;

      // Mass constraint for this step.
      double deltaUnits = std::min(unitsStation, thisStepKg / unitMass);
      if (deltaUnits <= 1e-9) continue;

      // Price / profit at the current step.
      const auto qFrom = quote(priceFrom, fromModel, cid, p.bidAskSpread);
      const auto qTo   = quote(priceTo, toModel, cid, p.bidAskSpread);

      const double netBuyPerUnit = qFrom.ask * (1.0 + p.fromFeeRate);
      const double netSellPerUnit = qTo.bid * (1.0 - p.toFeeRate);
      const double netProfitPerUnit = netSellPerUnit - netBuyPerUnit;
      if (netProfitPerUnit <= 1e-9) continue;

      if (useCredits) {
        const double remainingCr = std::max(0.0, maxBuy - totalBuy);
        if (remainingCr <= 1e-9) continue;
        const double affordable = remainingCr / std::max(1e-9, netBuyPerUnit);
        deltaUnits = std::min(deltaUnits, affordable);
        if (deltaUnits <= 1e-9) continue;
      }

      const double profitPerKg = netProfitPerUnit / unitMass;

      // Greedy selection: maximize marginal profit per kg.
      if (bestIdx < 0 || profitPerKg > bestProfitPerKg + 1e-12) {
        bestIdx = (int)i;
        bestProfitPerKg = profitPerKg;
        bestUnits = deltaUnits;
        bestNetBuy = netBuyPerUnit;
        bestNetSell = netSellPerUnit;
      } else if (std::abs(profitPerKg - bestProfitPerKg) <= 1e-12) {
        // Tie-breaker: prefer larger absolute per-unit profit.
        if (netProfitPerUnit > (bestNetSell - bestNetBuy) + 1e-9) {
          bestIdx = (int)i;
          bestProfitPerKg = profitPerKg;
          bestUnits = deltaUnits;
          bestNetBuy = netBuyPerUnit;
          bestNetSell = netSellPerUnit;
        }
      }
    }

    if (bestIdx < 0 || bestUnits <= 1e-9) break;

    const CommodityId bestCid = static_cast<CommodityId>(bestIdx);
    const double unitMass = std::max(1e-6, commodityDef(bestCid).massKg);
    const double massUsed = bestUnits * unitMass;
    if (massUsed <= 1e-9) break;

    // Advance "allocation state": always update copies so we don't exceed station limits.
    from.inventory[(std::size_t)bestIdx] = std::max(0.0, from.inventory[(std::size_t)bestIdx] - bestUnits);
    const double capTo = safePos(toModel.capacity[(std::size_t)bestIdx]);
    to.inventory[(std::size_t)bestIdx] = std::min(capTo, safePos(to.inventory[(std::size_t)bestIdx]) + bestUnits);

    const double buyCr = bestNetBuy * bestUnits;
    const double sellCr = bestNetSell * bestUnits;

    acc[(std::size_t)bestIdx].units += bestUnits;
    acc[(std::size_t)bestIdx].massKg += massUsed;
    acc[(std::size_t)bestIdx].buyCr += buyCr;
    acc[(std::size_t)bestIdx].sellCr += sellCr;

    filledKg += massUsed;
    totalBuy += buyCr;
    totalSell += sellCr;

    // Defensive: if we don't make progress, stop.
    if (massUsed < 1e-9) break;
  }

  plan.cargoFilledKg = filledKg;
  plan.netBuyCr = totalBuy;
  plan.netSellCr = totalSell;
  plan.netProfitCr = totalSell - totalBuy;

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    if (acc[i].units <= 1e-9) continue;

    const CommodityId cid = static_cast<CommodityId>(i);
    CargoManifestLine line{};
    line.commodity = cid;
    line.units = acc[i].units;
    line.unitMassKg = std::max(1e-6, commodityDef(cid).massKg);
    line.massKg = acc[i].massKg;

    line.netBuyCr = acc[i].buyCr;
    line.netSellCr = acc[i].sellCr;
    line.netProfitCr = line.netSellCr - line.netBuyCr;

    line.avgNetBuyPrice = (line.units > 1e-9) ? (line.netBuyCr / line.units) : 0.0;
    line.avgNetSellPrice = (line.units > 1e-9) ? (line.netSellCr / line.units) : 0.0;
    line.netProfitPerUnit = (line.units > 1e-9) ? (line.netProfitCr / line.units) : 0.0;
    line.netProfitPerKg = (line.massKg > 1e-9) ? (line.netProfitCr / line.massKg) : 0.0;

    plan.lines.push_back(std::move(line));
  }

  std::sort(plan.lines.begin(), plan.lines.end(), [](const CargoManifestLine& a, const CargoManifestLine& b) {
    if (a.netProfitCr != b.netProfitCr) return a.netProfitCr > b.netProfitCr;
    return a.netProfitPerKg > b.netProfitPerKg;
  });

  // If we ended up with negative profit (possible if numerical noise / tiny steps),
  // return an empty plan.
  if (plan.netProfitCr <= 1e-9) {
    plan = CargoManifestPlan{};
    plan.cargoCapacityKg = safePos(params.cargoCapacityKg);
  }

  return plan;
}


} // namespace stellar::econ
