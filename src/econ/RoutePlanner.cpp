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

struct ManifestAccum {
  double units{0.0};
  double massKg{0.0};
  double buyCr{0.0};
  double sellCr{0.0};
};

static CargoManifestPlan finalizeManifestPlan(double cargoCapacityKg,
                                             double filledKg,
                                             double totalBuy,
                                             double totalSell,
                                             const std::array<ManifestAccum, kCommodityCount>& acc,
                                             double originalCapacityKg) {
  CargoManifestPlan plan{};
  plan.cargoCapacityKg = cargoCapacityKg;
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
    plan.cargoCapacityKg = originalCapacityKg;
  }

  return plan;
}

// Lightweight quote helper for the manifest planner.
//
// We intentionally avoid copying StationEconomyState history vectors into the planner state.
// Only inventory levels affect prices, so we operate on inventory arrays directly.
struct QuoteLite {
  double mid{0.0};
  double ask{0.0};
  double bid{0.0};
};

static QuoteLite quoteLite(const std::array<double, kCommodityCount>& inv,
                           const StationEconomyModel& model,
                           CommodityId id,
                           double bidAskSpread) {
  QuoteLite q{};

  if (!std::isfinite(bidAskSpread)) bidAskSpread = 0.0;
  bidAskSpread = std::clamp(bidAskSpread, 0.0, 1.0);
  const double half = bidAskSpread * 0.5;

  const std::size_t i = idx(id);
  const auto& def = commodityDef(id);

  // Clamp inventory defensively.
  double cap = model.capacity[i];
  if (!std::isfinite(cap)) cap = 0.0;
  cap = std::max(0.0, cap);

  double cur = inv[i];
  if (!std::isfinite(cur)) cur = 0.0;
  cur = std::clamp(cur, 0.0, cap);

  // Match econ::midPrice(...) exactly (but without requiring StationEconomyState).
  const double desired = model.desiredStock[i];
  double mid = def.basePrice;

  if (desired > 1e-9) {
    const double ratio = cur / desired; // 1.0 means "normal"
    double factor = 1.0 + model.priceVolatility * (1.0 - ratio);

    // Scarcity spike if nearly empty
    if (cur < desired * 0.05) factor *= 1.4;

    factor = std::clamp(factor, 0.2, 5.0);
    mid = def.basePrice * factor;
  }

  q.mid = mid;
  q.ask = mid * (1.0 + half);
  q.bid = mid * (1.0 - half);
  return q;
}

static CargoManifestPlan bestManifestForCargoGreedy(const StationEconomyState& fromState,
                                                   const StationEconomyModel& fromModel,
                                                   const StationEconomyState& toState,
                                                   const StationEconomyModel& toModel,
                                                   const CargoManifestParams& params,
                                                   double stepKg,
                                                   double maxBuy,
                                                   bool useCredits) {
  // Mutable copies are used for:
  //  - enforcing per-commodity availability (inventory/capacity)
  //  - optionally simulating price impact (when params.simulatePriceImpact is true)
  StationEconomyState from{};
  StationEconomyState to{};
  from.inventory = fromState.inventory;
  to.inventory = toState.inventory;
  clampEconomyState(from, fromModel);
  clampEconomyState(to, toModel);

  // Base pricing snapshot (used when simulatePriceImpact==false).
  const StationEconomyState baseFrom = from;
  const StationEconomyState baseTo = to;

  std::array<ManifestAccum, kCommodityCount> acc{};

  double filledKg = 0.0;
  double totalBuy = 0.0;
  double totalSell = 0.0;

  // Hard guard: prevents accidental infinite loops when stepKg is tiny.
  const std::size_t maxIter = static_cast<std::size_t>(std::ceil(params.cargoCapacityKg / stepKg + 1e-6)) + 64;

  for (std::size_t iter = 0; iter < maxIter; ++iter) {
    const double remainingKg = std::max(0.0, params.cargoCapacityKg - filledKg);
    if (remainingKg <= 1e-9) break;
    if (useCredits && totalBuy + 1e-9 >= maxBuy) break;

    const double thisStepKg = std::min(stepKg, remainingKg);

    int bestIdx = -1;
    double bestProfitPerKg = 0.0;
    double bestUnits = 0.0;
    double bestNetBuy = 0.0;
    double bestNetSell = 0.0;

    const StationEconomyState& priceFrom = params.simulatePriceImpact ? from : baseFrom;
    const StationEconomyState& priceTo = params.simulatePriceImpact ? to : baseTo;

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
      const auto qFrom = quote(priceFrom, fromModel, cid, params.bidAskSpread);
      const auto qTo   = quote(priceTo, toModel, cid, params.bidAskSpread);

      const double netBuyPerUnit = qFrom.ask * (1.0 + params.fromFeeRate);
      const double netSellPerUnit = qTo.bid * (1.0 - params.toFeeRate);
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

  return finalizeManifestPlan(params.cargoCapacityKg,
                              filledKg,
                              totalBuy,
                              totalSell,
                              acc,
                              safePos(params.cargoCapacityKg));
}

static CargoManifestPlan bestManifestForCargoBeamSearch(const StationEconomyState& fromState,
                                                       const StationEconomyModel& fromModel,
                                                       const StationEconomyState& toState,
                                                       const StationEconomyModel& toModel,
                                                       const CargoManifestParams& params,
                                                       double stepKg,
                                                       double maxBuy,
                                                       bool useCredits) {
  // Clamp beam width aggressively: this planner may be used in UI scans.
  std::size_t beamWidth = params.beamWidth;
  if (beamWidth == 0) beamWidth = 1;
  beamWidth = std::clamp<std::size_t>(beamWidth, 1, 128);

  // For safety/perf: cap iterations by inflating stepKg if needed.
  std::size_t maxIter = static_cast<std::size_t>(std::ceil(params.cargoCapacityKg / stepKg + 1e-6)) + 64;
  constexpr std::size_t kHardMaxIter = 1024;
  if (maxIter > kHardMaxIter) {
    const double minStep = params.cargoCapacityKg / std::max<double>(1.0, (double)(kHardMaxIter - 64));
    stepKg = std::max(stepKg, minStep);
    maxIter = static_cast<std::size_t>(std::ceil(params.cargoCapacityKg / stepKg + 1e-6)) + 64;
  }

  // Clamp inventories (without copying history vectors).
  StationEconomyState baseFrom{};
  StationEconomyState baseTo{};
  baseFrom.inventory = fromState.inventory;
  baseTo.inventory = toState.inventory;
  clampEconomyState(baseFrom, fromModel);
  clampEconomyState(baseTo, toModel);

  const std::array<double, kCommodityCount> baseFromInv = baseFrom.inventory;
  const std::array<double, kCommodityCount> baseToInv = baseTo.inventory;

  // Precompute optimistic bounds (used for beam ranking).
  //
  // We rank partial states by: score = currentProfit + optimisticRemainingProfit,
  // where optimisticRemainingProfit is bounded by both:
  //  - remaining mass (bestProfitPerKgBound), and
  //  - remaining credits (bestProfitPerCrBound) when a credit cap is active.
  double bestProfitPerKgBound = 0.0;
  double bestProfitPerCrBound = 0.0;

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    const CommodityId cid = static_cast<CommodityId>(i);
    const double unitMass = std::max(1e-6, commodityDef(cid).massKg);

    const auto qFrom = quoteLite(baseFromInv, fromModel, cid, params.bidAskSpread);
    const auto qTo = quoteLite(baseToInv, toModel, cid, params.bidAskSpread);

    const double netBuy = qFrom.ask * (1.0 + params.fromFeeRate);
    const double netSell = qTo.bid * (1.0 - params.toFeeRate);
    const double netProfit = netSell - netBuy;
    if (netProfit <= 1e-12) continue;

    bestProfitPerKgBound = std::max(bestProfitPerKgBound, netProfit / unitMass);

    if (netBuy > 1e-12) {
      bestProfitPerCrBound = std::max(bestProfitPerCrBound, netProfit / netBuy);
    }
  }

  struct BeamState {
    std::array<double, kCommodityCount> fromInv{};
    std::array<double, kCommodityCount> toInv{};
    std::array<ManifestAccum, kCommodityCount> acc{};
    double filledKg{0.0};
    double buyCr{0.0};
    double sellCr{0.0};
    std::uint64_t seq{0};
  };

  auto profitCr = [](const BeamState& s) { return s.sellCr - s.buyCr; };

  auto score = [&](const BeamState& s) {
    const double curProfit = profitCr(s);
    const double remainingKg = std::max(0.0, params.cargoCapacityKg - s.filledKg);

    double boundMass = remainingKg * bestProfitPerKgBound;
    double boundCr = boundMass;

    if (useCredits) {
      const double remainingCr = std::max(0.0, maxBuy - s.buyCr);
      boundCr = remainingCr * bestProfitPerCrBound;
    }

    const double optimistic = std::min(boundMass, boundCr);
    return curProfit + optimistic;
  };

  // Initial beam.
  std::vector<BeamState> beam;
  beam.reserve(beamWidth);

  BeamState root{};
  root.fromInv = baseFromInv;
  root.toInv = baseToInv;
  beam.push_back(std::move(root));

  std::uint64_t seqCounter = 1;

  for (std::size_t iter = 0; iter < maxIter; ++iter) {
    std::vector<BeamState> next;
    next.reserve(beam.size() * (kCommodityCount + 1));

    bool anyExpanded = false;

    for (const auto& st : beam) {
      const double remainingKg = std::max(0.0, params.cargoCapacityKg - st.filledKg);
      const bool creditsExhausted = useCredits && (st.buyCr + 1e-9 >= maxBuy);

      if (remainingKg <= 1e-9 || creditsExhausted) {
        // Terminal: keep it around in case it's the best overall.
        next.push_back(st);
        continue;
      }

      const double thisStepKg = std::min(stepKg, remainingKg);

      for (std::size_t i = 0; i < kCommodityCount; ++i) {
        const CommodityId cid = static_cast<CommodityId>(i);
        const double unitMass = std::max(1e-6, commodityDef(cid).massKg);

        const double unitsFrom = safePos(st.fromInv[i]);
        const double capTo = safePos(toModel.capacity[i]);
        const double unitsToSpace = std::max(0.0, capTo - safePos(st.toInv[i]));
        const double unitsStation = std::min(unitsFrom, unitsToSpace);
        if (unitsStation <= 1e-9) continue;

        double deltaUnits = std::min(unitsStation, thisStepKg / unitMass);
        if (deltaUnits <= 1e-9) continue;

        // Price snapshot: either dynamic (price impact) or fixed (baseline).
        const auto& priceFromInv = params.simulatePriceImpact ? st.fromInv : baseFromInv;
        const auto& priceToInv = params.simulatePriceImpact ? st.toInv : baseToInv;

        const auto qFrom = quoteLite(priceFromInv, fromModel, cid, params.bidAskSpread);
        const auto qTo = quoteLite(priceToInv, toModel, cid, params.bidAskSpread);

        const double netBuyPerUnit = qFrom.ask * (1.0 + params.fromFeeRate);
        const double netSellPerUnit = qTo.bid * (1.0 - params.toFeeRate);
        const double netProfitPerUnit = netSellPerUnit - netBuyPerUnit;
        if (netProfitPerUnit <= 1e-9) continue;

        if (useCredits) {
          const double remainingCr = std::max(0.0, maxBuy - st.buyCr);
          if (remainingCr <= 1e-9) continue;
          const double affordable = remainingCr / std::max(1e-9, netBuyPerUnit);
          deltaUnits = std::min(deltaUnits, affordable);
          if (deltaUnits <= 1e-9) continue;
        }

        const double massUsed = deltaUnits * unitMass;
        if (massUsed <= 1e-9) continue;

        BeamState ns = st;
        ns.seq = seqCounter++;

        ns.fromInv[i] = std::max(0.0, ns.fromInv[i] - deltaUnits);
        ns.toInv[i] = std::min(capTo, std::max(0.0, ns.toInv[i] + deltaUnits));

        ns.acc[i].units += deltaUnits;
        ns.acc[i].massKg += massUsed;
        ns.acc[i].buyCr += netBuyPerUnit * deltaUnits;
        ns.acc[i].sellCr += netSellPerUnit * deltaUnits;

        ns.filledKg += massUsed;
        ns.buyCr += netBuyPerUnit * deltaUnits;
        ns.sellCr += netSellPerUnit * deltaUnits;

        next.push_back(std::move(ns));
        anyExpanded = true;
      }
    }

    if (!anyExpanded) break;

    std::sort(next.begin(), next.end(), [&](const BeamState& a, const BeamState& b) {
      const double sa = score(a);
      const double sb = score(b);
      if (sa != sb) return sa > sb;

      const double pa = profitCr(a);
      const double pb = profitCr(b);
      if (pa != pb) return pa > pb;

      if (a.filledKg != b.filledKg) return a.filledKg > b.filledKg;

      // Deterministic final tie-break.
      return a.seq < b.seq;
    });

    if (next.size() > beamWidth) next.resize(beamWidth);
    beam = std::move(next);
  }

  // Pick the best by realized profit (not heuristic score).
  const BeamState* best = nullptr;
  for (const auto& st : beam) {
    if (!best) {
      best = &st;
      continue;
    }
    const double pA = profitCr(st);
    const double pB = profitCr(*best);
    if (pA > pB + 1e-12) {
      best = &st;
    } else if (std::abs(pA - pB) <= 1e-12) {
      if (st.filledKg > best->filledKg + 1e-9) best = &st;
      else if (std::abs(st.filledKg - best->filledKg) <= 1e-9 && st.seq < best->seq) best = &st;
    }
  }

  if (!best) {
    CargoManifestPlan empty{};
    empty.cargoCapacityKg = safePos(params.cargoCapacityKg);
    return empty;
  }

  return finalizeManifestPlan(params.cargoCapacityKg,
                              best->filledKg,
                              best->buyCr,
                              best->sellCr,
                              best->acc,
                              safePos(params.cargoCapacityKg));
}

CargoManifestPlan bestManifestForCargo(const StationEconomyState& fromState,
                                       const StationEconomyModel& fromModel,
                                       const StationEconomyState& toState,
                                       const StationEconomyModel& toModel,
                                       const CargoManifestParams& params) {
  const double cap = safePos(params.cargoCapacityKg);

  CargoManifestPlan plan{};
  plan.cargoCapacityKg = cap;
  if (cap <= 1e-9) return plan;

  CargoManifestParams p = params;
  p.cargoCapacityKg = cap;
  p.fromFeeRate = clamp01(p.fromFeeRate);
  p.toFeeRate = clamp01(p.toFeeRate);

  if (!std::isfinite(p.bidAskSpread)) p.bidAskSpread = 0.10;
  p.bidAskSpread = std::clamp(p.bidAskSpread, 0.0, 1.0);

  double stepKg = p.stepKg;
  if (!std::isfinite(stepKg)) stepKg = 1.0;
  // Be conservative: overly small steps can be expensive in tools.
  stepKg = std::clamp(stepKg, 0.05, std::max(0.05, cap));

  double maxBuy = p.maxBuyCreditsCr;
  if (!std::isfinite(maxBuy)) maxBuy = 0.0;
  maxBuy = std::max(0.0, maxBuy);
  const bool useCredits = maxBuy > 0.0;

  // Default: existing greedy behavior.
  if (p.planner == CargoManifestPlanner::BeamSearch) {
    return bestManifestForCargoBeamSearch(fromState, fromModel, toState, toModel, p, stepKg, maxBuy, useCredits);
  }

  return bestManifestForCargoGreedy(fromState, fromModel, toState, toModel, p, stepKg, maxBuy, useCredits);
}



} // namespace stellar::econ
