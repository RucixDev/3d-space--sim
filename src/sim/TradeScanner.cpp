#include "stellar/sim/TradeScanner.h"

#include "stellar/econ/RoutePlanner.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace stellar::sim {

static double clampFee(double feeRate) {
  if (!std::isfinite(feeRate)) return 0.0;
  // Fees in this prototype are expected to be small (0..25%), but be generous.
  return std::clamp(feeRate, 0.0, 0.95);
}

static double effectiveCargoKg(const TradeScanParams& p) {
  double cap = p.cargoCapacityKg;
  if (!std::isfinite(cap)) cap = 0.0;
  cap = std::max(0.0, cap);

  if (!p.useFreeHold) return cap;

  double used = p.cargoUsedKg;
  if (!std::isfinite(used)) used = 0.0;
  used = std::max(0.0, used);

  return std::max(0.0, cap - used);
}


static double effectiveCargoKg(const TradeManifestScanParams& p) {
  double cap = p.cargoCapacityKg;
  if (!std::isfinite(cap)) cap = 0.0;
  cap = std::max(0.0, cap);

  if (!p.useFreeHold) return cap;

  double used = p.cargoUsedKg;
  if (!std::isfinite(used)) used = 0.0;
  used = std::max(0.0, used);

  return std::max(0.0, cap - used);
}

static void pushManifestOpportunity(std::vector<TradeManifestOpportunity>& out,
                                    const TradeManifestOpportunity& t,
                                    const TradeManifestScanParams& params) {
  if (t.cargoFilledKg <= 1e-9) return;
  if (t.netProfitCr <= 1e-9) return;
  if (t.netProfitCr + 1e-9 < params.minNetProfit) return;
  if (t.lines.empty()) return;
  out.push_back(t);
}


static TradeFeeRateFn defaultFeeFn() {
  return [](const Station& st) { return st.feeRate; };
}

static void pushOpportunity(std::vector<TradeOpportunity>& out,
                            const TradeOpportunity& t,
                            const TradeScanParams& params) {
  if (t.unitsPossible <= 1e-9) return;
  if (t.netProfitPerUnit <= 1e-9) return;
  if (t.netProfitTotal + 1e-9 < params.minNetProfit) return;
  out.push_back(t);
}

std::vector<TradeOpportunity> scanTradeOpportunities(Universe& u,
                                                     const SystemStub& originStub,
                                                     const Station& originStation,
                                                     double timeDays,
                                                     const std::vector<SystemStub>& candidates,
                                                     const TradeScanParams& params,
                                                     TradeFeeRateFn feeRate) {
  std::vector<TradeOpportunity> out;

  if (!std::isfinite(timeDays)) return out;
  if (originStub.id == 0 || originStation.id == 0) return out;
  if (candidates.empty()) return out;

  const double capKg = effectiveCargoKg(params);
  if (capKg <= 0.0) return out;

  if (!feeRate) feeRate = defaultFeeFn();

  // Pre-compute origin economy and fee.
  auto& fromEcon = u.stationEconomy(originStation, timeDays);
  const double feeFrom = clampFee(feeRate(originStation));

  // When filtering by commodity, query all commodities so we don't accidentally
  // miss the requested one due to the route planner's maxResults cutoff.
  const std::size_t routeQueryLimit = params.commodityFilterEnabled
                                      ? econ::kCommodityCount
                                      : std::max<std::size_t>(1, params.perStationLimit);

  out.reserve(std::min<std::size_t>(params.maxResults * 4, candidates.size() * 2));

  for (const auto& stub : candidates) {
    if (stub.id == 0) continue;
    if (!params.includeSameSystem && stub.id == originStub.id) continue;

    const auto& sys = u.getSystem(stub.id, &stub);
    const double distLy = (stub.posLy - originStub.posLy).length();

    for (const auto& toSt : sys.stations) {
      if (toSt.id == 0) continue;
      if (stub.id == originStub.id && toSt.id == originStation.id) continue;

      auto& toEcon = u.stationEconomy(toSt, timeDays);
      const double feeTo = clampFee(feeRate(toSt));

      const auto routes = econ::bestRoutesForCargo(fromEcon,
                                                   originStation.economyModel,
                                                   toEcon,
                                                   toSt.economyModel,
                                                   capKg,
                                                   feeFrom,
                                                   feeTo,
                                                   params.bidAskSpread,
                                                   routeQueryLimit);

      if (routes.empty()) continue;

      if (params.commodityFilterEnabled) {
        for (const auto& r : routes) {
          if (r.commodity != params.commodityFilter) continue;

          TradeOpportunity t;
          t.toSystem = sys.stub.id;
          t.toStation = toSt.id;
          t.toSystemName = sys.stub.name;
          t.toStationName = toSt.name;
          t.commodity = r.commodity;
          t.buyPrice = r.buyPrice;
          t.sellPrice = r.sellPrice;
          t.unitsFrom = r.unitsFrom;
          t.unitsToSpace = r.unitsToSpace;
          t.unitsPossible = r.unitsPossible;
          t.unitMassKg = r.unitMassKg;
          t.feeFrom = r.feeFrom;
          t.feeTo = r.feeTo;
          t.netProfitPerUnit = r.netProfitPerUnit;
          t.netProfitTotal = r.netProfitTotal;
          t.distanceLy = distLy;

          pushOpportunity(out, t, params);
          break; // only one idea for the requested commodity
        }
      } else {
        const std::size_t take = std::min<std::size_t>(params.perStationLimit, routes.size());
        for (std::size_t i = 0; i < take; ++i) {
          const auto& r = routes[i];

          TradeOpportunity t;
          t.toSystem = sys.stub.id;
          t.toStation = toSt.id;
          t.toSystemName = sys.stub.name;
          t.toStationName = toSt.name;
          t.commodity = r.commodity;
          t.buyPrice = r.buyPrice;
          t.sellPrice = r.sellPrice;
          t.unitsFrom = r.unitsFrom;
          t.unitsToSpace = r.unitsToSpace;
          t.unitsPossible = r.unitsPossible;
          t.unitMassKg = r.unitMassKg;
          t.feeFrom = r.feeFrom;
          t.feeTo = r.feeTo;
          t.netProfitPerUnit = r.netProfitPerUnit;
          t.netProfitTotal = r.netProfitTotal;
          t.distanceLy = distLy;

          pushOpportunity(out, t, params);
        }
      }
    }
  }

  std::sort(out.begin(), out.end(), [](const TradeOpportunity& a, const TradeOpportunity& b) {
    if (a.netProfitTotal != b.netProfitTotal) return a.netProfitTotal > b.netProfitTotal;
    if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
    if (a.toSystem != b.toSystem) return a.toSystem < b.toSystem;
    if (a.toStation != b.toStation) return a.toStation < b.toStation;
    return static_cast<std::size_t>(a.commodity) < static_cast<std::size_t>(b.commodity);
  });

  if (params.maxResults > 0 && out.size() > params.maxResults) {
    out.resize(params.maxResults);
  }

  return out;
}

std::vector<TradeOpportunity> scanTradeOpportunities(Universe& u,
                                                     const SystemStub& originStub,
                                                     const Station& originStation,
                                                     double timeDays,
                                                     const TradeScanParams& params,
                                                     TradeFeeRateFn feeRate) {
  const std::size_t maxSystems = std::max<std::size_t>(1, params.maxSystems);
  const auto candidates = u.queryNearby(originStub.posLy, params.radiusLy, maxSystems);
  return scanTradeOpportunities(u, originStub, originStation, timeDays, candidates, params, std::move(feeRate));
}



std::vector<TradeManifestOpportunity> scanTradeManifests(Universe& u,
                                                         const SystemStub& originStub,
                                                         const Station& originStation,
                                                         double timeDays,
                                                         const std::vector<SystemStub>& candidates,
                                                         const TradeManifestScanParams& params,
                                                         TradeFeeRateFn feeRate) {
  std::vector<TradeManifestOpportunity> out;

  if (!std::isfinite(timeDays)) return out;
  if (originStub.id == 0 || originStation.id == 0) return out;
  if (candidates.empty()) return out;

  const double capKg = effectiveCargoKg(params);
  if (capKg <= 0.0) return out;

  if (!feeRate) feeRate = defaultFeeFn();

  // Pre-compute origin economy and fee.
  auto& fromEcon = u.stationEconomy(originStation, timeDays);
  const double feeFrom = clampFee(feeRate(originStation));

  econ::CargoManifestParams mp;
  mp.cargoCapacityKg = capKg;
  mp.bidAskSpread = params.bidAskSpread;
  mp.fromFeeRate = feeFrom;
  mp.stepKg = params.stepKg;
  mp.maxBuyCreditsCr = params.maxBuyCreditsCr;
  mp.simulatePriceImpact = params.simulatePriceImpact;

  out.reserve(std::min<std::size_t>(params.maxResults * 3, candidates.size() * 2));

  for (const auto& stub : candidates) {
    if (stub.id == 0) continue;
    if (!params.includeSameSystem && stub.id == originStub.id) continue;

    const auto& sys = u.getSystem(stub.id, &stub);
    const double distLy = (stub.posLy - originStub.posLy).length();

    for (const auto& toSt : sys.stations) {
      if (toSt.id == 0) continue;
      if (stub.id == originStub.id && toSt.id == originStation.id) continue;

      auto& toEcon = u.stationEconomy(toSt, timeDays);
      const double feeTo = clampFee(feeRate(toSt));
      mp.toFeeRate = feeTo;

      const auto plan = econ::bestManifestForCargo(fromEcon,
                                                   originStation.economyModel,
                                                   toEcon,
                                                   toSt.economyModel,
                                                   mp);

      if (plan.lines.empty()) continue;

      TradeManifestOpportunity t;
      t.toSystem = sys.stub.id;
      t.toStation = toSt.id;
      t.toSystemName = sys.stub.name;
      t.toStationName = toSt.name;

      t.cargoFilledKg = plan.cargoFilledKg;
      t.netBuyCr = plan.netBuyCr;
      t.netSellCr = plan.netSellCr;
      t.netProfitCr = plan.netProfitCr;

      t.lines.reserve(plan.lines.size());
      for (const auto& l : plan.lines) {
        TradeManifestLine tl;
        tl.commodity = l.commodity;
        tl.units = l.units;
        tl.unitMassKg = l.unitMassKg;
        tl.massKg = l.massKg;
        tl.avgNetBuyPrice = l.avgNetBuyPrice;
        tl.avgNetSellPrice = l.avgNetSellPrice;
        tl.netBuyCr = l.netBuyCr;
        tl.netSellCr = l.netSellCr;
        tl.netProfitCr = l.netProfitCr;
        tl.netProfitPerUnit = l.netProfitPerUnit;
        tl.netProfitPerKg = l.netProfitPerKg;
        t.lines.push_back(std::move(tl));
      }

      t.distanceLy = distLy;

      pushManifestOpportunity(out, t, params);
    }
  }

  std::sort(out.begin(), out.end(), [](const TradeManifestOpportunity& a, const TradeManifestOpportunity& b) {
    if (a.netProfitCr != b.netProfitCr) return a.netProfitCr > b.netProfitCr;
    if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
    if (a.toSystem != b.toSystem) return a.toSystem < b.toSystem;
    return a.toStation < b.toStation;
  });

  if (params.maxResults > 0 && out.size() > params.maxResults) {
    out.resize(params.maxResults);
  }

  return out;
}

std::vector<TradeManifestOpportunity> scanTradeManifests(Universe& u,
                                                         const SystemStub& originStub,
                                                         const Station& originStation,
                                                         double timeDays,
                                                         const TradeManifestScanParams& params,
                                                         TradeFeeRateFn feeRate) {
  const std::size_t maxSystems = std::max<std::size_t>(1, params.maxSystems);
  const auto candidates = u.queryNearby(originStub.posLy, params.radiusLy, maxSystems);
  return scanTradeManifests(u, originStub, originStation, timeDays, candidates, params, std::move(feeRate));
}


} // namespace stellar::sim
