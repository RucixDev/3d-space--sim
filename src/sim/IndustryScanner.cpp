#include "stellar/sim/IndustryScanner.h"

#include "stellar/econ/Market.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace stellar::sim {

static double clampFee(double feeRate) {
  if (!std::isfinite(feeRate)) return 0.0;
  return std::clamp(feeRate, 0.0, 0.95);
}

static double effectiveCargoKg(const IndustryTradeScanParams& p) {
  double cap = p.cargoCapacityKg;
  if (!std::isfinite(cap)) cap = 0.0;
  cap = std::max(0.0, cap);

  if (!p.useFreeHold) return cap;

  double used = p.cargoUsedKg;
  if (!std::isfinite(used)) used = 0.0;
  used = std::max(0.0, used);
  return std::max(0.0, cap - used);
}

static IndustryFeeRateFn defaultFeeFn() {
  return [](const Station& st) { return st.feeRate; };
}

static void pushOpportunity(std::vector<IndustryTradeOpportunity>& out,
                            const IndustryTradeOpportunity& t,
                            const IndustryTradeScanParams& params) {
  if (t.batches <= 1e-9) return;
  if (t.outputUnits <= 1e-9) return;
  if (t.netProfitCr <= 1e-9) return;
  if (t.netProfitCr + 1e-9 < params.minNetProfit) return;
  out.push_back(t);
}

// Compute max batches allowed by origin input inventory.
static int maxBatchesByInputInventory(econ::StationEconomyState& fromEcon,
                                      const econ::StationEconomyModel& fromModel,
                                      const IndustryQuote& q1,
                                      double bidAskSpread) {
  int maxA = 1'000'000;
  int maxB = 1'000'000;

  if (q1.inputAUnits > 1e-9) {
    const auto qa = econ::quote(fromEcon, fromModel, q1.inputA, bidAskSpread);
    maxA = (int)std::floor(std::max(0.0, qa.inventory) / q1.inputAUnits + 1e-9);
  }

  if (q1.inputBUnits > 1e-9) {
    const auto qb = econ::quote(fromEcon, fromModel, q1.inputB, bidAskSpread);
    maxB = (int)std::floor(std::max(0.0, qb.inventory) / q1.inputBUnits + 1e-9);
  }

  return std::max(0, std::min(maxA, maxB));
}

static int maxBatchesByOutputCargo(double capKg, const IndustryQuote& q1) {
  if (capKg <= 1e-9) return 0;
  if (q1.outputUnits <= 1e-9) return 0;

  const double unitMassKg = econ::commodityDef(q1.output).massKg;
  const double massPerBatch = q1.outputUnits * unitMassKg;
  if (massPerBatch <= 1e-9) {
    // Massless commodity (unlikely): allow a lot.
    return 1'000'000;
  }

  return (int)std::floor(capKg / massPerBatch + 1e-9);
}

static int maxBatchesByDestinationSpace(econ::StationEconomyState& toEcon,
                                        const econ::StationEconomyModel& toModel,
                                        econ::CommodityId output,
                                        double outUnitsPerBatch,
                                        double bidAskSpread) {
  if (outUnitsPerBatch <= 1e-9) return 0;

  const auto q = econ::quote(toEcon, toModel, output, bidAskSpread);
  const std::size_t idx = static_cast<std::size_t>(output);
  const double cap = std::max(0.0, toModel.capacity[idx]);
  const double cur = std::clamp(q.inventory, 0.0, cap);
  const double space = std::max(0.0, cap - cur);

  return (int)std::floor(space / outUnitsPerBatch + 1e-9);
}

std::vector<IndustryTradeOpportunity> scanIndustryTradeOpportunities(Universe& u,
                                                                     const SystemStub& originStub,
                                                                     const Station& originStation,
                                                                     double timeDays,
                                                                     const std::vector<SystemStub>& candidates,
                                                                     const IndustryTradeScanParams& params,
                                                                     IndustryFeeRateFn feeRate) {
  std::vector<IndustryTradeOpportunity> out;

  if (!std::isfinite(timeDays)) return out;
  if (originStub.id == 0 || originStation.id == 0) return out;
  if (candidates.empty()) return out;

  const double capKg = effectiveCargoKg(params);
  if (capKg <= 0.0) return out;

  if (!feeRate) feeRate = defaultFeeFn();

  // Origin station: economy + fee.
  auto& fromEcon = u.stationEconomy(originStation, timeDays);
  const double feeFrom = clampFee(feeRate(originStation));

  const auto recipes = availableIndustryRecipes(originStation.type);
  if (recipes.empty()) return out;

  // Precompute per-recipe caps that are independent of destination.
  struct RecipeCaps {
    const IndustryRecipeDef* def{nullptr};
    IndustryQuote q1{}; // 1 batch quote (includes yield/speed/fee mods)
    int maxByInputs{0};
    int maxByCargo{0};
    int maxOverall{0};
  };

  std::vector<RecipeCaps> caps;
  caps.reserve(recipes.size());

  for (const auto* r : recipes) {
    if (!r) continue;

    RecipeCaps rc;
    rc.def = r;

    // Use the same effective fee rate for industry service fees.
    rc.q1 = quoteIndustryOrder(*r,
                              originStation.id,
                              originStation.type,
                              1.0,
                              feeFrom,
                              params.processingRep);

    // Defensive: skip broken recipes.
    if (rc.q1.outputUnits <= 1e-9) continue;

    rc.maxByInputs = maxBatchesByInputInventory(fromEcon, originStation.economyModel, rc.q1, params.bidAskSpread);
    rc.maxByCargo = maxBatchesByOutputCargo(capKg, rc.q1);
    rc.maxOverall = std::max(0, std::min(rc.maxByInputs, rc.maxByCargo));
    if (rc.maxOverall <= 0) continue;

    caps.push_back(rc);
  }

  if (caps.empty()) return out;

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

      std::vector<IndustryTradeOpportunity> stationIdeas;

      for (const auto& rc : caps) {
        if (!rc.def) continue;
        const double outPerBatch = rc.q1.outputUnits;
        const int maxBySpace = maxBatchesByDestinationSpace(toEcon,
                                                            toSt.economyModel,
                                                            rc.q1.output,
                                                            outPerBatch,
                                                            params.bidAskSpread);
        const int batches = std::max(0, std::min(rc.maxOverall, maxBySpace));
        if (batches <= 0) continue;

        const auto q = quoteIndustryOrder(*rc.def,
                                          originStation.id,
                                          originStation.type,
                                          (double)batches,
                                          feeFrom,
                                          params.processingRep);

        if (q.outputUnits <= 1e-9) continue;

        IndustryTradeOpportunity t;
        t.toSystem = sys.stub.id;
        t.toStation = toSt.id;
        t.toSystemName = sys.stub.name;
        t.toStationName = toSt.name;
        t.recipe = rc.def->id;
        t.batches = q.batches;
        t.inputA = q.inputA;
        t.inputAUnits = q.inputAUnits;
        t.inputB = q.inputB;
        t.inputBUnits = q.inputBUnits;
        t.output = q.output;
        t.outputUnits = q.outputUnits;
        t.serviceFeeCr = q.serviceFeeCr;
        t.timeDays = q.timeDays;
        t.feeFrom = feeFrom;
        t.feeTo = feeTo;
        t.distanceLy = distLy;

        // Market cost/revenue.
        if (t.inputAUnits > 1e-9) {
          const auto qa = econ::quote(fromEcon, originStation.economyModel, t.inputA, params.bidAskSpread);
          t.inputAAsk = qa.ask;
          t.inputACostCr = qa.ask * t.inputAUnits * (1.0 + feeFrom);
        }
        if (t.inputBUnits > 1e-9) {
          const auto qb = econ::quote(fromEcon, originStation.economyModel, t.inputB, params.bidAskSpread);
          t.inputBAsk = qb.ask;
          t.inputBCostCr = qb.ask * t.inputBUnits * (1.0 + feeFrom);
        }
        {
          const auto qo = econ::quote(toEcon, toSt.economyModel, t.output, params.bidAskSpread);
          t.outputBid = qo.bid;
          t.outputRevenueCr = qo.bid * t.outputUnits * (1.0 - feeTo);
        }

        const double unitMassKg = econ::commodityDef(t.output).massKg;
        t.outputMassKg = std::max(0.0, t.outputUnits * unitMassKg);

        t.netProfitCr = t.outputRevenueCr - t.inputACostCr - t.inputBCostCr - t.serviceFeeCr;
        if (!std::isfinite(t.netProfitCr)) t.netProfitCr = 0.0;

        t.netProfitPerKg = (t.outputMassKg > 1e-9) ? (t.netProfitCr / t.outputMassKg) : t.netProfitCr;
        t.netProfitPerDay = (t.timeDays > 1e-9) ? (t.netProfitCr / t.timeDays) : t.netProfitCr;

        pushOpportunity(stationIdeas, t, params);
      }

      if (stationIdeas.empty()) continue;

      std::sort(stationIdeas.begin(), stationIdeas.end(), [](const IndustryTradeOpportunity& a,
                                                            const IndustryTradeOpportunity& b) {
        if (a.netProfitCr != b.netProfitCr) return a.netProfitCr > b.netProfitCr;
        if (a.netProfitPerDay != b.netProfitPerDay) return a.netProfitPerDay > b.netProfitPerDay;
        if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
        if (a.toSystem != b.toSystem) return a.toSystem < b.toSystem;
        if (a.toStation != b.toStation) return a.toStation < b.toStation;
        return static_cast<std::size_t>(a.recipe) < static_cast<std::size_t>(b.recipe);
      });

      const std::size_t take = std::min<std::size_t>(params.perStationLimit, stationIdeas.size());
      for (std::size_t i = 0; i < take; ++i) {
        out.push_back(std::move(stationIdeas[i]));
      }
    }
  }

  std::sort(out.begin(), out.end(), [](const IndustryTradeOpportunity& a, const IndustryTradeOpportunity& b) {
    if (a.netProfitCr != b.netProfitCr) return a.netProfitCr > b.netProfitCr;
    if (a.netProfitPerDay != b.netProfitPerDay) return a.netProfitPerDay > b.netProfitPerDay;
    if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
    if (a.toSystem != b.toSystem) return a.toSystem < b.toSystem;
    if (a.toStation != b.toStation) return a.toStation < b.toStation;
    return static_cast<std::size_t>(a.recipe) < static_cast<std::size_t>(b.recipe);
  });

  if (params.maxResults > 0 && out.size() > params.maxResults) {
    out.resize(params.maxResults);
  }

  return out;
}

std::vector<IndustryTradeOpportunity> scanIndustryTradeOpportunities(Universe& u,
                                                                     const SystemStub& originStub,
                                                                     const Station& originStation,
                                                                     double timeDays,
                                                                     const IndustryTradeScanParams& params,
                                                                     IndustryFeeRateFn feeRate) {
  const std::size_t maxSystems = std::max<std::size_t>(1, params.maxSystems);
  const auto candidates = u.queryNearby(originStub.posLy, params.radiusLy, maxSystems);
  return scanIndustryTradeOpportunities(u, originStub, originStation, timeDays, candidates, params, std::move(feeRate));
}

} // namespace stellar::sim
