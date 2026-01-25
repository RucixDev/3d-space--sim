#include "stellar/sim/SmugglingScanner.h"

#include "stellar/econ/Market.h"
#include "stellar/sim/BlackMarket.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Law.h"
#include "stellar/sim/PoliceScan.h"
#include "stellar/sim/SecurityModel.h"
#include "stellar/sim/SystemConditions.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double clampFee(double feeRate) {
  if (!std::isfinite(feeRate)) return 0.0;
  // Fees in this prototype are expected to be small (0..25%), but be generous.
  return std::clamp(feeRate, 0.0, 0.95);
}

static double effectiveCargoKg(const SmugglingScanParams& p) {
  double cap = p.cargoCapacityKg;
  if (!std::isfinite(cap)) cap = 0.0;
  cap = std::max(0.0, cap);

  if (!p.useFreeHold) return cap;

  double used = p.cargoUsedKg;
  if (!std::isfinite(used)) used = 0.0;
  used = std::max(0.0, used);
  return std::max(0.0, cap - used);
}

static SmugglingFeeRateFn defaultFeeFn() {
  return [](const Station& st) { return st.feeRate; };
}

static double heatMult(double playerHeat) {
  if (!std::isfinite(playerHeat)) return 1.0;
  const double h = std::clamp(playerHeat, 0.0, 100.0);
  return std::clamp(1.0 + 0.012 * h, 1.0, 2.20);
}

static double availabilityMul(const SmugglingScanParams& p, const BlackMarketProfile& bm) {
  switch (p.availability) {
    case SmugglingAvailabilityMode::TodayOnly: return bm.available ? 1.0 : 0.0;
    case SmugglingAvailabilityMode::Expected:  return std::clamp(bm.access01, 0.0, 1.0);
    case SmugglingAvailabilityMode::Ignore:    return 1.0;
    default:                                   return bm.available ? 1.0 : 0.0;
  }
}

static double scoreOne(const SmugglingScanParams& p,
                       double availMul,
                       double distanceLy,
                       double expectedProfitCr,
                       double cleanProfitCr,
                       double stdDevCr) {
  double score = 0.0;
  switch (p.scoreMode) {
    case SmugglingScoreMode::ExpectedProfit:
      score = expectedProfitCr;
      break;
    case SmugglingScoreMode::RiskAdjusted:
      score = expectedProfitCr - std::max(0.0, p.riskLambda) * std::max(0.0, stdDevCr);
      break;
    case SmugglingScoreMode::CleanProfit:
      score = cleanProfitCr;
      break;
    case SmugglingScoreMode::ExpectedProfitPerLy:
      score = expectedProfitCr / std::max(1e-6, std::max(0.0, distanceLy));
      break;
    default:
      score = expectedProfitCr;
      break;
  }

  return score * std::max(0.0, availMul);
}


std::vector<SmugglingOpportunity> scanSmugglingOpportunities(Universe& u,
                                                            const SystemStub& originStub,
                                                            const Station& originStation,
                                                            double timeDays,
                                                            const std::vector<SystemStub>& candidates,
                                                            const SmugglingScanParams& params,
                                                            SmugglingFeeRateFn feeRate) {
  std::vector<SmugglingOpportunity> out;

  if (!std::isfinite(timeDays)) return out;
  if (originStub.id == 0 || originStation.id == 0) return out;
  if (candidates.empty()) return out;

  const double capKg = effectiveCargoKg(params);
  if (capKg <= 0.0) return out;

  if (!feeRate) feeRate = defaultFeeFn();
  const double feeFrom = clampFee(feeRate(originStation));

  const core::u64 seed = u.seed();

  auto& fromEcon = u.stationEconomy(originStation, timeDays);

  out.reserve(std::min<std::size_t>(params.maxResults * 4, candidates.size() * 2));

  // Pre-compute the origin station illegality mask so we can quickly check origin legality.
  const core::u32 originIllegalMask = illegalCommodityMaskForStation(seed,
                                                                     originStation.factionId,
                                                                     originStation.id,
                                                                     originStation.type);

  const std::size_t maxBits = std::min<std::size_t>(econ::kCommodityCount, 32);

  for (const auto& stub : candidates) {
    if (stub.id == 0) continue;
    if (!params.includeSameSystem && stub.id == originStub.id) continue;

    const auto& sys = u.getSystem(stub.id, &stub);

    const double distLy = (stub.posLy - originStub.posLy).length();
    SystemEvent sysEv{};
    SystemSecurityProfile sec{};
    if (params.useLiveSystemConditions) {
      const auto* deltaMap = u.systemSecurityDeltaMap();
      const SystemSecurityDeltaState* deltaState = nullptr;
      if (deltaMap) {
        const auto it = deltaMap->find(sys.stub.id);
        if (it != deltaMap->end()) deltaState = &it->second;
      }
      sec = effectiveSystemSecurityProfile(seed,
                                         sys,
                                         timeDays,
                                         deltaState,
                                         u.systemSecurityDynamicsParams(),
                                         u.systemEventParams(),
                                         &sysEv);
    } else {
      sec = systemSecurityProfile(seed, sys);
    }

    for (const auto& toSt : sys.stations) {
      if (toSt.id == 0) continue;
      if (stub.id == originStub.id && toSt.id == originStation.id) continue;

      // Resolve black market conditions at the destination.
      const LawProfile law = lawProfile(seed, toSt.factionId);
      double rep = params.playerRep;
      if (params.repForFaction) {
        rep = params.repForFaction(toSt.factionId);
      }
      const BlackMarketProfile bm = blackMarketProfile(seed,
                                                       sys.stub.id,
                                                       toSt,
                                                       sec,
                                                       law,
                                                       timeDays,
                                                       rep);
      const double aMul = availabilityMul(params, bm);
      if (aMul <= 1e-9) continue;

      // Destination illegality mask defines which goods are smuggle-worthy here.
      const core::u32 illegalMask = illegalCommodityMaskForStation(seed, toSt.factionId, toSt.id, toSt.type);
      if (illegalMask == 0u) continue;

      auto& toEcon = u.stationEconomy(toSt, timeDays);

      std::vector<SmugglingOpportunity> stationOpps;
      stationOpps.reserve(8);

      for (std::size_t i = 0; i < maxBits; ++i) {
        if ((illegalMask & ((core::u32)1u << (core::u32)i)) == 0u) continue;
        const auto cid = static_cast<econ::CommodityId>(i);

        if (params.requireOriginLegal) {
          if ((originIllegalMask & ((core::u32)1u << (core::u32)i)) != 0u) {
            continue;
          }
        }

        const auto fromQ = econ::quote(fromEcon, originStation.economyModel, cid, params.bidAskSpread);
        const auto toQ = econ::quote(toEcon, toSt.economyModel, cid, params.bidAskSpread);
        const auto blackQ = applyBlackMarketQuote(toQ, bm);

        double buyAsk = fromQ.ask;
        if (!std::isfinite(buyAsk) || buyAsk <= 1e-9) continue;
        double sellBid = blackQ.bid;
        if (!std::isfinite(sellBid) || sellBid <= 1e-9) continue;

        const double buyAskNet = buyAsk * (1.0 + feeFrom);
        if (!std::isfinite(buyAskNet) || buyAskNet <= 1e-9) continue;

        const double unitsFrom = std::max(0.0, fromQ.inventory);

        const double capTo = std::max(0.0, toSt.economyModel.capacity[i]);
        double invTo = toEcon.inventory[i];
        if (!std::isfinite(invTo)) invTo = 0.0;
        invTo = std::clamp(invTo, 0.0, capTo);
        const double unitsToSpace = std::max(0.0, capTo - invTo);

        const double unitMassKg = econ::commodityDef(cid).massKg;
        if (!(unitMassKg > 1e-12) || !std::isfinite(unitMassKg)) continue;

        const double unitsByCargo = capKg / unitMassKg;
        const double unitsPossible = std::max(0.0, std::min({unitsFrom, unitsToSpace, unitsByCargo}));
        if (unitsPossible <= 1e-9) continue;

        const double buyCost = buyAskNet * unitsPossible;
        const double payout = sellBid * unitsPossible;
        const double cleanProfit = payout - buyCost;
        if (!(cleanProfit > 1e-6) || !std::isfinite(cleanProfit)) continue;

        // Effective sting probability.
        double pSting = bm.stingChance;
        if (!std::isfinite(pSting)) pSting = 0.0;
        pSting = std::clamp(pSting, 0.0, 1.0);

        pSting *= heatMult(params.playerHeat);
        pSting *= smuggleHoldScanRateMult(params.smuggleHoldMk);
        pSting = std::clamp(pSting, 0.0, 1.0);

        double mid = toQ.mid;
        if (!std::isfinite(mid) || mid <= 1e-9) mid = econ::commodityDef(cid).basePrice;
        const double illegalValue = unitsPossible * std::max(0.0, mid);
        const double fine = law.contrabandFineCr(illegalValue);

        const double stungProfit = -buyCost - fine;
        const double expectedProfit = (1.0 - pSting) * cleanProfit + pSting * stungProfit;

        // For a Bernoulli outcome distribution, var = p(1-p)*(a-b)^2.
        const double diff = cleanProfit - stungProfit;
        const double var = pSting * (1.0 - pSting) * diff * diff;
        const double stdDev = std::sqrt(std::max(0.0, var));

        const double score = scoreOne(params, aMul, distLy, expectedProfit, cleanProfit, stdDev);
        if (!std::isfinite(score)) continue;
        if (score + 1e-9 < params.minScoreCr) continue;

        SmugglingOpportunity t{};
        t.toSystem = sys.stub.id;
        t.toStation = toSt.id;
        t.toSystemName = sys.stub.name;
        t.toStationName = toSt.name;
        t.commodity = cid;

        t.buyAskCr = buyAsk;
        t.buyAskNetCr = buyAskNet;
        t.sellBidCr = sellBid;
        t.officialSellBidCr = toQ.bid;

        t.unitsFrom = unitsFrom;
        t.unitsToSpace = unitsToSpace;
        t.unitsPossible = unitsPossible;
        t.unitMassKg = unitMassKg;
        t.distanceLy = distLy;

        t.systemEventKind = (params.useLiveSystemConditions && sysEv.active) ? sysEv.kind : SystemEventKind::None;
        t.systemEventSeverity01 = (params.useLiveSystemConditions && sysEv.active) ? sysEv.severity01 : 0.0;
        t.systemSecurity01 = sec.security01;
        t.systemPiracy01 = sec.piracy01;
        t.systemTraffic01 = sec.traffic01;

        t.blackMarketAvailable = bm.available;
        t.blackMarketAccess01 = bm.access01;
        t.stingChance = pSting;
        t.illegalValueCr = illegalValue;
        t.fineCr = fine;

        t.buyCostCr = buyCost;
        t.payoutCr = payout;
        t.cleanProfitCr = cleanProfit;
        t.stungProfitCr = stungProfit;
        t.expectedProfitCr = expectedProfit;
        t.stdDevCr = stdDev;
        t.scoreCr = score;

        stationOpps.push_back(std::move(t));
      }

      if (stationOpps.empty()) continue;

      std::sort(stationOpps.begin(), stationOpps.end(), [](const SmugglingOpportunity& a, const SmugglingOpportunity& b) {
        if (std::fabs(a.scoreCr - b.scoreCr) > 1e-6) return a.scoreCr > b.scoreCr;
        if (std::fabs(a.expectedProfitCr - b.expectedProfitCr) > 1e-6) return a.expectedProfitCr > b.expectedProfitCr;
        if (std::fabs(a.cleanProfitCr - b.cleanProfitCr) > 1e-6) return a.cleanProfitCr > b.cleanProfitCr;
        if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
        if (a.toSystem != b.toSystem) return a.toSystem < b.toSystem;
        if (a.toStation != b.toStation) return a.toStation < b.toStation;
        return static_cast<std::size_t>(a.commodity) < static_cast<std::size_t>(b.commodity);
      });

      const std::size_t take = std::max<std::size_t>(1, params.perStationLimit);
      if (stationOpps.size() > take) stationOpps.resize(take);

      for (auto& t : stationOpps) out.push_back(std::move(t));
    }
  }

  if (out.empty()) return out;

  std::sort(out.begin(), out.end(), [](const SmugglingOpportunity& a, const SmugglingOpportunity& b) {
    if (std::fabs(a.scoreCr - b.scoreCr) > 1e-6) return a.scoreCr > b.scoreCr;
    if (std::fabs(a.expectedProfitCr - b.expectedProfitCr) > 1e-6) return a.expectedProfitCr > b.expectedProfitCr;
    if (std::fabs(a.cleanProfitCr - b.cleanProfitCr) > 1e-6) return a.cleanProfitCr > b.cleanProfitCr;
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


std::vector<SmugglingOpportunity> scanSmugglingOpportunities(Universe& u,
                                                            const SystemStub& originStub,
                                                            const Station& originStation,
                                                            double timeDays,
                                                            const SmugglingScanParams& params,
                                                            SmugglingFeeRateFn feeRate) {
  const std::size_t maxSystems = std::max<std::size_t>(1, params.maxSystems);
  const auto candidates = u.queryNearby(originStub.posLy, params.radiusLy, maxSystems);
  return scanSmugglingOpportunities(u, originStub, originStation, timeDays, candidates, params, std::move(feeRate));
}

} // namespace stellar::sim
