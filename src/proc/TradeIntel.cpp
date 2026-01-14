#include "stellar/proc/TradeIntel.h"

#include "stellar/econ/Commodity.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {

static const sim::SystemStub* findStub(sim::SystemId id,
                                      const std::vector<sim::SystemStub>& stubs) {
  for (const auto& s : stubs) {
    if (s.id == id) return &s;
  }
  return nullptr;
}

static const TradeProfile* findProfile(sim::SystemId id,
                                      const std::vector<sim::SystemStub>& stubs,
                                      const std::vector<TradeProfile>& profiles) {
  const std::size_t n = std::min(stubs.size(), profiles.size());
  for (std::size_t i = 0; i < n; ++i) {
    if (stubs[i].id == id) return &profiles[i];
  }
  return nullptr;
}

static TradeRouteEconomics enrichRoute(TradeIntelDirection dir,
                                       const sim::SystemStub& origin,
                                       const TradeProfile& originP,
                                       const sim::SystemStub& other,
                                       const TradeProfile& otherP,
                                       const TradeRoute& r,
                                       const TradeIntelParams& ip,
                                       const TradePriceModelParams& pp) {
  TradeRouteEconomics e{};
  e.direction = dir;
  e.otherSystem = r.otherSystem;
  e.commodity = r.commodity;
  e.score = r.score;
  e.distanceLy = r.distanceLy;
  e.directDistanceLy = r.directDistanceLy;
  e.laneDistanceLy = r.laneDistanceLy;
  e.laneHops = r.laneHops;
  e.laneBottleneckBandwidth01 = r.laneBottleneckBandwidth01;
  e.risk01 = r.risk01;

  const auto& def = econ::commodityDef(r.commodity);
  e.unitMassKg = std::max(1e-9, def.massKg);
  e.unitsForCargo = std::max(0.0, ip.cargoCapacityKg) / e.unitMassKg;

  // Direction determines which side we buy/sell.
  TradePriceQuote buyQ{};
  TradePriceQuote sellQ{};

  if (dir == TradeIntelDirection::Export) {
    buyQ = estimateTradePriceQuote(origin.seed, originP, r.commodity, ip.bidAskSpread, pp);
    sellQ = estimateTradePriceQuote(other.seed, otherP, r.commodity, ip.bidAskSpread, pp);
  } else {
    buyQ = estimateTradePriceQuote(other.seed, otherP, r.commodity, ip.bidAskSpread, pp);
    sellQ = estimateTradePriceQuote(origin.seed, originP, r.commodity, ip.bidAskSpread, pp);
  }

  e.buyMidCr = buyQ.mid;
  e.buyAskCr = buyQ.ask * (1.0 + std::max(0.0, ip.buyFeeRate));
  e.sellMidCr = sellQ.mid;
  e.sellBidCr = sellQ.bid * (1.0 - std::max(0.0, ip.sellFeeRate));

  e.profitPerUnitCr = e.sellBidCr - e.buyAskCr;
  e.profitPerKgCr = e.profitPerUnitCr / e.unitMassKg;
  e.profitForCargoCr = e.profitPerUnitCr * e.unitsForCargo;

  return e;
}

TradeIntelReport buildTradeIntel(const sim::SystemStub& origin,
                                 const TradeProfile& originProfile,
                                 const std::vector<sim::SystemStub>& candidates,
                                 const std::vector<TradeProfile>& candidateProfiles,
                                 const TradeRouteParams& routeParams,
                                 const TradeIntelParams& intelParams,
                                 const TradePriceModelParams& priceParams) {
  TradeIntelReport rep{};

  if (candidates.size() != candidateProfiles.size()) return rep;

  // Clamp some user-driven params defensively.
  TradeIntelParams ip = intelParams;
  if (!std::isfinite(ip.bidAskSpread)) ip.bidAskSpread = 0.10;
  ip.bidAskSpread = std::clamp(ip.bidAskSpread, 0.0, 1.0);
  if (!std::isfinite(ip.buyFeeRate)) ip.buyFeeRate = 0.0;
  if (!std::isfinite(ip.sellFeeRate)) ip.sellFeeRate = 0.0;
  ip.buyFeeRate = std::clamp(ip.buyFeeRate, 0.0, 1.0);
  ip.sellFeeRate = std::clamp(ip.sellFeeRate, 0.0, 1.0);
  if (!std::isfinite(ip.cargoCapacityKg)) ip.cargoCapacityKg = 0.0;
  ip.cargoCapacityKg = std::max(0.0, ip.cargoCapacityKg);

  auto sugg = suggestTradeRoutes(origin, originProfile, candidates, candidateProfiles, routeParams);

  rep.exports.reserve(sugg.exports.size());
  rep.imports.reserve(sugg.imports.size());

  for (const auto& r : sugg.exports) {
    const sim::SystemStub* other = findStub(r.otherSystem, candidates);
    const TradeProfile* otherP = findProfile(r.otherSystem, candidates, candidateProfiles);
    if (!other || !otherP) continue;
    rep.exports.push_back(enrichRoute(TradeIntelDirection::Export, origin, originProfile, *other, *otherP, r, ip, priceParams));
  }

  for (const auto& r : sugg.imports) {
    const sim::SystemStub* other = findStub(r.otherSystem, candidates);
    const TradeProfile* otherP = findProfile(r.otherSystem, candidates, candidateProfiles);
    if (!other || !otherP) continue;
    rep.imports.push_back(enrichRoute(TradeIntelDirection::Import, origin, originProfile, *other, *otherP, r, ip, priceParams));
  }

  // ---- 2-leg loops: match exports/imports by otherSystem ----
  for (const auto& ex : rep.exports) {
    for (const auto& im : rep.imports) {
      if (ex.otherSystem != im.otherSystem) continue;

      TradeLoop2 loop{};
      loop.otherSystem = ex.otherSystem;
      loop.outLeg = ex;
      loop.backLeg = im;
      loop.roundTripProfitCr = ex.profitForCargoCr + im.profitForCargoCr;
      loop.roundTripDistanceLy = 2.0 * ex.distanceLy;
      loop.avgRisk01 = 0.5 * (ex.risk01 + im.risk01);

      rep.loops.push_back(loop);
      break;
    }
  }

  std::sort(rep.loops.begin(), rep.loops.end(), [](const TradeLoop2& a, const TradeLoop2& b) {
    if (a.roundTripProfitCr != b.roundTripProfitCr) return a.roundTripProfitCr > b.roundTripProfitCr;
    if (a.roundTripDistanceLy != b.roundTripDistanceLy) return a.roundTripDistanceLy < b.roundTripDistanceLy;
    return a.otherSystem < b.otherSystem;
  });

  if (rep.loops.size() > ip.maxLoops) rep.loops.resize(ip.maxLoops);

  return rep;
}

TradeIntelReport buildTradeIntel(const sim::SystemStub& origin,
                                 const TradeProfile& originProfile,
                                 const std::vector<sim::SystemStub>& candidates,
                                 const std::vector<TradeProfile>& candidateProfiles,
                                 const HyperlaneNetwork& hyperlanes,
                                 const TradeRouteParams& routeParams,
                                 const HyperlaneTravelParams& travelParams,
                                 const TradeIntelParams& intelParams,
                                 const TradePriceModelParams& priceParams) {
  TradeIntelReport rep{};

  if (candidates.size() != candidateProfiles.size()) return rep;

  // Clamp some user-driven params defensively.
  TradeIntelParams ip = intelParams;
  if (!std::isfinite(ip.bidAskSpread)) ip.bidAskSpread = 0.10;
  ip.bidAskSpread = std::clamp(ip.bidAskSpread, 0.0, 1.0);
  if (!std::isfinite(ip.buyFeeRate)) ip.buyFeeRate = 0.0;
  if (!std::isfinite(ip.sellFeeRate)) ip.sellFeeRate = 0.0;
  ip.buyFeeRate = std::clamp(ip.buyFeeRate, 0.0, 1.0);
  ip.sellFeeRate = std::clamp(ip.sellFeeRate, 0.0, 1.0);
  if (!std::isfinite(ip.cargoCapacityKg)) ip.cargoCapacityKg = 0.0;
  ip.cargoCapacityKg = std::max(0.0, ip.cargoCapacityKg);

  auto sugg = suggestTradeRoutes(origin, originProfile, candidates, candidateProfiles, hyperlanes, routeParams, travelParams);

  rep.exports.reserve(sugg.exports.size());
  rep.imports.reserve(sugg.imports.size());

  for (const auto& r : sugg.exports) {
    const sim::SystemStub* other = findStub(r.otherSystem, candidates);
    const TradeProfile* otherP = findProfile(r.otherSystem, candidates, candidateProfiles);
    if (!other || !otherP) continue;
    rep.exports.push_back(enrichRoute(TradeIntelDirection::Export, origin, originProfile, *other, *otherP, r, ip, priceParams));
  }

  for (const auto& r : sugg.imports) {
    const sim::SystemStub* other = findStub(r.otherSystem, candidates);
    const TradeProfile* otherP = findProfile(r.otherSystem, candidates, candidateProfiles);
    if (!other || !otherP) continue;
    rep.imports.push_back(enrichRoute(TradeIntelDirection::Import, origin, originProfile, *other, *otherP, r, ip, priceParams));
  }

  // ---- 2-leg loops: match exports/imports by otherSystem ----
  for (const auto& ex : rep.exports) {
    for (const auto& im : rep.imports) {
      if (ex.otherSystem != im.otherSystem) continue;

      TradeLoop2 loop{};
      loop.otherSystem = ex.otherSystem;
      loop.outLeg = ex;
      loop.backLeg = im;
      loop.roundTripProfitCr = ex.profitForCargoCr + im.profitForCargoCr;
      loop.roundTripDistanceLy = 2.0 * ex.distanceLy;
      loop.avgRisk01 = 0.5 * (ex.risk01 + im.risk01);

      rep.loops.push_back(loop);
      break;
    }
  }

  std::sort(rep.loops.begin(), rep.loops.end(), [](const TradeLoop2& a, const TradeLoop2& b) {
    if (a.roundTripProfitCr != b.roundTripProfitCr) return a.roundTripProfitCr > b.roundTripProfitCr;
    if (a.roundTripDistanceLy != b.roundTripDistanceLy) return a.roundTripDistanceLy < b.roundTripDistanceLy;
    return a.otherSystem < b.otherSystem;
  });

  if (rep.loops.size() > ip.maxLoops) rep.loops.resize(ip.maxLoops);

  return rep;
}

} // namespace stellar::proc
