#include "stellar/sim/BlackMarket.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double stationOpenness(econ::StationType t) {
  // Higher => fences are easier to find / more common.
  using econ::StationType;
  switch (t) {
    case StationType::Outpost:       return 0.78;
    case StationType::Mining:        return 0.70;
    case StationType::Refinery:      return 0.66;
    case StationType::Agricultural:  return 0.58;
    case StationType::Industrial:    return 0.55;
    case StationType::TradeHub:      return 0.48;
    case StationType::Research:      return 0.40;
    case StationType::Shipyard:      return 0.35;
    default:                         return 0.55;
  }
}

static double stationRiskBias(econ::StationType t) {
  // Higher => more cameras/patrols, harder to operate.
  using econ::StationType;
  switch (t) {
    case StationType::Shipyard:      return 1.20;
    case StationType::Research:      return 1.15;
    case StationType::TradeHub:      return 1.05;
    case StationType::Industrial:    return 1.00;
    case StationType::Refinery:      return 0.98;
    case StationType::Agricultural:  return 0.95;
    case StationType::Mining:        return 0.92;
    case StationType::Outpost:       return 0.88;
    default:                         return 1.00;
  }
}

BlackMarketProfile blackMarketProfile(core::u64 universeSeed,
                                     SystemId systemId,
                                     const Station& station,
                                     const SystemSecurityProfile& sec,
                                     const LawProfile& law,
                                     double timeDays,
                                     double playerRep) {
  BlackMarketProfile bm{};

  // Independent jurisdiction: allow the profile math to run, but availability
  // will generally be low and eligibility will still be checked per-commodity.
  const double openness = stationOpenness(station.type);

  const double corruption = std::clamp(law.corruption, 0.0, 1.0);
  const double piracy = std::clamp(sec.piracy01, 0.0, 1.0);
  const double security = std::clamp(sec.security01, 0.0, 1.0);
  const double contest = std::clamp(sec.contest01, 0.0, 1.0);
  const double control = std::clamp(sec.controlFrac, 0.0, 1.0);

  // ---- Access -------------------------------------------------------------
  // Intuition: fences thrive when:
  //  - stations are less formal (openness)
  //  - the local regime is corrupt
  //  - piracy is present and control is contested
  //  - security is not hyper-effective
  double baseAccess = openness;
  baseAccess *= (0.35 + 0.65 * corruption);
  baseAccess *= (0.45 + 0.55 * piracy);
  baseAccess *= (0.40 + 0.60 * (1.0 - security));
  baseAccess *= (0.70 + 0.30 * contest);
  baseAccess *= (0.85 + 0.15 * (1.0 - control));

  baseAccess = std::clamp(baseAccess * 1.50, 0.0, 1.0);
  const double repMul = std::clamp(0.85 + playerRep / 250.0, 0.50, 1.40);
  bm.access01 = std::clamp(baseAccess * repMul, 0.0, 1.0);

  // Daily availability roll. We use floor(timeDays) so that availability doesn't
  // flicker within a day.
  const core::u64 day = (timeDays > 0.0 && std::isfinite(timeDays))
                          ? (core::u64)std::floor(timeDays)
                          : 0ull;

  core::u64 s = core::hashCombine(universeSeed, core::seedFromText("black_market_v1"));
  s = core::hashCombine(s, (core::u64)systemId);
  s = core::hashCombine(s, (core::u64)station.id);
  s = core::hashCombine(s, day);

  core::SplitMix64 rng(s);
  bm.available = rng.nextUnit() < bm.access01;

  // ---- Risk ---------------------------------------------------------------
  // More security and stricter scanning makes deals riskier.
  const double strict = std::clamp(law.scanStrictness, 0.5, 2.0);
  double risk = 0.22;
  risk += 0.55 * security;
  risk += 0.25 * std::clamp(strict - 1.0, 0.0, 1.0);
  risk += 0.20 * (1.0 - corruption);
  risk *= stationRiskBias(station.type);
  risk *= (0.90 + 0.20 * control);

  // Good rep reduces scrutiny; bad rep increases it.
  const double repRisk = std::clamp(1.0 - playerRep / 400.0, 0.70, 1.30);
  bm.risk01 = std::clamp(risk * repRisk, 0.0, 1.0);

  // ---- Pricing ------------------------------------------------------------
  // Scarcity is higher when security is high and piracy is low.
  double scarcity = 0.55 * security;
  scarcity += 0.25 * std::clamp(strict - 1.0, 0.0, 1.0);
  scarcity += 0.20 * (1.0 - piracy);
  scarcity = std::clamp(scarcity, 0.0, 1.0);

  bm.bidMul = std::clamp(1.05 + 0.55 * scarcity, 0.80, 2.00);
  bm.askMul = std::clamp(1.10 + 0.65 * scarcity + 0.15 * bm.risk01, 0.85, 2.50);

  // Fence cut: grows with risk and low corruption.
  bm.fenceCut = 0.08 + 0.22 * bm.risk01 + 0.10 * (1.0 - corruption);
  bm.fenceCut = std::clamp(bm.fenceCut, 0.05, 0.45);

  // ---- Sting chance -------------------------------------------------------
  // Intuition: stings are more likely when enforcement is competent and you are
  // not well-connected.
  double sting = 0.05;
  sting += 0.35 * bm.risk01;
  sting += 0.10 * (security - piracy);
  sting *= (0.65 + 0.35 * (1.0 - corruption));
  sting *= (1.15 - 0.65 * bm.access01);
  sting *= repRisk;
  bm.stingChance = std::clamp(sting, 0.0, 0.85);

  return bm;
}

econ::MarketQuote applyBlackMarketQuote(const econ::MarketQuote& official,
                                       const BlackMarketProfile& bm) {
  econ::MarketQuote out = official;

  // Treat mid as the average of (ask,bid) so that consumers can display one
  // consistent "mid".
  const double bid = std::max(0.0, official.bid * bm.bidMul * (1.0 - bm.fenceCut));
  const double ask = std::max(0.0, official.ask * bm.askMul * (1.0 + 0.25 * bm.fenceCut));

  out.bid = bid;
  out.ask = ask;
  out.mid = 0.5 * (bid + ask);

  // Inventory is not modeled for the fence (off-books). We keep the official
  // inventory for UI convenience.
  return out;
}

bool rollBlackMarketSting(core::u64 eventSeed,
                          const BlackMarketProfile& bm,
                          double playerHeat) {
  double p = bm.stingChance;
  if (!std::isfinite(p)) p = 0.0;
  p = std::clamp(p, 0.0, 1.0);

  if (std::isfinite(playerHeat)) {
    // 0 heat -> 1.00x, 100 heat -> ~2.20x (clamped).
    const double h = std::clamp(playerHeat, 0.0, 100.0);
    p *= std::clamp(1.0 + 0.012 * h, 1.0, 2.20);
    p = std::clamp(p, 0.0, 1.0);
  }

  core::SplitMix64 rng(core::hashCombine(eventSeed, core::seedFromText("bm_sting")));
  return rng.nextUnit() < p;
}

BlackMarketSellResult sellToBlackMarket(core::u64 universeSeed,
                                        core::u64 eventSeed,
                                        const Station& station,
                                        const BlackMarketProfile& bm,
                                        const LawProfile& law,
                                        double playerHeat,
                                        econ::CommodityId commodity,
                                        double units,
                                        double bidAskSpread,
                                        const std::array<double, econ::kCommodityCount>* midPriceOverrideCr,
                                        double& credits,
                                        std::array<double, econ::kCommodityCount>& cargoUnits) {
  BlackMarketSellResult out{};
  out.commodity = commodity;
  out.intendedUnits = units;

  if (!(units > 0.0) || !std::isfinite(units)) {
    out.reason = "units<=0";
    return out;
  }

  if (!bm.available) {
    out.reason = "no black market";
    return out;
  }

  if (!blackMarketEligibleCommodity(universeSeed, station, commodity)) {
    out.reason = "commodity not eligible";
    return out;
  }

  const std::size_t idx = (std::size_t)commodity;
  if (idx >= econ::kCommodityCount) {
    out.reason = "invalid commodity";
    return out;
  }

  double have = cargoUnits[idx];
  if (!std::isfinite(have)) have = 0.0;
  if (have <= 1e-6) {
    out.reason = "no cargo";
    return out;
  }

  const double sellUnits = std::min(units, have);
  if (sellUnits <= 1e-6) {
    out.reason = "no cargo";
    return out;
  }
  out.unitsSold = sellUnits;

  const double creditsBefore = credits;

  // Roll the sting before moving cargo/credits.
  out.stung = rollBlackMarketSting(core::hashCombine(eventSeed, (core::u64)idx), bm, playerHeat);

  if (out.stung) {
    const core::u32 mask = illegalCommodityMaskForStation(universeSeed, station.factionId, station.id, station.type);
    out.scan = scanIllegalCargoMask(mask, cargoUnits, midPriceOverrideCr);
    out.enforcement = enforceContraband(law, credits, cargoUnits, out.scan.scannedIllegalUnits, out.scan.illegalValueCr);

    credits = out.enforcement.creditsAfter;
    cargoUnits = out.enforcement.cargoAfter;

    out.creditsDelta = credits - creditsBefore;
    out.ok = true;
    return out;
  }

  // Price reference: mid -> bid (spread/2), then apply bm multipliers.
  double mid = econ::commodityDef(commodity).basePrice;
  if (midPriceOverrideCr) {
    const double m = (*midPriceOverrideCr)[idx];
    if (std::isfinite(m) && m > 0.0) mid = m;
  }
  if (!std::isfinite(mid) || mid < 0.0) mid = 0.0;

  if (!std::isfinite(bidAskSpread)) bidAskSpread = 0.10;
  bidAskSpread = std::clamp(bidAskSpread, 0.0, 1.0);

  const double officialBid = mid * (1.0 - 0.5 * bidAskSpread);
  const double unitPrice = std::max(0.0, officialBid * bm.bidMul * (1.0 - bm.fenceCut));
  const double payout = unitPrice * sellUnits;

  credits += payout;
  cargoUnits[idx] = std::max(0.0, have - sellUnits);

  out.pricePerUnitCr = unitPrice;
  out.payoutCr = payout;
  out.creditsDelta = credits - creditsBefore;
  out.ok = true;
  return out;
}


BlackMarketBuyResult buyFromBlackMarket(core::u64 universeSeed,
                                        core::u64 eventSeed,
                                        const Station& station,
                                        const BlackMarketProfile& bm,
                                        const LawProfile& law,
                                        double playerHeat,
                                        econ::CommodityId commodity,
                                        double units,
                                        double bidAskSpread,
                                        const std::array<double, econ::kCommodityCount>* midPriceOverrideCr,
                                        double& credits,
                                        std::array<double, econ::kCommodityCount>& cargoUnits) {
  BlackMarketBuyResult out{};
  out.commodity = commodity;
  out.intendedUnits = units;

  if (units <= 1e-9) {
    out.reason = "units<=0";
    return out;
  }
  if (!bm.available) {
    out.reason = "no_black_market";
    return out;
  }

  const std::size_t idx = (std::size_t)commodity;
  if (idx >= econ::kCommodityCount) {
    out.reason = "invalid_commodity";
    return out;
  }

  if (!blackMarketEligibleCommodity(universeSeed, station, commodity)) {
    out.reason = "not_eligible";
    return out;
  }

  const double creditsBefore = credits;

  const double spread = std::clamp(bidAskSpread, 0.0, 0.95);

  // Resolve a reference mid price.
  double mid = econ::commodityDef(commodity).basePrice;
  if (midPriceOverrideCr) {
    mid = std::max(0.0, (*midPriceOverrideCr)[idx]);
    if (mid <= 1e-9) mid = econ::commodityDef(commodity).basePrice;
  }

  const double officialAsk = mid * (1.0 + 0.5 * spread);
  const double fenceCut = std::clamp(bm.fenceCut, 0.0, 1.0);

  double unitPrice = officialAsk * std::max(0.01, bm.askMul) * (1.0 + 0.25 * fenceCut);
  if (!(unitPrice > 1e-9) || !std::isfinite(unitPrice)) {
    out.reason = "invalid_price";
    return out;
  }

  // Clamp purchase to available credits. (Cargo mass limits are owned by the caller/UI.)
  const double affordable = std::max(0.0, credits) / unitPrice;
  const double buyUnits = std::min(std::max(0.0, units), affordable);

  if (buyUnits <= 1e-9) {
    out.reason = "need_credits";
    return out;
  }

  out.unitsBought = buyUnits;
  const double cost = unitPrice * buyUnits;

  out.pricePerUnitCr = unitPrice;
  out.costCr = cost;

  // Roll sting using a per-commodity salted seed so buy/sell outcomes don't correlate.
  const core::u64 salt = core::seedFromText("bm_buy");
  out.stung = rollBlackMarketSting(core::hashCombine(eventSeed, core::hashCombine(salt, (core::u64)idx)), bm, playerHeat);

  // Apply the purchase first (even on a sting). In a sting narrative, the player is
  // assumed to complete the exchange and then enforcement moves in.
  credits = std::max(0.0, credits) - cost;
  cargoUnits[idx] = std::max(0.0, cargoUnits[idx]) + buyUnits;

  if (out.stung) {
    // Enforce contraband on the *full* cargo under this station's illegality mask.
    const core::u32 mask = illegalCommodityMaskForStation(universeSeed, station.factionId, station.id, station.type);

    out.scan = scanIllegalCargoMask(mask, cargoUnits, midPriceOverrideCr);
    out.enforcement = enforceContraband(law, credits, cargoUnits, out.scan.scannedIllegalUnits, out.scan.illegalValueCr);

    credits = out.enforcement.creditsAfter;
    cargoUnits = out.enforcement.cargoAfter;

    out.creditsDelta = credits - creditsBefore;
    out.ok = true;
    return out;
  }

  out.creditsDelta = credits - creditsBefore;
  out.ok = true;
  return out;
}


} // namespace stellar::sim
