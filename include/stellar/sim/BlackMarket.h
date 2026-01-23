#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/econ/Market.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Law.h"
#include "stellar/sim/PoliceScan.h"
#include "stellar/sim/SecurityModel.h"
#include "stellar/sim/System.h"

#include <array>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Black Market (core, headless)
// -----------------------------------------------------------------------------
//
// The prototype already supports:
//   - per-faction + station-local contraband legality masks (Contraband.h)
//   - deterministic scans, fines, bribes, and confiscation math (PoliceScan.h)
//
// What's been missing is a *profitable* and *system-sensitive* way to actually
// move illegal goods, turning "contraband" from a pure penalty into a game loop.
//
// This module provides deterministic, bounded, faction/system flavored black
// market mechanics:
//   - whether a station has an accessible fence today
//   - how good the prices are (bid/ask multipliers + fence cut)
//   - how risky a deal is (sting chance)
//
// The goal is to make "lawless" space meaningfully different without requiring
// any renderer/UI dependencies.

struct BlackMarketProfile {
  bool available{false};

  // How easy it is to find/use a fence here (0..1).
  double access01{0.0};

  // Enforcement risk (0..1). Roughly: higher security + stricter law + lower
  // corruption => higher risk.
  double risk01{0.0};

  // Fence cut fraction (0..1). Applied to player sell payouts.
  double fenceCut{0.0};

  // Price multipliers applied to the *official* market quote.
  //
  // - bidMul increases what the fence pays relative to official bid.
  // - askMul increases what the fence charges relative to official ask.
  double bidMul{1.0};
  double askMul{1.0};

  // Baseline per-transaction sting chance (0..1).
  // Callers can optionally modulate this with ship heat.
  double stingChance{0.0};
};

// Compute a deterministic black market profile for a station.
//
// Inputs:
//  - universeSeed: world seed.
//  - systemId: used only for stable daily availability rolls.
//  - station: the station being queried.
//  - sec: system-level security/piracy/traffic metrics.
//  - law: law profile of the jurisdiction (usually station.factionId).
//  - timeDays: used to vary daily availability deterministically (floor(timeDays)).
//  - playerRep: reputation with the station faction; affects access/risk.
BlackMarketProfile blackMarketProfile(core::u64 universeSeed,
                                     SystemId systemId,
                                     const Station& station,
                                     const SystemSecurityProfile& sec,
                                     const LawProfile& law,
                                     double timeDays,
                                     double playerRep);

// True if `cid` is illegal at this station under deterministic contraband rules
// (and thus eligible for black market trading).
inline bool blackMarketEligibleCommodity(core::u64 universeSeed,
                                        const Station& station,
                                        econ::CommodityId cid) {
  return isIllegalCommodityAtStation(universeSeed, station.factionId, station.id, station.type, cid);
}

// Apply black market multipliers to an official market quote.
//
// NOTE: This does *not* check legality or availability; it's purely arithmetic.
econ::MarketQuote applyBlackMarketQuote(const econ::MarketQuote& official,
                                       const BlackMarketProfile& bm);

// Deterministic roll for whether a transaction is "stung" by enforcement.
//
// `eventSeed` should come from the caller's RNG/event stream so that
// transactions don't all share the same outcome.
//
// `playerHeat` (0..~100) increases sting probability in a bounded way.
bool rollBlackMarketSting(core::u64 eventSeed,
                          const BlackMarketProfile& bm,
                          double playerHeat = 0.0);

// Result of attempting to sell contraband via the black market.
struct BlackMarketSellResult {
  bool ok{false};
  bool stung{false};
  const char* reason{nullptr};

  econ::CommodityId commodity{econ::CommodityId::Food};
  double intendedUnits{0.0};
  double unitsSold{0.0};

  // For successful (not stung) sales.
  double pricePerUnitCr{0.0};
  double payoutCr{0.0};

  // Delta applied to credits by this attempt.
  double creditsDelta{0.0};

  // If stung, we return the underlying scan + enforcement outcome so the caller
  // can apply rep/ledger side effects and show messaging.
  IllegalCargoScanResult scan{};
  ContrabandEnforcementResult enforcement{};
};

// Sell `units` of `commodity` via a station black market.
//
// Behavior:
//  - Requires bm.available == true
//  - Requires commodity to be illegal at this station
//  - Does not modify station economy inventory (off-books)
//  - On a sting, enforces contraband on the *full* cargo under the station's
//    illegality mask (confiscation + fine schedule via PoliceScan)
//
// Pricing:
//  - Uses midPriceOverrideCr[cid] if provided; otherwise uses commodity basePrice.
//  - Converts mid -> bid using bidAskSpread/2, then applies bm.bidMul and bm.fenceCut.
//
// The caller owns reputation and law-ledger application. This function only
// mutates credits/cargo and returns the suggested penalties.
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
                                        std::array<double, econ::kCommodityCount>& cargoUnits);



// Result of attempting to buy contraband via the black market.
struct BlackMarketBuyResult {
  bool ok{false};
  bool stung{false};
  const char* reason{nullptr};

  econ::CommodityId commodity{econ::CommodityId::Food};
  double intendedUnits{0.0};
  double unitsBought{0.0};

  // For successful (not stung) purchases.
  double pricePerUnitCr{0.0};
  double costCr{0.0};

  // Delta applied to credits by this attempt.
  double creditsDelta{0.0};

  // If stung, we return the underlying scan + enforcement outcome so the caller
  // can apply rep/ledger side effects and show messaging.
  IllegalCargoScanResult scan{};
  ContrabandEnforcementResult enforcement{};
};

// Buy `units` of `commodity` via a station black market.
//
// Behavior:
//  - Requires bm.available == true
//  - Requires commodity to be illegal at this station
//  - Does not modify station economy inventory (off-books)
//  - On a sting, the player is assumed to have completed the purchase before enforcement:
//      - credits are deducted
//      - bought goods are added to cargo
//      - contraband is enforced on the full cargo under the station's illegality mask (confiscation + fine)
//
// Pricing:
//  - Uses midPriceOverrideCr[cid] if provided; otherwise uses commodity basePrice.
//  - Converts mid -> ask using bidAskSpread/2, then applies bm.askMul and bm.fenceCut.
//
// The caller owns reputation and law-ledger application. This function only
// mutates credits/cargo and returns the suggested penalties.
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
                                        std::array<double, econ::kCommodityCount>& cargoUnits);

} // namespace stellar::sim
