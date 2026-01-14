#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/proc/HyperlaneRouter.h"
#include "stellar/proc/TradePrices.h"
#include "stellar/proc/TradeProfile.h"
#include "stellar/proc/TradeRoutes.h"

#include <vector>

namespace stellar::proc {

// -----------------------------------------------------------------------------
// Trade intel: enrich macro TradeRoutes with cheap price/profit estimates.
// -----------------------------------------------------------------------------
//
// TradeRoutes returns a *directional affinity signal* (what systems are likely
// trade partners, and what commodities fit the flow), but it intentionally does
// not talk about prices.
//
// This module adds a deterministic, profile-driven macro price field
// (TradePrices) to turn routes into actionable "intel": estimated buy/sell
// prices and rough profit numbers.
//
// IMPORTANT: This is still a macro model. It does not inspect station inventory
// or live market state.

enum class TradeIntelDirection : core::u8 {
  Export = 0, // origin -> other
  Import = 1, // other -> origin
};

struct TradeRouteEconomics {
  TradeIntelDirection direction{TradeIntelDirection::Export};

  // Mirrors proc::TradeRoute.
  sim::SystemId otherSystem{0};
  econ::CommodityId commodity{econ::CommodityId::Food};
  double score{0.0};
  double distanceLy{0.0};
  double directDistanceLy{0.0};
  double laneDistanceLy{0.0};
  int laneHops{0};
  double laneBottleneckBandwidth01{0.0};
  double risk01{0.0};

  // Macro price estimates (credits / unit).
  double buyMidCr{0.0};
  double buyAskCr{0.0};
  double sellMidCr{0.0};
  double sellBidCr{0.0};

  // Cargo scaling.
  double unitMassKg{0.0};
  double unitsForCargo{0.0};

  // Profit estimates.
  double profitPerUnitCr{0.0};
  double profitPerKgCr{0.0};
  double profitForCargoCr{0.0};
};

// Round-trip (2-leg) trade loop: origin -> other -> origin.
struct TradeLoop2 {
  sim::SystemId otherSystem{0};
  TradeRouteEconomics outLeg{};
  TradeRouteEconomics backLeg{};

  double roundTripProfitCr{0.0};
  double roundTripDistanceLy{0.0};
  double avgRisk01{0.0};
};

struct TradeIntelParams {
  // Market assumptions.
  double bidAskSpread{0.10};
  double buyFeeRate{0.0};  // fee fraction applied to buy cost (ask)
  double sellFeeRate{0.0}; // fee fraction applied to sell revenue (bid)

  // Cargo sizing for profit/trip estimates.
  double cargoCapacityKg{420.0};

  // Loop output shaping.
  std::size_t maxLoops{8};
};

struct TradeIntelReport {
  std::vector<TradeRouteEconomics> exports;
  std::vector<TradeRouteEconomics> imports;
  std::vector<TradeLoop2> loops;
};

// Build a small macro trade report for a selected system.
//
// Requirements:
//  - candidates.size() == candidateProfiles.size()
//
// Notes:
//  - The origin stub may appear in candidates; it is ignored.
//  - Output ordering for exports/imports matches TradeRoutes ordering.
TradeIntelReport buildTradeIntel(const sim::SystemStub& origin,
                                 const TradeProfile& originProfile,
                                 const std::vector<sim::SystemStub>& candidates,
                                 const std::vector<TradeProfile>& candidateProfiles,
                                 const TradeRouteParams& routeParams = {},
                                 const TradeIntelParams& intelParams = {},
                                 const TradePriceModelParams& priceParams = {});

// Hyperlane-aware overload: routes are computed using the supplied hyperlane network.
TradeIntelReport buildTradeIntel(const sim::SystemStub& origin,
                                 const TradeProfile& originProfile,
                                 const std::vector<sim::SystemStub>& candidates,
                                 const std::vector<TradeProfile>& candidateProfiles,
                                 const HyperlaneNetwork& hyperlanes,
                                 const TradeRouteParams& routeParams = {},
                                 const HyperlaneTravelParams& travelParams = {},
                                 const TradeIntelParams& intelParams = {},
                                 const TradePriceModelParams& priceParams = {});

} // namespace stellar::proc
