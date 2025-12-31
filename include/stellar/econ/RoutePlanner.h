#pragma once

#include "stellar/econ/Market.h"

#include <vector>

namespace stellar::econ {

struct RouteOpportunity {
  CommodityId commodity{};

  // Raw (no fees): toBid - fromAsk.
  double profitPerUnit{0.0};

  // Raw prices (no fees):
  //  - buyPrice:  source station ask (station sells to player)
  //  - sellPrice: destination station bid (station buys from player)
  double buyPrice{0.0};
  double sellPrice{0.0};

  // Feasibility / sizing helpers.
  // These are populated by the route planner so UI/gameplay can avoid
  // recommending trades that the stations (or cargo hold) cannot support.
  double unitsFrom{0.0};      // available at source station (inventory units)
  double unitsToSpace{0.0};   // free capacity at destination station (units)
  double unitsPossible{0.0};  // min(unitsFrom, unitsToSpace[, cargo cap])

  // Commodity sizing.
  double unitMassKg{0.0};
  double profitPerKg{0.0};

  // Total raw profit for `unitsPossible` (no fees).
  double profitTotal{0.0};

  // Optional fees (only filled by bestRoutesForCargo).
  double feeFrom{0.0};
  double feeTo{0.0};
  double netProfitPerUnit{0.0};
  double netProfitTotal{0.0};
};

std::vector<RouteOpportunity> bestRoutes(const StationEconomyState& fromState,
                                        const StationEconomyModel& fromModel,
                                        const StationEconomyState& toState,
                                        const StationEconomyModel& toModel,
                                        double bidAskSpread = 0.10,
                                        std::size_t maxResults = 5);

// Cargo-aware route planner.
//
// Computes feasible trade opportunities given:
//  - station inventories/capacities
//  - an optional cargo mass limit (kg)
//  - optional station fee rates (0..1)
//
// Results are sorted by netProfitTotal (trip profit) when fees/cargo are provided.
std::vector<RouteOpportunity> bestRoutesForCargo(const StationEconomyState& fromState,
                                                const StationEconomyModel& fromModel,
                                                const StationEconomyState& toState,
                                                const StationEconomyModel& toModel,
                                                double cargoCapacityKg,
                                                double fromFeeRate = 0.0,
                                                double toFeeRate = 0.0,
                                                double bidAskSpread = 0.10,
                                                std::size_t maxResults = 5);



// ------------------------------
// Cargo manifest planner (multi-commodity trade mix)
//
// Unlike bestRoutesForCargo(), which ranks *single commodity* routes, the
// manifest planner can recommend a *mix* of commodities that jointly maximizes
// trip profit under a cargo mass limit.
//
// The default planner is a greedy "marginal profit per kg" allocator that can
// optionally simulate price impact by updating station inventories as it fills
// the hold. This approximates the best mix when prices respond to inventory.

struct CargoManifestParams {
  // Cargo constraint
  double cargoCapacityKg{0.0};

  // Market assumptions
  double bidAskSpread{0.10};

  // Station fees applied on buy/sell (0..1).
  double fromFeeRate{0.0};
  double toFeeRate{0.0};

  // Optimization knobs
  double stepKg{1.0};              // resolution of the greedy allocator
  double maxBuyCreditsCr{0.0};     // <= 0 => ignore player credits
  bool simulatePriceImpact{true};  // if true, inventory changes affect prices while planning
};

struct CargoManifestLine {
  CommodityId commodity{};

  double units{0.0};
  double unitMassKg{0.0};
  double massKg{0.0};

  // Net prices (after station fees), averaged over the planned volume.
  double avgNetBuyPrice{0.0};
  double avgNetSellPrice{0.0};

  // Net totals (after station fees).
  double netBuyCr{0.0};
  double netSellCr{0.0};
  double netProfitCr{0.0};

  double netProfitPerUnit{0.0};
  double netProfitPerKg{0.0};
};

struct CargoManifestPlan {
  double cargoCapacityKg{0.0};
  double cargoFilledKg{0.0};

  // Totals are net-of-fees.
  double netBuyCr{0.0};
  double netSellCr{0.0};
  double netProfitCr{0.0};

  std::vector<CargoManifestLine> lines;
};

// Compute a best-effort profitable cargo manifest for a single trip.
//
// Notes:
//  - The planner does not mutate the input economy states.
//  - If params.simulatePriceImpact is enabled, the planner uses internal copies
//    of the economy states and updates inventories as it fills the hold so
//    prices "move" with volume.
//  - If no profitable trade exists, the returned plan will have 0 profit and an
//    empty lines list.
CargoManifestPlan bestManifestForCargo(const StationEconomyState& fromState,
                                       const StationEconomyModel& fromModel,
                                       const StationEconomyState& toState,
                                       const StationEconomyModel& toModel,
                                       const CargoManifestParams& params);


} // namespace stellar::econ
