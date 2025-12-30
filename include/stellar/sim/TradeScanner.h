#pragma once

#include "stellar/econ/Commodity.h"
#include "stellar/sim/System.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace stellar::sim {

class Universe;

// Trade scanner / route suggestion helper.
//
// This module is intentionally UI-agnostic: it computes a ranked list of
// profitable, cargo-feasible trade opportunities between stations.
//
// It is designed to be shared by:
//  - the in-game Trade Helper window
//  - headless tooling (stellar_sandbox --trade)
//
// The scanner delegates fee computation to a caller-provided function so
// reputation / tariffs can be applied differently by each front-end.

struct TradeScanParams {
  // Universe query settings (only used by scanTradeOpportunities()).
  double radiusLy{200.0};
  std::size_t maxSystems{192};

  // Output shaping.
  std::size_t maxResults{12};
  std::size_t perStationLimit{1}; // number of top opportunities to consider per destination station

  // Cargo sizing.
  double cargoCapacityKg{420.0};
  double cargoUsedKg{0.0};
  bool useFreeHold{true}; // if true, effective capacity = max(0, cargoCapacityKg - cargoUsedKg)

  // Market assumptions.
  double bidAskSpread{0.10};

  // Filtering.
  double minNetProfit{0.0}; // min net profit per trip (after fees)
  bool includeSameSystem{true};

  bool commodityFilterEnabled{false};
  econ::CommodityId commodityFilter{econ::CommodityId::Food};
};

struct TradeOpportunity {
  SystemId toSystem{0};
  StationId toStation{0};

  // Optional convenience labels (filled by the scanner).
  std::string toSystemName;
  std::string toStationName;

  econ::CommodityId commodity{econ::CommodityId::Food};

  // Raw prices (no fees).
  double buyPrice{0.0};
  double sellPrice{0.0};

  // Feasibility.
  double unitsFrom{0.0};
  double unitsToSpace{0.0};
  double unitsPossible{0.0};
  double unitMassKg{0.0};

  // Fees applied.
  double feeFrom{0.0};
  double feeTo{0.0};

  // Net profits (after fees).
  double netProfitPerUnit{0.0};
  double netProfitTotal{0.0};

  // Geometry.
  double distanceLy{0.0};
};

using TradeFeeRateFn = std::function<double(const Station&)>;

// Scan trade opportunities inside a sphere around the origin system.
//
// `feeRate` may be empty; if so, Station::feeRate is used.
std::vector<TradeOpportunity> scanTradeOpportunities(Universe& u,
                                                     const SystemStub& originStub,
                                                     const Station& originStation,
                                                     double timeDays,
                                                     const TradeScanParams& params,
                                                     TradeFeeRateFn feeRate = {});

// Scan trade opportunities using a pre-selected candidate system list.
//
// This is useful for tooling/tests that already performed a queryNearby() call.
std::vector<TradeOpportunity> scanTradeOpportunities(Universe& u,
                                                     const SystemStub& originStub,
                                                     const Station& originStation,
                                                     double timeDays,
                                                     const std::vector<SystemStub>& candidates,
                                                     const TradeScanParams& params,
                                                     TradeFeeRateFn feeRate = {});

} // namespace stellar::sim
