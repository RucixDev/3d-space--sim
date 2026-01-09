#pragma once

#include "stellar/sim/Industry.h"
#include "stellar/sim/System.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace stellar::sim {

class Universe;

// Scan "industrial trade" opportunities:
//   - buy recipe inputs at an origin station
//   - process them via an industry recipe at that origin station
//   - haul the output commodity to another station and sell it
//
// This is analogous to TradeScanner, but includes the cost/time/fees of a
// processing step.

struct IndustryTradeScanParams {
  // Universe query settings (only used by scanIndustryTradeOpportunities() overload).
  double radiusLy{200.0};
  std::size_t maxSystems{192};

  // Output shaping.
  std::size_t maxResults{12};
  std::size_t perStationLimit{1};

  // Cargo sizing (used to cap *output* hauled after processing).
  double cargoCapacityKg{420.0};
  double cargoUsedKg{0.0};
  bool useFreeHold{true};

  // Market assumptions.
  double bidAskSpread{0.10};

  // Capital limit (optional). If > 0, the scanner will cap batches so the
  // estimated upfront spend (inputs + service fee) fits within this budget.
  // This uses quoted ask prices and does not simulate price impact.
  double maxBuyCreditsCr{0.0};

  // Optional player rep at the *processing* station (affects quote yield/speed).
  double processingRep{0.0};

  // Filtering.
  double minNetProfit{0.0};
  bool includeSameSystem{true};
};

struct IndustryTradeOpportunity {
  SystemId toSystem{0};
  StationId toStation{0};

  // Optional convenience labels (filled by the scanner).
  std::string toSystemName;
  std::string toStationName;

  IndustryRecipeId recipe{IndustryRecipeId::SmeltOre};

  double batches{0.0};

  econ::CommodityId inputA{econ::CommodityId::Food};
  double inputAUnits{0.0};
  double inputAAsk{0.0};
  double inputACostCr{0.0}; // includes fees

  econ::CommodityId inputB{econ::CommodityId::Food};
  double inputBUnits{0.0};
  double inputBAsk{0.0};
  double inputBCostCr{0.0}; // includes fees

  econ::CommodityId output{econ::CommodityId::Food};
  double outputUnits{0.0};
  double outputBid{0.0};
  double outputRevenueCr{0.0}; // includes fees

  double outputMassKg{0.0};
  double serviceFeeCr{0.0};
  double timeDays{0.0};

  // Fees used.
  double feeFrom{0.0}; // applied to market buys at origin
  double feeTo{0.0};   // applied to market sells at destination

  // Final economics.
  double netProfitCr{0.0};
  double netProfitPerKg{0.0};
  double netProfitPerDay{0.0};

  // Geometry.
  double distanceLy{0.0};
};

using IndustryFeeRateFn = std::function<double(const Station&)>;

// Scan industrial-trade opportunities inside a candidate system list.
std::vector<IndustryTradeOpportunity> scanIndustryTradeOpportunities(Universe& u,
                                                                     const SystemStub& originStub,
                                                                     const Station& originStation,
                                                                     double timeDays,
                                                                     const std::vector<SystemStub>& candidates,
                                                                     const IndustryTradeScanParams& params,
                                                                     IndustryFeeRateFn feeRate = {});

// Scan industrial-trade opportunities inside a sphere around the origin system.
std::vector<IndustryTradeOpportunity> scanIndustryTradeOpportunities(Universe& u,
                                                                     const SystemStub& originStub,
                                                                     const Station& originStation,
                                                                     double timeDays,
                                                                     const IndustryTradeScanParams& params,
                                                                     IndustryFeeRateFn feeRate = {});

} // namespace stellar::sim
