#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/System.h"
#include "stellar/sim/SystemEvents.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace stellar::sim {

class Universe;

// Smuggling scanner / suggestion helper.
//
// This module bridges:
//  - deterministic contraband legality masks (Contraband.h)
//  - deterministic station black market profiles (BlackMarket.h)
//
// ...into a UI-agnostic list of profitable *illegal* trade opportunities.
//
// Conceptual model:
//  - Buy a commodity legally at the origin station (official market ask).
//  - Sell it at the destination station via the black market (black-market bid).
//  - A "sting" event may occur with probability p, resulting in confiscation and a fine.
//
// The scanner outputs both the clean profit and a risk-aware score.

enum class SmugglingAvailabilityMode : core::u8 {
  // Only return opportunities where the fence is available today.
  TodayOnly = 0,

  // Treat access01 as a probability of "today" availability and scale the score by it.
  // (This is useful for planning where you're willing to wait or re-roll days.)
  Expected = 1,

  // Ignore daily availability (treat as always available).
  // Useful for tooling/unit tests.
  Ignore = 2,
};

enum class SmugglingScoreMode : core::u8 {
  // Expected profit (mean) after sting probability.
  ExpectedProfit = 0,

  // Mean-variance style: expectedProfit - riskLambda * stdDev.
  RiskAdjusted = 1,

  // Profit assuming no sting (pure upside).
  CleanProfit = 2,

  // Expected profit per light-year traveled.
  ExpectedProfitPerLy = 3,
};

using SmugglingRepFn = std::function<double(core::u32 /*factionId*/) >;

struct SmugglingScanParams {
  // Universe query settings (only used by scanSmugglingOpportunities()).
  double radiusLy{200.0};
  std::size_t maxSystems{192};

  // Output shaping.
  std::size_t maxResults{12};
  std::size_t perStationLimit{2};

  // Cargo sizing.
  double cargoCapacityKg{420.0};
  double cargoUsedKg{0.0};
  bool useFreeHold{true};

  // Market assumptions.
  double bidAskSpread{0.10};

  // Behavior.
  bool includeSameSystem{true};
  bool requireOriginLegal{true};

  // Player context.
  //
  // `playerRep` is used as a fallback value; if `repForFaction` is provided it will
  // be queried per-destination station using that station's factionId.
  double playerRep{0.0};
  SmugglingRepFn repForFaction{};
  double playerHeat{0.0};
  int smuggleHoldMk{0};

  // If true, incorporate the Universe's live system security dynamics and system events
  // (if configured) when modelling black market access/risk.
  //
  // This does *not* change legality masks (those remain deterministic).
  bool useLiveSystemConditions{false};

  // Scoring.
  SmugglingAvailabilityMode availability{SmugglingAvailabilityMode::TodayOnly};
  SmugglingScoreMode scoreMode{SmugglingScoreMode::RiskAdjusted};

  // Risk penalty factor used when scoreMode==RiskAdjusted.
  //  - 0 => no penalty (reduces to ExpectedProfit).
  //  - 1 => subtract one standard deviation.
  double riskLambda{0.65};

  // Filtering.
  double minScoreCr{0.0};
};

struct SmugglingOpportunity {
  SystemId toSystem{0};
  StationId toStation{0};

  std::string toSystemName;
  std::string toStationName;

  econ::CommodityId commodity{econ::CommodityId::Food};

  // Prices.
  double buyAskCr{0.0};          // official ask at origin
  double buyAskNetCr{0.0};       // ask with origin fee
  double sellBidCr{0.0};         // black market bid at destination (includes fence cut)
  double officialSellBidCr{0.0}; // official destination bid (for debugging/UI)

  // Feasibility.
  double unitsFrom{0.0};
  double unitsToSpace{0.0};
  double unitsPossible{0.0};
  double unitMassKg{0.0};

  // Geometry.
  double distanceLy{0.0};

  // System conditions context (filled for every row; systemEventKind is only meaningful
  // when params.useLiveSystemConditions==true).
  SystemEventKind systemEventKind{SystemEventKind::None};
  double systemEventSeverity01{0.0};
  double systemSecurity01{0.0};
  double systemPiracy01{0.0};
  double systemTraffic01{0.0};

  // Black market profile summary.
  bool blackMarketAvailable{false};
  double blackMarketAccess01{0.0};

  // Enforcement/risk model.
  double stingChance{0.0};
  double illegalValueCr{0.0};
  double fineCr{0.0};

  // Outcome metrics (for the planned cargo size).
  double buyCostCr{0.0};
  double payoutCr{0.0};

  double cleanProfitCr{0.0};    // profit if not stung
  double stungProfitCr{0.0};    // profit if stung (confiscation + fine)
  double expectedProfitCr{0.0}; // mean
  double stdDevCr{0.0};         // sqrt(var)

  // Final score used for ranking (includes availability handling).
  double scoreCr{0.0};
};

using SmugglingFeeRateFn = std::function<double(const Station&)>;

// Scan smuggling opportunities inside a sphere around the origin system.
std::vector<SmugglingOpportunity> scanSmugglingOpportunities(Universe& u,
                                                            const SystemStub& originStub,
                                                            const Station& originStation,
                                                            double timeDays,
                                                            const SmugglingScanParams& params,
                                                            SmugglingFeeRateFn feeRate = {});

// Scan smuggling opportunities using a pre-selected candidate system list.
std::vector<SmugglingOpportunity> scanSmugglingOpportunities(Universe& u,
                                                            const SystemStub& originStub,
                                                            const Station& originStation,
                                                            double timeDays,
                                                            const std::vector<SystemStub>& candidates,
                                                            const SmugglingScanParams& params,
                                                            SmugglingFeeRateFn feeRate = {});

} // namespace stellar::sim
