#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/RoutePlanner.h"
#include "stellar/sim/System.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace stellar::sim {

class Universe;

// Fee-rate callback for trade scanning.
// Return a station fee rate in [0..1]. Implementations should be robust to out-of-range values.
using TradeRunFeeRateFn = std::function<double(const Station&)>;

// How to rank trade runs.
//
// TotalProfit is simplest and tends to produce intuitive "best raw money" chains.
// The ratio modes are helpful when travel constraints are enabled and you want
// more efficiency-oriented suggestions.
enum class TradeRunScoreMode : core::u8 {
  TotalProfit  = 0,
  ProfitPerLy  = 1,
  ProfitPerHop = 2,
  ProfitPerCost = 3,
};

// Parameters for multi-leg trade-run scanning.
//
// A "trade run" is a sequence of legs:
//   origin -> A -> B -> ...
//
// Each leg uses the cargo-manifest planner from `stellar/econ/RoutePlanner.h`.
// NOTE: The scanner does not mutate station economy states.
struct TradeRunScanParams {
  // Manifest planner settings used for each leg (fees are overridden per-leg via TradeRunFeeRateFn).
  econ::CargoManifestParams manifest{};

  // Filtering.
  double minLegProfitCr{0.0}; // minimum net profit per leg (after fees)
  double minRunProfitCr{0.0}; // minimum net profit for the entire run

  // Whether to include other stations in the origin system as candidates.
  bool includeSameSystem{true};

  // If true, a run cannot visit the same station twice.
  // (This avoids silly A->B->A->B oscillation runs when scanning with 3+ legs.)
  bool loopless{true};

  // Search / limits.
  std::size_t legs{2};              // number of legs (>= 1)
  std::size_t beamWidth{32};        // beam-search width (controls CPU vs quality)
  std::size_t maxLegCandidates{16}; // how many best outgoing legs to consider per stop
  std::size_t maxResults{12};       // how many runs to return
  std::size_t maxStations{256};     // cap on candidate station nodes (across all candidate systems)

  // Navigation / reachability.
  //
  // If jumpRangeLy > 0, each leg must be reachable under that maximum single-jump
  // distance using the provided candidate system graph. Each leg stores the computed
  // hop route in TradeRunLeg::route.
  //
  // If jumpRangeLy <= 0, reachability is not checked and legs assume a direct jump.
  double jumpRangeLy{0.0};

  // Route cost model (used for ProfitPerCost ranking and per-leg bookkeeping).
  // The default is a hop-only "fuel-like" cost where each jump costs 1 unit.
  double routeCostPerJump{1.0};
  double routeCostPerLy{0.0};

  TradeRunScoreMode scoreMode{TradeRunScoreMode::TotalProfit};
};

struct TradeRunLeg {
  SystemId fromSystem{0};
  StationId fromStation{0};
  SystemId toSystem{0};
  StationId toStation{0};

  std::string fromSystemName;
  std::string fromStationName;
  std::string toSystemName;
  std::string toStationName;

  // Jump route between fromSystem and toSystem (includes start + goal).
  std::vector<SystemId> route{};
  int routeHops{0};
  double routeDistanceLy{0.0};
  double routeCost{0.0};

  double feeFrom{0.0};
  double feeTo{0.0};

  // Leg manifest and net profit (after fees).
  econ::CargoManifestPlan manifest{};
};

struct TradeRun {
  std::vector<TradeRunLeg> legs{};

  double totalProfitCr{0.0};
  double totalRouteDistanceLy{0.0};
  double totalRouteCost{0.0};
  int totalHops{0};

  // Convenience precomputed ratios.
  double profitPerLy{0.0};
  double profitPerHop{0.0};
  double profitPerCost{0.0};
};

// Plan profitable multi-leg trade runs using a pre-selected candidate system list
// (typically from Universe::queryNearby()).
std::vector<TradeRun> planTradeRuns(Universe& u,
                                    const SystemStub& originStub,
                                    const Station& originStation,
                                    double timeDays,
                                    const std::vector<SystemStub>& candidates,
                                    const TradeRunScanParams& params,
                                    TradeRunFeeRateFn feeRate = {});

// Convenience overload: performs a queryNearby() centered at originStub.posLy.
std::vector<TradeRun> planTradeRuns(Universe& u,
                                    const SystemStub& originStub,
                                    const Station& originStation,
                                    double timeDays,
                                    double radiusLy,
                                    std::size_t maxSystems,
                                    const TradeRunScanParams& params,
                                    TradeRunFeeRateFn feeRate = {});

} // namespace stellar::sim
