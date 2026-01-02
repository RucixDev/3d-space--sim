#pragma once

#include "stellar/econ/RoutePlanner.h"
#include "stellar/sim/System.h"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace stellar::core { class JobSystem; }

namespace stellar::sim {

class Universe;

// Fee-rate callback for trade scanning.
// Return a station fee rate in [0..1]. Implementations should be robust to out-of-range values.
using TradeLoopFeeRateFn = std::function<double(const Station&)>;

// Parameters for multi-leg trade-loop scanning.
//
// The scanner looks for profitable closed loops that start and end at the origin station.
// Currently supported loop lengths:
//  - 2 legs: origin -> B -> origin
//  - 3 legs: origin -> B -> C -> origin
//
// Each leg uses the cargo-manifest planner from `stellar/econ/RoutePlanner.h`.
// NOTE: The scanner does not mutate station economy states.
struct TradeLoopScanParams {
  // Manifest planner settings used for each leg (fees are overridden per-leg via TradeLoopFeeRateFn).
  econ::CargoManifestParams manifest{};

  // Filtering
  double minLegProfitCr{0.0};   // minimum net profit per leg (after fees)
  double minLoopProfitCr{0.0};  // minimum total net profit for the whole loop

  // Whether to include other stations in the origin system as candidates.
  bool includeSameSystem{true};

  // Search / limits
  std::size_t legs{2};              // 2 or 3
  std::size_t maxLegCandidates{16}; // how many best outgoing legs to expand at each step
  std::size_t maxResults{12};       // how many loops to return
  std::size_t maxStations{256};     // cap on candidate station nodes (across all candidate systems)
};

struct TradeLoopLeg {
  SystemId fromSystem{0};
  StationId fromStation{0};
  SystemId toSystem{0};
  StationId toStation{0};

  std::string fromSystemName;
  std::string fromStationName;
  std::string toSystemName;
  std::string toStationName;

  double distanceLy{0.0}; // system-to-system distance (in-system distance is ignored)
  double feeFrom{0.0};
  double feeTo{0.0};

  // Leg manifest and net profit (after fees).
  econ::CargoManifestPlan manifest{};
};

struct TradeLoop {
  std::vector<TradeLoopLeg> legs;

  double totalProfitCr{0.0};
  double totalDistanceLy{0.0};
  double profitPerLy{0.0};
};

// Scan trade loops using a pre-selected candidate system list (typically from Universe::queryNearby()).
std::vector<TradeLoop> scanTradeLoops(Universe& u,
                                      const SystemStub& originStub,
                                      const Station& originStation,
                                      double timeDays,
                                      const std::vector<SystemStub>& candidates,
                                      const TradeLoopScanParams& params,
                                      TradeLoopFeeRateFn feeRate = {});

// Convenience overload: performs a queryNearby() centered at originStub.posLy.
std::vector<TradeLoop> scanTradeLoops(Universe& u,
                                      const SystemStub& originStub,
                                      const Station& originStation,
                                      double timeDays,
                                      double radiusLy,
                                      std::size_t maxSystems,
                                      const TradeLoopScanParams& params,
                                      TradeLoopFeeRateFn feeRate = {});


// Parallel variant of scanTradeLoops() using a caller-provided JobSystem.
//
// This pre-resolves station data and snapshots station economy states on the
// calling thread (to avoid thread-safety issues with Universe caches), then
// computes the heavy cargo-manifest legs in parallel.
std::vector<TradeLoop> scanTradeLoopsParallel(core::JobSystem& jobs,
                                              Universe& u,
                                              const SystemStub& originStub,
                                              const Station& originStation,
                                              double timeDays,
                                              const std::vector<SystemStub>& candidates,
                                              const TradeLoopScanParams& params,
                                              TradeLoopFeeRateFn feeRate = {});

// Convenience overload: performs a queryNearbyParallel() centered at originStub.posLy.
std::vector<TradeLoop> scanTradeLoopsParallel(core::JobSystem& jobs,
                                              Universe& u,
                                              const SystemStub& originStub,
                                              const Station& originStation,
                                              double timeDays,
                                              double radiusLy,
                                              std::size_t maxSystems,
                                              const TradeLoopScanParams& params,
                                              TradeLoopFeeRateFn feeRate = {});

} // namespace stellar::sim
