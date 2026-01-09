#pragma once

#include "stellar/core/JobSystem.h"
#include "stellar/core/Types.h"
#include "stellar/sim/System.h"
#include "stellar/sim/IndustryScanner.h"
#include "stellar/sim/TradeLoopScanner.h"
#include "stellar/sim/TradeRunPlanner.h"
#include "stellar/sim/Universe.h"

#include <functional>
#include <memory>
#include <string_view>
#include <vector>

namespace stellar::game {

// Multi-leg trade route planner UI.
//
// This window exposes the headless trade search algorithms (TradeRunPlanner
// and TradeLoopScanner) to the in-game UI.

struct TradePlannerWindowState {
  bool open{false};

  enum class Mode : int {
    Runs = 0,
    Loops = 1,
    Industry = 2,
  };

  // Window mode (multi-leg runs vs loops).
  Mode mode{Mode::Runs};

  // Query scope.
  double radiusLy{350.0};
  int maxSystems{256};
  int maxStations{256};
  bool includeSameSystem{true};

  // Cargo / market assumptions.
  bool useFreeHold{true};
  bool limitByCredits{true};
  double bidAskSpread{0.10};
  double stepKg{1.0};
  bool simulatePriceImpact{true};

  // Profit filters.
  double minLegProfitCr{0.0};
  double minTotalProfitCr{0.0};
  int maxResults{12};

  // Runs settings.
  int runLegs{2};
  int beamWidth{32};
  int maxLegCandidates{16};
  bool loopless{true};

  bool enforceJumpRange{false};
  double jumpRangeLy{15.0};
  double routeCostPerJump{1.0};
  double routeCostPerLy{0.0};
  int scoreMode{0};

  // Loops settings.
  int loopLegs{2};

  // Industry settings.
  int industryPerStationLimit{1};

  // Performance.
  bool useParallel{true};
  int threads{0};      // 0 = auto
  bool autoRefresh{false};

  // Export.
  char exportPath[256]{"trade_planner.csv"};

  // --- Internal state / caches ---
  stellar::sim::SystemId cacheSystemId{0};
  stellar::sim::StationId fromStationId{0};
  int cacheStamp{-999999};
  double cacheCargoCapKg{0.0};
  double cacheCargoUsedKg{0.0};
  double cacheCreditsCr{0.0};
  Mode cacheMode{Mode::Runs};

  double lastComputeMs{0.0};

  std::vector<stellar::sim::TradeRun> runs;
  std::vector<stellar::sim::TradeLoop> loops;
  std::vector<stellar::sim::IndustryTradeOpportunity> industryOps;

  std::unique_ptr<stellar::core::JobSystem> jobs;
  std::size_t jobsThreadCount{0};
};

struct TradePlannerContext {
  stellar::sim::Universe& universe;
  const stellar::sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};

  // Suggested default origin station (docked/targeted), if any.
  stellar::sim::StationId preferredFromStationId{0};

  // Player context.
  double cargoCapacityKg{0.0};
  double cargoUsedKg{0.0};
  double playerCreditsCr{0.0};

  // Fee model used to produce realistic net profits.
  std::function<double(const stellar::sim::Station&)> effectiveFeeRate;

  // Optional: faction reputation lookup (used by industry quotes).
  std::function<double(stellar::core::u32)> reputationForFaction;

  // UI callbacks.
  std::function<void(stellar::sim::SystemId, stellar::sim::StationId)> routeToStation;
  std::function<void(std::string_view, double)> toast;
};

void drawTradePlannerWindow(TradePlannerWindowState& st, const TradePlannerContext& ctx);

} // namespace stellar::game
