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
#include <string>
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

  // ---------------------------------------------------------------------------
  // Trade Route Runner (cross-integrates TradePlanner -> Nav/Docking/Comms)
  // ---------------------------------------------------------------------------
  //
  // A lightweight, UI-owned plan that can be executed by the main game loop
  // (auto-run to the next stop, auto-dock when in range). This intentionally
  // stores only IDs + UI-friendly labels + the precomputed manifest.
  //
  // Notes:
  //  - We copy legs when arming a route so the runner stays stable even if the
  //    user re-scans trade ideas.
  //  - Persistence is intentionally left to the save system; the runner is
  //    designed to be safe to drop on load.

  struct TradeRouteLeg {
    stellar::sim::SystemId fromSystem{0};
    stellar::sim::StationId fromStation{0};
    stellar::sim::SystemId toSystem{0};
    stellar::sim::StationId toStation{0};

    std::string fromSystemName;
    std::string fromStationName;
    std::string toSystemName;
    std::string toStationName;

    // Geometry / scoring hints.
    double distanceLy{0.0};
    int hops{0};

    // Per-leg economics (net-of-fees).
    double netProfitCr{0.0};
    stellar::econ::CargoManifestPlan manifest{};
  };

  struct TradeRouteRunner {
    bool active{false};

    // Route type: loops can optionally repeat.
    bool isLoop{false};
    bool repeat{false};

    // If true, the runner will arm nav auto-run on each leg transition.
    bool armAutoRun{true};

    // If true, execute the leg's manifest automatically when docking at each stop:
    //  - Sell inbound cargo from the just-completed leg.
    //  - Buy outbound cargo for the next leg.
    // This mirrors the tracked trade-loop automation and makes multi-leg runs feel
    // like a coherent "dockside turn".
    bool autoTradeOnDock{true};
    bool autoTradeToast{true};
    bool autoTradeAllowIllegalViaBlackMarket{true};

    // If true, automatically plot the next leg when arriving at a stop.
    // (When armAutoRun is also true, this will arm auto-run/auto-dock.)
    bool autoPlotNextLeg{true};

    // UI -> main loop bridge: request a one-shot dockside execution.
    bool requestExecuteDockedTurn{false};

    // Index of the current leg we are traveling *towards* (destination).
    int legIndex{0};

    // Integration markers.
    bool startEventEmitted{false};
    bool endEventPending{false};
    enum class EndReason : int { None = 0, Cancelled = 1 };
    EndReason endReason{EndReason::None};

    double startedAtDays{0.0};
    double lastAdvanceAtDays{0.0};

    // Last auto-trade summary (for UI feedback).
    double lastTradeTimeDays{0.0};
    stellar::sim::StationId lastTradeStationId{0};
    int lastTradeSoldLines{0};
    int lastTradeBoughtLines{0};
    int lastTradeIllegalSkips{0};
    double lastTradeCreditsDelta{0.0};
    double lastTradeSoldCr{0.0};
    double lastTradeBoughtCr{0.0};

    std::vector<TradeRouteLeg> legs;
  };

  TradeRouteRunner route{};

  std::unique_ptr<stellar::core::JobSystem> jobs;
  std::size_t jobsThreadCount{0};
};

struct TradePlannerContext {
  stellar::sim::Universe& universe;
  const stellar::sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};

  // Docking context (enables route runner UI + optional auto-trade).
  bool docked{false};
  stellar::sim::StationId dockedStationId{0};

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
  // Like routeToStation, but can also arm nav auto-run + auto-docking.
  // If unset, the UI will gracefully fall back to routeToStation.
  std::function<void(stellar::sim::SystemId, stellar::sim::StationId, bool armAutoRun)> goToStation;
  std::function<void(std::string_view, double)> toast;
};

void drawTradePlannerWindow(TradePlannerWindowState& st, const TradePlannerContext& ctx);

} // namespace stellar::game
