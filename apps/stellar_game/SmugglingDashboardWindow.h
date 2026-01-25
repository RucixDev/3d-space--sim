#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/SmugglingScanner.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Universe.h"

#include <functional>
#include <string_view>
#include <vector>

namespace stellar::game {

// A UI window that scans nearby space for profitable smuggling opportunities.
// This is a thin wrapper around sim::scanSmugglingOpportunities(...), with
// enough player context (cargo, heat, smuggling hold upgrades, faction rep) to
// produce risk-adjusted, actionable results.
struct SmugglingDashboardWindowState {
  bool open{false};

  // Origin station selection (within the current system).
  stellar::sim::StationId originStationId{0};

  // Query scope.
  double radiusLy{300.0};
  int maxSystems{256};
  bool includeSameSystem{true};

  // Cargo / market assumptions.
  bool useFreeHold{true};
  double bidAskSpread{0.10};

  // Smuggling constraints.
  bool requireOriginLegal{true};
  bool useLiveSystemConditions{true}; // include security dynamics + system events in BM modelling

  // Scoring / filtering.
  stellar::sim::SmugglingAvailabilityMode availability{stellar::sim::SmugglingAvailabilityMode::Expected};
  stellar::sim::SmugglingScoreMode scoreMode{stellar::sim::SmugglingScoreMode::RiskAdjusted};
  double riskLambda{0.65};
  double minScoreCr{0.0};

  int maxResults{24};
  int perStationLimit{2};

  // Performance.
  bool autoRefreshDaily{false};

  // UI.
  bool showAdvanced{false};
  char filter[96]{};

  // --- Cache ---
  int cacheDay{-999999};
  stellar::sim::StationId cacheOriginStationId{0};
  double cacheRadiusLy{0.0};
  int cacheMaxSystems{0};
  bool cacheIncludeSameSystem{true};
  bool cacheUseFreeHold{true};
  double cacheBidAskSpread{0.0};
  bool cacheRequireOriginLegal{true};
  bool cacheUseLiveSystemConditions{false};
  stellar::sim::SmugglingAvailabilityMode cacheAvailability{stellar::sim::SmugglingAvailabilityMode::Expected};
  stellar::sim::SmugglingScoreMode cacheScoreMode{stellar::sim::SmugglingScoreMode::RiskAdjusted};
  double cacheRiskLambda{0.0};
  double cacheMinScoreCr{0.0};
  int cacheMaxResults{0};
  int cachePerStationLimit{0};

  double lastComputeMs{0.0};
  std::vector<stellar::sim::SmugglingOpportunity> results;
};

struct SmugglingDashboardContext {
  stellar::sim::Universe& universe;
  const stellar::sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};

  // Suggested default origin station (docked/targeted), if any.
  stellar::sim::StationId preferredOriginStationId{0};

  // Player context.
  double cargoCapacityKg{0.0};
  double cargoUsedKg{0.0};
  double playerCreditsCr{0.0};
  double playerHeat{0.0};
  int smuggleHoldMk{0};

  // Fee model used for origin buy quotes.
  std::function<double(const stellar::sim::Station&)> effectiveFeeRate;

  // Per-faction reputation lookup (optional but recommended).
  std::function<double(stellar::core::u32)> reputationForFaction;

  // UI callbacks.
  std::function<void(stellar::sim::SystemId, stellar::sim::StationId)> routeToStation;
  std::function<void(stellar::sim::SystemId, stellar::sim::StationId, bool armAutoRun)> goToStation;
  std::function<void(std::string_view, double)> toast;
};

void drawSmugglingDashboardWindow(SmugglingDashboardWindowState& st, const SmugglingDashboardContext& ctx);

} // namespace stellar::game
