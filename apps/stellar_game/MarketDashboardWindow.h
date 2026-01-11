#pragma once

#include "stellar/econ/Commodity.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Universe.h"

#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::game {

// A row in the market dashboard table.
struct MarketDashboardRow {
  stellar::sim::SystemId systemId{0};
  stellar::sim::StationId stationId{0};

  std::string systemName;
  std::string stationName;
  stellar::econ::StationType stationType{stellar::econ::StationType::Outpost};

  double distanceLy{0.0};

  // Prices are in credits per unit.
  double mid{0.0};
  double askEff{0.0}; // effective ask incl. station fees
  double bidEff{0.0}; // effective bid incl. station fees

  double inventory{0.0};
  double capacity{0.0};
  double feeRate{0.0};

  // Price trend computed over a recent window (see MarketDashboardWindowState::trendWindowDays).
  // This is a regression slope (credits/unit/day), not a single-step derivative.
  double trendPerDay{0.0};
  double trendR2{0.0}; // goodness-of-fit for regression (0..1)
  double changePct{0.0}; // percent change over the trend window
  double volatilityPerDay{0.0}; // stddev of log returns per day (e.g., 0.02 ~= 2%/day)
  double lastSampleDay{0.0};
  bool trendValid{false};
};

struct MarketDashboardWindowState {
  bool open{false};

  int commodityIndex{0};

  // Scope: current system, or scan nearby systems.
  bool scopeNearby{true};
  double radiusLy{350.0};
  int maxSystems{256};
  bool includeCurrentSystem{true};

  // 0=Best Buy (ask asc), 1=Best Sell (bid desc), 2=Distance, 3=Name
  int sortMode{0};

  // Analytics tuning.
  double trendWindowDays{7.0};    // window used to compute trend in table rows
  double historyWindowDays{30.0}; // window used for the selected station's plot/stats

  bool showHistoryPanel{true};
  char filter[96]{};

  // Selection
  stellar::sim::SystemId selectedSystem{0};
  stellar::sim::StationId selectedStation{0};

  // Cache
  int cacheStamp{-999999};
  int cacheCommodity{-1};
  bool cacheScopeNearby{false};
  double cacheRadiusLy{0.0};
  int cacheMaxSystems{0};
  bool cacheIncludeCurrentSystem{true};
  double cacheTrendWindowDays{0.0};

  std::vector<MarketDashboardRow> rows;
};

// Runtime context the dashboard needs from the game.
struct MarketDashboardContext {
  stellar::sim::Universe& universe;
  const stellar::sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};

  // Player context (optional; used for profit estimates).
  double cargoCapacityKg{0.0};
  double cargoUsedKg{0.0};
  double playerCreditsCr{0.0};

  // Fee model used for quotes.
  std::function<double(const stellar::sim::Station&)> effectiveFeeRate;

  // UI callbacks.
  std::function<void(stellar::sim::SystemId, stellar::sim::StationId)> routeToStation;
  std::function<void(stellar::sim::StationId)> targetStation; // only meaningful for current system
  std::function<void(std::string_view, double)> toast;
};

// Draw the market dashboard window. (No-op if st.open == false.)
void drawMarketDashboardWindow(MarketDashboardWindowState& st, const MarketDashboardContext& ctx);

} // namespace stellar::game
