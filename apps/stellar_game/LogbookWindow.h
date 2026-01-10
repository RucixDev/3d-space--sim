#pragma once

#include "stellar/sim/Logbook.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Universe.h"

#include <functional>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace stellar::game {

struct ExplorationDataBrokerState {
  // Premium model (optional): payout multiplier based on distance and jurisdiction.
  bool enablePremium{true};
  float distanceScaleLy{250.0f};
  float maxDistancePremium{0.30f}; // up to +30% for very distant data
  float sameFactionBonus{0.10f};   // +10% when selling data from same faction space
  float otherFactionPenalty{0.05f}; // -5% when selling foreign-space data

  // UI
  bool showDetails{false};
  bool groupBySystem{true};

  // Selection for partial sales.
  std::unordered_set<stellar::sim::SystemId> selectedSystems;
};

struct LogbookWindowState {
  bool open{false};

  bool showUnsold{true};
  bool showSold{true};

  // -1 = all
  int kindFilter{-1};

  // 0=newest first, 1=highest value first, 2=system name
  int sortMode{0};

  char filter[96]{0};

  // Selection.
  stellar::core::u64 selectedKey{0};

  // Embedded broker state (also used by station UI).
  ExplorationDataBrokerState broker;
};

struct LogbookContext {
  stellar::sim::Universe& universe;
  const stellar::sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};

  std::vector<stellar::sim::LogbookEntry>& logbook;

  // Bank is the pre-multiplier exploration data value (kept for backward compatibility).
  double explorationDataBankCr{0.0};

  // UI callbacks.
  std::function<bool(stellar::sim::SystemId)> plotRouteToSystem;
  std::function<void(stellar::sim::StationId)> targetStation;
  std::function<void(std::string_view, double)> toast;
};

// Draw the persistent logbook window.
void drawLogbookWindow(LogbookWindowState& st, const LogbookContext& ctx);

// Draw a station-side exploration data brokerage panel (supports partial sales).
//
// `explorationDataBankCr` is reduced by the BASE value sold (so legacy HUD stays consistent).
// `creditsCr` is increased by the PAYOUT after premiums.
void drawExplorationDataBrokerPanel(ExplorationDataBrokerState& st,
                                    stellar::sim::Universe& universe,
                                    const stellar::sim::StarSystem& saleSystem,
                                    const stellar::sim::Station& saleStation,
                                    double timeDays,
                                    bool canSellHere,
                                    std::vector<stellar::sim::LogbookEntry>& logbook,
                                    double& explorationDataBankCr,
                                    double& creditsCr,
                                    std::function<void(std::string_view, double)> toast,
                                    std::function<void(stellar::core::u32, double)> addRep);

} // namespace stellar::game
