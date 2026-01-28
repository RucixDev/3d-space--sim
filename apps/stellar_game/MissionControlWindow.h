#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/MissionPlanner.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Universe.h"

#include <functional>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace stellar::game {

// Mission Control (UI)
//
// A high-level mission planning dashboard:
//  - builds a greedy itinerary across active missions (via sim::MissionPlanner)
//  - presents per-stop risk/reward/deadline analytics
//  - optional "runner" that can auto-plot the next stop as you complete legs

struct MissionControlWindowState {
  bool open{false};

  // --- Itinerary planner settings ---
  bool groupBySystem{false};
  int routeMode{0}; // 0=hops, 1=distance

  // When > 0, prefer safer inter-system routes by adding a risk penalty per LY.
  // (See sim::MissionItineraryParams::navRiskWeightPerLy)
  double navRiskWeightPerLy{0.0};

  double maxJumpLyOverride{0.0}; // 0 => use ctx.maxJumpLy
  double queryRadiusLy{0.0};     // 0 => auto
  int maxSystems{1800};
  int maxStops{12};

  double rewardWeight{1.0};
  double riskWeight{0.45};
  double urgencyWeight{0.35};

  // ETA model (used for deadline-aware urgency + UI analytics).
  bool etaAwareUrgency{true};
  double etaSecondsPerJump{45.0};
  double etaSecondsPerLy{0.0};
  double etaSecondsPerStop{600.0};
  double etaSecondsPerSite{420.0};

  bool autoRebuild{true};

  // Mission selection.
  bool selectionInitialized{false};
  bool includeFailed{false};
  bool includeCompleted{false};
  std::unordered_set<stellar::core::u64> selectedMissionIds;

  stellar::core::u64 focusedMissionId{0};

  // Cached plan.
  stellar::core::u64 cacheKey{0};
  stellar::sim::MissionItineraryResult plan;

  // --- Runner ---
  struct Runner {
    bool active{false};
    bool armAutoRun{true};
    bool autoPlotNextStop{true};
    bool autoTrackOnAdvance{true};

    int stopIndex{0};
    double startedAtDays{0.0};
    double lastAdvanceAtDays{0.0};

    bool startToastEmitted{false};
    bool replotPending{false};

    std::vector<stellar::sim::MissionItineraryStop> stops;
  } runner;
};

struct MissionControlContext {
  stellar::sim::Universe& universe;
  const stellar::sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};
  double timeRealSec{0.0};

  // Player state.
  bool docked{false};
  stellar::sim::StationId dockedStationId{0};
  double maxJumpLy{0.0};

  // Missions.
  const std::vector<stellar::sim::Mission>* missions{nullptr};
  stellar::core::u64 trackedMissionId{0};
  std::span<const stellar::sim::FactionReputation> playerRepWithFaction{};
  std::span<const stellar::sim::SystemSecurityDeltaState> securityDeltas{};

  // UI callbacks.
  std::function<void(stellar::core::u64)> trackMission;

  std::function<void(stellar::sim::SystemId)> plotRouteToSystem;
  std::function<void(stellar::sim::SystemId, stellar::sim::StationId)> routeToStation;
  std::function<void(stellar::sim::SystemId, stellar::sim::StationId, bool armAutoRun)> goToStation;

  std::function<void(std::string_view, double)> toast;
};

// Runner tick (call even when the window is closed).
void tickMissionControl(MissionControlWindowState& st, const MissionControlContext& ctx);

void drawMissionControlWindow(MissionControlWindowState& st, const MissionControlContext& ctx);

} // namespace stellar::game
