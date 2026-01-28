#pragma once

#include "stellar/econ/Commodity.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h" // SystemId / StationId

#include <array>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace stellar::sim {
class Universe;
struct StarSystem;
struct Mission;
struct FactionReputation;
struct SystemSecurityDeltaState;
} // namespace stellar::sim

namespace stellar::game {

// Small, pre-computed hint from the progression system.
// (We keep it lightweight and pass it in from main.cpp to avoid tightly coupling
// this window to progression internals.)
struct CopilotProgressionHint {
  std::string title;
  std::string detail;
  float progress01{0.0f};
};

// Runtime context for the Copilot window.
//
// This is intentionally "thin": it is mostly read-only data, plus optional callbacks
// that the window can invoke when the player clicks action buttons.
struct CopilotContext {
  sim::Universe& universe;
  const sim::StarSystem* currentSystem{nullptr};

  double timeDays{0.0};
  double timeRealSec{0.0};

  // Estimated max jump range (LY) used by planning helpers.
  double maxJumpLy{0.0};

  // Optional inputs for risk-aware planning.
  std::span<const sim::FactionReputation> playerRepWithFaction{};
  std::span<const sim::SystemSecurityDeltaState> securityDeltas{};

  // Player position (km) for in-system distance checks.
  math::Vec3d shipPosKm{0.0, 0.0, 0.0};

  // Ship status.
  double hull{0.0};
  double hullMax{0.0};
  double shield{0.0};
  double shieldMax{0.0};
  double fuel{0.0};
  double fuelMax{0.0};
  double heat{0.0};
  int fuelScoopMk{0};

  // Docking / navigation state.
  bool docked{false};
  sim::StationId dockedStationId{0};

  bool navAutoRun{false};
  bool supercruiseActive{false};
  bool fsdBusy{false};
  bool dockingAutopilot{false};

  const std::vector<sim::SystemId>* navRoute{nullptr};
  std::size_t navRouteHop{0};

  // Cargo + missions.
  const std::array<double, econ::kCommodityCount>* cargo{nullptr};
  const std::vector<sim::Mission>* missions{nullptr};
  core::u64 trackedMissionId{0};

  CopilotProgressionHint progression{};

  // --- Optional callbacks ---
  std::function<void(core::u64 missionId)> trackMission;

  // Plot a route to a system (does not necessarily arm auto-run).
  std::function<void(sim::SystemId systemId)> plotRouteToSystem;

  // Plot a route to a system and set a preferred arrival station.
  std::function<void(sim::SystemId systemId, sim::StationId stationId)> routeToStation;

  // High-level travel action (typically via Integration Hub so it is traceable).
  std::function<void(sim::SystemId systemId, sim::StationId stationId, bool armAutoRun)> goToStation;

  // In-system targeting / docking.
  std::function<void(sim::StationId stationId)> targetStation;
  std::function<void(sim::StationId stationId)> requestDocking;
  std::function<void()> engageDockingComputer;

  // UI helpers.
  std::function<void()> openStationServices;
  std::function<void()> openMissionsWindow;
  std::function<void()> openMissionControl;
  std::function<void()> openTradePlanner;
  std::function<void()> openProgressionWindow;

  std::function<void(std::string_view msg, double ttlSec)> toast;
};

struct CopilotItineraryStopSummary {
  sim::SystemId systemId{0};
  sim::StationId stationId{0};
  bool isSite{false};
  bool reachable{true};

  int missionCount{0};
  double totalRewardCr{0.0};
  double avgRisk01{0.0};

  double etaDays{0.0};
  double earliestDeadlineDay{0.0};
  double etaSlackHours{0.0};

  std::vector<core::u64> missionIds{};
};

struct CopilotItineraryCache {
  core::u64 key{0};
  double builtAtRealSec{-1.0};
  bool ok{false};
  std::string status{};

  double totalRewardCr{0.0};
  int unreachableStops{0};

  std::vector<CopilotItineraryStopSummary> stops{};
};

struct CopilotWindowState {
  bool open{false};

  bool showRecommendations{true};
  bool showShipStatus{true};
  bool showMissions{true};
  bool showPlaybook{true};
  bool showProgression{true};

  bool compact{false};

  // --- Playbook / itinerary planner ---
  bool playbookGroupBySystem{true};
  int playbookMaxStops{8};
  float playbookRecomputeSec{2.0f};

  CopilotItineraryCache playbook{};


  float dismissTtlSec{120.0f};
  bool showDismissed{false};

  char filter[128]{};

  // suggestionId -> epoch seconds until which it is hidden.
  std::unordered_map<std::string, double> dismissedUntilRealSec;
};

void drawCopilotWindow(CopilotWindowState& st, const CopilotContext& ctx);

} // namespace stellar::game
