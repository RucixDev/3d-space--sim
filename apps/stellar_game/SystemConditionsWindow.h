#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/System.h"
#include "stellar/sim/SystemEvents.h"
#include "stellar/sim/SystemSecurityDynamics.h"

#include <functional>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace stellar::sim {
class Universe;
struct SystemStub;
} // namespace stellar::sim

namespace stellar::game {

struct SystemConditionsWindowState {
  bool open{false};

  // GalNet bulletin settings.
  // These are gameplay-facing toggles that control when the game posts system-event
  // bulletins to the player's Comms inbox.
  bool galNetAutoBroadcastLocal{true};
  bool galNetAutoBroadcastWatched{true};
  bool galNetBroadcastEventEnds{false};

  // Presentation toggles for auto-broadcasts.
  bool galNetLocalShowOverlay{true};
  bool galNetLocalShowToast{true};
  bool galNetWatchedShowOverlay{false};
  bool galNetWatchedShowToast{false};

  // Ignore events below this severity when auto-broadcasting.
  float galNetMinSeverity{0.0f};

  // Query controls.
  float radiusLy{220.0f};
  int maxResults{96};

  bool onlyActiveEvents{false};
  bool includeCurrentSystem{true};

  enum class SortKey : int {
    Distance = 0,
    Security = 1,
    Piracy = 2,
    Traffic = 3,
    EventSeverity = 4,
  };

  SortKey sortKey{SortKey::Distance};
  bool sortDescending{false};

  // Cache (updated when the day stamp or settings change).
  stellar::sim::SystemId cacheFromSystem{0};
  int cacheDayStamp{-999999};
  float cacheRadiusLy{0.0f};
  int cacheMaxResults{0};
  bool cacheOnlyActiveEvents{false};
  bool cacheIncludeCurrentSystem{true};
  SortKey cacheSortKey{SortKey::Distance};
  bool cacheSortDescending{false};

  double lastComputeMs{0.0};

  struct Row {
    stellar::sim::SystemId id{0};
    double distanceLy{0.0};
    stellar::core::u32 controllingFactionId{0};

    stellar::sim::SystemEvent event{};
    stellar::sim::SystemSecurityProfile base{};
    stellar::sim::SystemSecurityProfile effective{};
  };

  std::vector<Row> cacheRows;
};

struct SystemConditionsWindowContext {
  stellar::sim::Universe& universe;
  const stellar::sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};

  // Optional: player-driven security delta states (from SaveGame).
  const std::unordered_map<stellar::sim::SystemId, stellar::sim::SystemSecurityDeltaState>* systemSecurityDeltaBySystem{
      nullptr};

  // Tuning to match gameplay.
  stellar::sim::SystemSecurityDynamicsParams dynamicsParams{};
  stellar::sim::SystemEventParams eventParams{};

  // Optional UI callbacks.
  // Optional: GalNet integration.
  // The window can manage a player watchlist and publish system-event bulletins to Comms.
  std::vector<stellar::sim::SystemId>* watchSystems{nullptr};
  std::function<void(stellar::sim::SystemId systemId, bool showOverlay, bool showToast)> postGalNetBulletin;

  std::function<bool(stellar::sim::SystemId)> plotRouteToSystem;
  std::function<void(stellar::sim::SystemId)> targetSystem;
  std::function<void(std::string_view msg, double ttlSec)> toast;
};

// Draw the System Conditions window. (No-op if st.open == false.)
void drawSystemConditionsWindow(SystemConditionsWindowState& st, const SystemConditionsWindowContext& ctx);

} // namespace stellar::game
