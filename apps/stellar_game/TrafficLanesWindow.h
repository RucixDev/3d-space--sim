#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/System.h"
#include "stellar/sim/TrafficLanes.h"
#include "stellar/sim/TrafficLedger.h"

#include <functional>
#include <string_view>

namespace stellar::sim {
class Universe; // fwd
}

namespace stellar::game {

// -----------------------------------------------------------------------------
// TrafficLanesWindow
// -----------------------------------------------------------------------------
//
// A developer-facing / gameplay helper UI for visualizing the deterministic
// traffic convoy layer (trade lanes) in the current system.
//
// This window is intentionally lightweight: it does not own convoy simulation.
// Instead, it regenerates a "view" list on demand from the TrafficLedger (if
// available) and TrafficLaneParams. Targeting / supercruise actions are
// delegated to callbacks provided by the game runtime.
//

struct TrafficLanesWindowState {
  bool open{false};

  // List/map UI
  stellar::core::u64 selectedConvoyId{0};
  bool showInactive{false};
  bool showMap{true};

  // Map display
  int pathSegments{48}; // polyline resolution for lane arcs
  float mapHeightPx{280.0f};
};

struct TrafficLanesContext {
  stellar::sim::Universe& universe;
  const stellar::sim::StarSystem* currentSystem{nullptr};
  const stellar::sim::TrafficLedger* trafficLedger{nullptr};

  // Time + player kinematics for relative readouts.
  double timeDays{0.0};
  stellar::math::Vec3d playerPosKm{0, 0, 0};
  stellar::math::Vec3d playerVelKmS{0, 0, 0};

  // Lane generation parameters (usually the game's live tuning values).
  stellar::sim::TrafficLaneParams laneParams{};

  // Optional filters/callbacks.
  std::function<bool(stellar::core::u64 convoyId)> isConvoySuppressed; // e.g., interdicted/destroyed
  std::function<void(stellar::core::u64 convoyId)> targetTrafficConvoy; // set current target to this convoy
  std::function<void()> requestSupercruise; // (optional) start supercruise toward current target

  // UI feedback.
  std::function<void(std::string_view msg, double ttlSec)> toast;
};

// Draw the Traffic Lanes window. (No-op if st.open == false.)
void drawTrafficLanesWindow(TrafficLanesWindowState& st, const TrafficLanesContext& ctx);

} // namespace stellar::game
