#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Comms.h"
#include "stellar/sim/System.h"

#include <functional>
#include <string_view>
#include <vector>

namespace stellar::sim {
class Universe;
} // namespace stellar::sim

namespace stellar::game {

// A focused view of the player's Comms inbox that surfaces GalNet bulletins,
// provides quick actions (watch / plot / target) and helps reduce notification
// overload.
struct GalNetInboxWindowState {
  bool open{false};

  bool newestFirst{true};
  bool watchedOnly{false};
  bool activeOnly{false};
  bool wrapBody{true};
  bool autoMarkReadOnSelect{true};

  char filter[128]{};
  stellar::core::u64 selectedMsgId{0};

  float leftPaneWidth{380.0f};
};

struct GalNetInboxWindowContext {
  stellar::sim::Universe& universe;
  stellar::sim::CommsLog& comms;

  const stellar::sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};

  // Player watchlist of systems.
  std::vector<stellar::sim::SystemId>* watchSystems{nullptr};

  // Optional integrations.
  std::function<void(stellar::sim::SystemId systemId, bool showOverlay, bool showToast)> postGalNetBulletin;
  std::function<bool(stellar::sim::SystemId)> plotRouteToSystem;
  std::function<void(stellar::sim::SystemId)> targetSystem;
  std::function<void(std::string_view msg, double ttlSec)> toast;
};

// Draw the GalNet inbox window. (No-op if st.open == false.)
void drawGalNetInboxWindow(GalNetInboxWindowState& st, const GalNetInboxWindowContext& ctx);

} // namespace stellar::game
