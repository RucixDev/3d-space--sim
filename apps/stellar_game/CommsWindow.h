#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Comms.h"

#include <functional>

namespace stellar::sim {
class Universe;
struct StarSystem;
} // namespace stellar::sim

namespace stellar::game {

struct CommsWindowState {
  bool open{false};

  // UX
  bool autoMarkReadOnSelect{true};
  bool unreadOnly{false};
  bool newestFirst{true};
  bool pinnedFirst{true};
  bool wrapBody{true};

  // Filters
  int channelFilter{-1}; // -1 = all, else cast to sim::CommsChannel
  char filter[128]{};

  // Selection
  core::u64 selectedId{0};
  double selectedOpenedSec{0.0};
};

// Non-interactive, diegetic "incoming transmission" overlay.
struct CommsOverlayState {
  bool enabled{true};

  // Timing
  double baseHoldSec{5.5};
  double fadeOutSec{0.6};

  // Layout
  float widthPx{520.0f};
  float padPx{10.0f};
  float marginPx{18.0f};

  // Queue
  std::vector<core::u64> queue;
  core::u64 activeId{0};
  double activeStartSec{0.0};
  double activeUntilSec{0.0};
};

struct CommsWindowContext {
  sim::Universe* universe{nullptr};
  const sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};
  double timeRealSec{0.0};

  sim::CommsLog* log{nullptr};

  std::function<void(const std::string& msg, double ttlSec)> toast;
  std::function<void(sim::SystemId sysId, sim::StationId stId)> plotTo;
};

void enqueueCommsOverlay(CommsOverlayState& ov, core::u64 messageId);
void tickAndDrawCommsOverlay(CommsOverlayState& ov, const sim::CommsLog& log, double timeRealSec);

void drawCommsWindow(CommsWindowState& st, CommsOverlayState& overlay, CommsWindowContext& ctx);

} // namespace stellar::game
