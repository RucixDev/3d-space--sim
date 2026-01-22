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



// Optional: actionable "live" security demands (authority scans / bounty submission).
// This bridges Comms messages back into the moment-to-moment gameplay without turning the
// sim::Comms layer into a quest scripting system.
struct SecurityDemandUi {
  enum class Kind : core::u8 { None = 0, BribeOffer, BountyDemand };
  Kind kind{Kind::None};

  core::u32 factionId{0};
  std::string authorityName;
  std::string detail;

  // Kind::BribeOffer
  double bribeCr{0.0};
  double fineCr{0.0};

  // Kind::BountyDemand
  double bountyCr{0.0};

  // Timing (seconds)
  double secondsLeft{0.0};
  double secondsTotal{0.0};

  // Capability flags for UI disable/tooltip.
  bool canPay{false};
  bool actionAllowed{true};
};

struct CommsWindowContext {
  sim::Universe* universe{nullptr};
  const sim::StarSystem* currentSystem{nullptr};
  double timeDays{0.0};
  double timeRealSec{0.0};

  sim::CommsLog* log{nullptr};

  std::function<void(const std::string& msg, double ttlSec)> toast;
  std::function<void(sim::SystemId sysId, sim::StationId stId)> plotTo;

  // Optional high-level travel hook.
  // When provided, the Comms UI can offer one-click travel actions that
  // optionally arm auto-run/docking.
  std::function<void(sim::SystemId sysId, sim::StationId stId, bool armAutoRun)> goTo;

  // Mission deep-links. Mission briefings carry the mission id in sim::CommsMessage::relatedId.
  // These hooks let the Comms UI bridge into the mission tracker and nav/autopilot.
  std::function<void(core::u64 missionId)> trackMission;
  std::function<void(core::u64 missionId, bool armAutoRun)> syncNavToMission;

  // Optional: make Security transmissions actionable (pay bribe / comply / submit bounty).
  // When provided, the Comms UI can surface live response buttons on the selected message.
  std::function<SecurityDemandUi(const sim::CommsMessage& selected)> querySecurityDemand;
  std::function<void()> actSecurityBribe;
  std::function<void()> actSecurityComplyOrSubmit;

  // Hotkey display strings (e.g. "C", "I"). Purely cosmetic.
  std::string securityBribeChord;
  std::string securityComplyChord;
};

void enqueueCommsOverlay(CommsOverlayState& ov, core::u64 messageId);
void tickAndDrawCommsOverlay(CommsOverlayState& ov, const sim::CommsLog& log, double timeRealSec);

void drawCommsWindow(CommsWindowState& st, CommsOverlayState& overlay, CommsWindowContext& ctx);

} // namespace stellar::game
