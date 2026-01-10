#pragma once

#include "stellar/core/Types.h"

#include "stellar/sim/Celestial.h"
#include "stellar/sim/SecurityModel.h"

namespace stellar::sim {

// Lightweight, deterministic "system event" layer.
//
// The goal is to add some BGS-like texture (booms, busts, raids, crackdowns)
// without requiring a persistent simulation step. Events are generated
// deterministically from:
//   (universeSeed, systemId, timeDays)
// and the system's *current* security profile (typically baseline + player-driven
// SystemSecurityDeltaState impulses).
//
// Events contribute a temporary additive delta to the core security metrics:
//   security01, piracy01, traffic01
//
// The caller decides how to combine this with persistent SystemSecurityDynamics.

enum class SystemEventKind : core::u8 {
  None = 0,
  TradeBoom,
  TradeBust,
  PirateRaid,
  SecurityCrackdown,
  CivilUnrest,
  ResearchBreakthrough,
};

struct SystemEventParams {
  // Length of an "event cycle" in days. Events are stable within a cycle.
  // (Example: with cycleDays=6, the active event only changes every 6 in-game days.)
  double cycleDays{6.0};

  // Chance that a cycle has an event at all.
  double eventChance{0.80};

  // Severity range for generated events.
  double minSeverity{0.35};
  double maxSeverity{1.00};

  // Clamp each channel's absolute delta.
  double maxAbsDelta{0.22};
};

struct SystemEvent {
  bool active{false};
  SystemEventKind kind{SystemEventKind::None};
  SystemId systemId{0};

  // Inclusive start, exclusive end.
  double startDay{0.0};
  double endDay{0.0};

  // 0..1, higher means larger deltas.
  double severity01{0.0};

  // Additive deltas to apply to a SystemSecurityProfile.
  double securityDelta{0.0};
  double piracyDelta{0.0};
  double trafficDelta{0.0};
};

// Human-friendly label for UI/debug.
// Returns a stable string literal.
const char* systemEventKindName(SystemEventKind k);

// Generate the currently-active event for a system.
//
// `baseProfile` should typically be the system's baseline security profile with
// any persistent (player-driven) deltas already applied. The returned event is
// generated *from* baseProfile but is NOT applied to it.
SystemEvent generateSystemEvent(core::u64 universeSeed,
                                SystemId systemId,
                                double timeDays,
                                const SystemSecurityProfile& baseProfile,
                                const SystemEventParams& params = {});

// Apply an event's deltas to a security profile (clamps to [0,1]).
SystemSecurityProfile applySystemEventToProfile(const SystemSecurityProfile& base,
                                                const SystemEvent& ev);

} // namespace stellar::sim
