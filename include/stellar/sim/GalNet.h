#pragma once

#include "stellar/core/Types.h"

#include "stellar/sim/Comms.h"
#include "stellar/sim/Faction.h"
#include "stellar/sim/SystemConditions.h"

#include <string>
#include <string_view>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// GalNet bulletins (headless)
// -----------------------------------------------------------------------------
//
// The SDL/OpenGL prototype historically built system-event bulletins inline in
// apps/stellar_game/main.cpp.
//
// This module extracts the deterministic bulletin composition + auto-broadcast
// decision logic into the core library so it can be:
//   - used by UI layers without duplicating strings / logic
//   - tested headlessly
//   - reused by tools (e.g., batch scanning watched systems)

// Persistent per-system state used to avoid repeating bulletins within a cycle.
//
// This mirrors the runtime bookkeeping that was historically kept in the game
// app (GalNetAnnounceState).
struct GalNetAnnounceState {
  double lastCycleStartDay{-1.0};
  bool lastCycleHadActiveEvent{false};
};

struct GalNetAutoBroadcastDecision {
  bool shouldPublish{false};
  // If true, the caller should publish even if the current cycle has no active
  // event (typically used for "event ended" notifications).
  bool allowWhenNoEvent{false};
};

// Decide whether a system should auto-publish a bulletin for the current cycle.
//
// This function is intended to be called once per cycle boundary (or whenever
// the caller detects that the system's event startDay changed).
//
// - ioState is updated to the current cycle regardless of whether publishing is
//   enabled, preventing retroactive "spam" when the player toggles settings.
// - autoEnabled controls whether publishing is allowed (local or watched).
// - minSeverity01 filters low-severity active events.
// - broadcastEventEnds allows a bulletin to be emitted when an event transitions
//   from active -> inactive.
GalNetAutoBroadcastDecision galNetMaybeAutoBroadcast(GalNetAnnounceState& ioState,
                                                     const SystemEvent& ev,
                                                     double minSeverity01,
                                                     bool autoEnabled,
                                                     bool broadcastEventEnds);

// Friendly labels and short tips for system event kinds.
// Returned strings are stable string literals.
const char* galNetSystemEventKindLabel(SystemEventKind kind);
const char* galNetSystemEventKindTip(SystemEventKind kind);

// Resolve a faction name from a faction list.
//
// - factionId==0 returns "Independent".
// - Unknown ids return "Faction <id>".
std::string factionNameOrId(core::u32 factionId, const std::vector<Faction>& factions);

struct GalNetBulletin {
  CommsMessage msg{};
  // Rough "importance" score for UI filtering. Conventionally:
  //   - active event: severity01
  //   - no event: 0.25
  double importance01{0.25};
  bool hasActiveEvent{false};
};

struct GalNetBulletinResult {
  bool ok{false};
  const char* reason{nullptr};
  GalNetBulletin bulletin{};
};

// Compose a GalNet bulletin from a precomputed system conditions snapshot.
//
// controllingFactionName is optional; if empty the bulletin omits the name and
// only includes the controlling faction id (if any).
GalNetBulletinResult makeGalNetBulletin(const StarSystem& sys,
                                        const SystemConditionsSnapshot& snap,
                                        double timeDays,
                                        std::string_view controllingFactionName = {},
                                        std::string_view contextTag = {},
                                        bool allowWhenNoEvent = false);

class Universe;

// Convenience: compute the snapshot from Universe state and build the bulletin.
//
// This respects Universe hooks for security deltas and event params.
GalNetBulletinResult makeGalNetBulletin(const Universe& universe,
                                        SystemId systemId,
                                        double timeDays,
                                        std::string_view contextTag = {},
                                        bool allowWhenNoEvent = false);

} // namespace stellar::sim
