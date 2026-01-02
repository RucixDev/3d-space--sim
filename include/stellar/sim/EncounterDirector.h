#pragma once

#include "stellar/core/Random.h"
#include "stellar/core/Types.h"

namespace stellar::sim {

// EncounterDirector
// -----------------
// A small, headless "spawn scheduler" for local-space encounters.
//
// The SDL prototype historically embedded this logic directly in
// apps/stellar_game/main.cpp. This module makes the scheduling rules reusable
// in tooling/tests and keeps the game TU smaller.
//
// Design goals:
//  - Deterministic per universe seed.
//  - Frame-rate independent scheduling (based on timeDays thresholds).
//  - Lightweight: this module only *decides when* to request a spawn.
//    The app decides *how* to instantiate contacts (loadouts, positions, etc.).

struct EncounterDirectorCounts {
  int alivePirates{0};
  int aliveTraders{0};
  int alivePolice{0};
  int aliveTotal{0};
};

struct EncounterDirectorContext {
  // Simulation time (days).
  double timeDays{0.0};
  // If false, no local-space encounters should be spawned.
  bool combatSimEnabled{true};

  // For trader spawns (requires at least one station to anchor behavior).
  int stationCount{0};

  // Local law / reputation context.
  core::u32 localFactionId{0};
  double localRep{0.0};       // roughly [-100, +100]
  double localBountyCr{0.0};  // credits
  double policeHeat{0.0};     // soft pursuit intensity (prototype tuning: ~[0..6])
  bool policeAlert{false};    // short-term "crime response" alert window
  bool playerWantedHere{false};

  // Optional system-level knobs (0..1). These let the app feed in a compact
  // notion of "how dangerous" or "how secure" a system feels.
  //
  // Defaults are 0.5 so existing tuning remains unchanged unless a caller
  // explicitly provides values.
  double security01{0.5}; // higher => more effective security response
  double piracy01{0.5};   // higher => pirates spawn more frequently
  double traffic01{0.5};  // higher => trader traffic spawns more frequently
};

struct EncounterDirectorState {
  // Next scheduled spawn times.
  double nextPirateSpawnDays{0.01};
  double nextTraderSpawnDays{0.008};
  double nextPoliceSpawnDays{0.006};

  // Dedicated RNG for scheduling. Keeps encounter cadence stable even if the
  // app consumes its own RNG for other gameplay events.
  core::SplitMix64 rng{};
};

// Returned by planPoliceSpawn.
struct PoliceSpawnPlan {
  int desiredPolice{0};   // how many police should exist in local space
  int spawnMaxCount{0};   // if > 0, request a pack spawn up to this size
};

// Create a new director state for a universe seed.
//
// This is intentionally *not* saved/loaded in SaveGame yet; it is used as an
// ephemeral spawn scheduler for the local-space sim.
EncounterDirectorState makeEncounterDirector(core::u64 universeSeed, double timeDays = 0.0);

// Spawn planners (call in order).
//
// These functions update the state's next*SpawnDays when a spawn is requested.
// The app should then instantiate contacts and update EncounterDirectorCounts
// before calling the next planner.
//
// Returns a max pack size to pass to the app's spawn routine (0 = no spawn).
int planPirateSpawn(EncounterDirectorState& state,
                    const EncounterDirectorContext& ctx,
                    const EncounterDirectorCounts& counts);

// Returns the number of traders to spawn this frame (0 or 1 currently).
int planTraderSpawn(EncounterDirectorState& state,
                    const EncounterDirectorContext& ctx,
                    const EncounterDirectorCounts& counts);

// Police plan includes a desired count (used by the app for AI/UI) and an
// optional pack spawn request.
PoliceSpawnPlan planPoliceSpawn(EncounterDirectorState& state,
                                const EncounterDirectorContext& ctx,
                                const EncounterDirectorCounts& counts);

// Helper exposed for UI/tests.
int computeDesiredPolice(const EncounterDirectorContext& ctx,
                         const EncounterDirectorCounts& counts);

} // namespace stellar::sim
