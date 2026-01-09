#pragma once

#include "stellar/core/Types.h"

namespace stellar::sim {

// Lightweight squad tactics helpers used by the prototype game.
//
// Motivation:
//  - The game spawns pirates/police in small packs (groupId + leaderId).
//  - Historically each ship acted independently.
//  - These helpers add *group-level* morale and reinforcement logic in a deterministic,
//    unit-testable way.

struct SquadSnapshot {
  int totalCount{0};
  int aliveCount{0};
  int underFireCount{0};

  bool leaderAlive{true};

  // Averages over alive members (0..1).
  double avgHullFrac{1.0};
  double avgShieldFrac{1.0};
};

struct PirateSquadContext {
  // System knobs (0..1). Higher security discourages pirates.
  double security01{0.5};
  double piracy01{0.5};

  // Estimate of how "worth it" the target is.
  double playerCargoValueCr{0.0};

  // True if this squad is committed to a raid (e.g., convoy ambush).
  // Committed squads break less easily.
  bool committedRaid{false};
};

struct PirateSquadDecision {
  double morale01{0.0};
  bool retreat{false};
  double retreatDurationSec{0.0};
};

PirateSquadDecision evaluatePirateSquad(const SquadSnapshot& snap,
                                       const PirateSquadContext& ctx,
                                       core::u64 universeSeed,
                                       core::u64 groupId);

struct PoliceSquadContext {
  double security01{0.5};

  // Player legal status in the local system.
  bool playerWanted{false};
  double playerBountyCr{0.0};

  // Current pursuit intensity (prototype: ~[0..6]).
  double policeHeat{0.0};

  // True if this squad is actively fighting the player (not just engaging pirates).
  bool fightingPlayer{false};
};

struct PoliceSquadDecision {
  double morale01{0.0};

  // Police rarely retreat, but can in extreme cases.
  bool retreat{false};
  double retreatDurationSec{0.0};

  // Request additional patrols/interceptors.
  bool callReinforcements{false};
  double reinforcementUrgency01{0.0};
};

PoliceSquadDecision evaluatePoliceSquad(const SquadSnapshot& snap,
                                       const PoliceSquadContext& ctx,
                                       core::u64 universeSeed,
                                       core::u64 groupId);

} // namespace stellar::sim
