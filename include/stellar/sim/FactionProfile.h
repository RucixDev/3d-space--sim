#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Faction.h"

namespace stellar::sim {

// Procedural traits for factions.
//
// These are lightweight, deterministic "flavor" parameters that gameplay code can
// use to derive other systems (law strictness, corruption/bribe frequency, trade
// bias, mission rewards, diplomacy, etc.).
//
// All fields are normalized to [0, 1] and computed deterministically from
// (universeSeed, factionId).
struct FactionProfile {
  core::u32 factionId{0};

  // 0 = anarchic / lax, 1 = strict / authoritarian.
  double authority{0.5};

  // 0 = honest, 1 = corrupt (bribes common).
  double corruption{0.5};

  // 0 = poor, 1 = wealthy.
  double wealth{0.5};

  // 0 = unstable, 1 = stable.
  double stability{0.5};

  // 0 = low tech, 1 = high tech.
  double tech{0.5};

  // 0 = pacifist, 1 = militaristic.
  double militarism{0.5};
};

// Deterministic per-faction profile.
FactionProfile factionProfile(core::u64 universeSeed, core::u32 factionId);

// A deterministic, symmetric relation score between two factions in [-1, 1].
// Positive => more allied, Negative => more hostile.
//
// This intentionally mixes "ideological" similarity (FactionProfile) and a rough
// "border pressure" heuristic (home distance relative to influence radii).
double factionRelation(core::u64 universeSeed, const Faction& a, const Faction& b);

enum class FactionRelationKind : core::u8 {
  Hostile = 0,
  Neutral = 1,
  Allied = 2,
};

FactionRelationKind classifyFactionRelation(double relation);
const char* factionRelationKindName(FactionRelationKind k);

} // namespace stellar::sim
