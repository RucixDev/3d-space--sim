#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/PowerDistributor.h"

namespace stellar::sim {

// A tiny, deterministic helper for making NPC combat behavior feel less "flat".
//
// This module doesn't move ships directly; it provides a suggested pip allocation and
// a small set of behavioral knobs (boost allowance + lateral strafe amount) that the
// game layer can feed into its flight controllers.

enum class NpcCombatRole : int {
  Pirate = 0,
  Police = 1,
  Trader = 2,
};

struct NpcPipContext {
  NpcCombatRole role{NpcCombatRole::Pirate};

  // [0,1]
  double aiSkill{0.5};

  // Normalized health/energy inputs.
  double shieldFrac{1.0};
  double hullFrac{1.0};

  double engCapFrac{1.0};
  double wepCapFrac{1.0};
  double sysCapFrac{1.0};

  // Situation flags.
  bool inCombat{false};
  bool wantsToFire{false};
  bool wantsToBoost{false};
  bool fleeing{false};
  bool underFire{false};
  bool missileThreat{false};
};

struct NpcPipDecision {
  // Chosen [0..4] pips that sum to sim::kPipTotal.
  Pips pips{2, 2, 2};

  // Suggested additional lateral strafe magnitude in [0..1].
  // The caller decides the direction (orbiting / jinking / etc.).
  double orbitStrafe01{0.0};

  // If true, it's reasonable to allow boost in guidance for this frame.
  bool allowBoost{false};
};

// Deterministic pip + behavior decision.
//
// seed:
//   Any stable seed (e.g. hash(universeSeed, npcId)). Used only to break ties.
NpcPipDecision decideNpcPips(const NpcPipContext& ctx, core::u64 seed);

} // namespace stellar::sim
