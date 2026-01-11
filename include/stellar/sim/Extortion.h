#pragma once

#include "stellar/core/Types.h"

namespace stellar::sim {

// Lightweight helper for planning ship-to-ship cargo extortion demands.
//
// Design goals:
//  - Deterministic from (universeSeed, targetId) so a given target produces a stable "demand" profile.
//  - Minimal inputs: caller supplies rough strength and local security.
//  - No gameplay state here; this is pure math and is unit-testable.
struct ExtortionDemandPlan {
  bool offer{false};           // If false, the caller should not offer an extortion option.
  double demandedValueCr{0.0}; // How much value the attacker demands (credits, approximate).
  double complyChance01{0.0};  // Probability the victim complies (0..1).

  // Suggested timings for game-layer UX.
  double suggestedCooldownSec{0.0}; // Cooldown before the same target will entertain another demand.
  double suggestedFleeSec{0.0};     // How long the victim will try to flee after complying/refusing.
};

// Plan an extortion demand.
//
// Inputs:
//  - targetCargoValueCr: rough total cargo value on the victim ship.
//  - attackerStrength/defenderStrength: caller-provided relative power measures.
//    These don't need to be in any particular units; only the ratio matters.
//  - systemSecurity01: effective local security (0..1).
//  - policeNearby: whether an authority response is likely/immediate.
ExtortionDemandPlan planExtortionDemand(core::u64 universeSeed,
                                        core::u64 targetId,
                                        double targetCargoValueCr,
                                        double attackerStrength,
                                        double defenderStrength,
                                        double systemSecurity01,
                                        bool policeNearby);

} // namespace stellar::sim
