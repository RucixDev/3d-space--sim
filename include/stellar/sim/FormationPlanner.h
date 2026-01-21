#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Formation.h"

#include <span>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// FormationPlanner
// -----------------------------------------------------------------------------
//
// Formation.h can generate formation slot offsets, but callers often assign
// slots by a stable id ordering. In practice, that can cause wingmen to
// "cross" through each other (e.g., when two escorts already sit on opposite
// sides of the leader).
//
// FormationPlanner fixes this by:
//  1) Generating the candidate slot positions for the current leader pose.
//  2) Solving a min-cost assignment (Hungarian / Kuhn–Munkres) to match
//     members to slots, minimizing travel distance.
//
// The planner is headless and purely geometric: it does not know about ships
// beyond (id, position). Callers decide how to steer toward the returned slots.
//
// Determinism:
//  - Costs are computed using quantized integer distances.
//  - Members are sorted by id before solving.
//  - Ties are broken by deterministic iteration order in the solver.
//
// Performance:
//  - The Hungarian solver is O(n^3). This is intended for small groups
//    (escort wings / squads), not large swarms.

struct FormationMember {
  core::u64 id{0};
  math::Vec3d posKm{0, 0, 0};
};

struct FormationSlotHint {
  core::u64 memberId{0};
  // -1 means no hint. Otherwise, the slot index the member previously occupied.
  int priorSlot{-1};
};

struct FormationPlanParams {
  // Formation geometry.
  FormationParams formation{};

  // Quantization for the assignment cost.
  // Cost is computed as round(distanceMeters / costQuantizeMeters).
  // Smaller values give more fidelity but larger integers.
  double costQuantizeMeters{1.0};

  // Optional stickiness to avoid slot thrash when two assignments are near-equal.
  // When a prior slot is provided for a member, switching away from it adds:
  //   stickySlotPenaltyMeters / costQuantizeMeters
  // to the cost.
  double stickySlotPenaltyMeters{0.0};
};

struct FormationAssignment {
  core::u64 memberId{0};
  int slotIndex{0};

  // Slot target position in world space (km) for the given leader pose.
  math::Vec3d slotPosKm{0, 0, 0};
};

struct FormationPlan {
  FormationPattern pattern{FormationPattern::Trail};
  Basis basis{};

  // One entry per member, sorted by member id.
  std::vector<FormationAssignment> assignments;
};

// Plan an escort-style formation around a leader pose.
//
// Inputs:
//  - leaderPosKm: leader world position (km)
//  - leaderVelKmS: leader world velocity (km/s) (used as a forward hint)
//  - leaderForwardWorld: fallback forward hint when velocity is near zero
//  - leaderUpWorld: up hint (for basis generation)
//  - members: group members to be placed into formation slots
//  - seed: deterministic seed used for slot jitter
//  - params: tuning
//  - priorSlots: optional prior slot hints (for stickiness)
//
// Returns:
//  - A FormationPlan with stable slot assignments.
FormationPlan planFormation(const math::Vec3d& leaderPosKm,
                            const math::Vec3d& leaderVelKmS,
                            const math::Vec3d& leaderForwardWorld,
                            const math::Vec3d& leaderUpWorld,
                            std::span<const FormationMember> members,
                            core::u64 seed,
                            const FormationPlanParams& params = {},
                            std::span<const FormationSlotHint> priorSlots = {});

} // namespace stellar::sim
