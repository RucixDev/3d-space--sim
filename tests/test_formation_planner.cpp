#include "test_harness.h"

#include "stellar/sim/FormationPlanner.h"

#include <cmath>
#include <vector>

int test_formation_planner() {
  int failures = 0;

  using namespace stellar;

  const math::Vec3d leaderPos{0, 0, 0};
  const math::Vec3d leaderVel{0, 0, 20};
  const math::Vec3d leaderFwd{0, 0, 1};
  const math::Vec3d leaderUp{0, 1, 0};

  sim::FormationPlanParams pp{};
  pp.formation.radiusKm = 100.0;
  pp.formation.lateralScale = 1.0;
  pp.formation.verticalScale = 0.0;
  pp.formation.behindScale = 1.0;
  pp.formation.rowSpacingScale = 0.60;
  pp.formation.jitterKm = 0.0;

  const core::u64 seed = 0xBADC0FFEE0DDF00Dull;

  // --- Case 1: swapped left/right positions should map to nearest slots (avoid crossovers). ---
  {
    std::vector<sim::FormationMember> members;
    members.push_back(sim::FormationMember{10ull, math::Vec3d{+60, 0, -90}}); // visually "right"
    members.push_back(sim::FormationMember{20ull, math::Vec3d{-60, 0, -90}}); // visually "left"

    const sim::FormationPlan plan = sim::planFormation(
      leaderPos,
      leaderVel,
      leaderFwd,
      leaderUp,
      std::span<const sim::FormationMember>(members.data(), members.size()),
      seed,
      pp
    );

    CHECK(plan.pattern == sim::FormationPattern::Wedge);

    auto slotFor = [&](const sim::FormationPlan& p, core::u64 id) {
      for (const auto& a : p.assignments) {
        if (a.memberId == id) return a.slotIndex;
      }
      return -1;
    };

    const int slot10 = slotFor(plan, 10ull);
    const int slot20 = slotFor(plan, 20ull);

    // In wedge formation, slot 0 is left (negative X), slot 1 is right (positive X).
    CHECK(slot10 == 1);
    CHECK(slot20 == 0);

    // Determinism.
    const sim::FormationPlan plan2 = sim::planFormation(
      leaderPos,
      leaderVel,
      leaderFwd,
      leaderUp,
      std::span<const sim::FormationMember>(members.data(), members.size()),
      seed,
      pp
    );
    CHECK(plan2.pattern == plan.pattern);
    CHECK(slotFor(plan2, 10ull) == 1);
    CHECK(slotFor(plan2, 20ull) == 0);
  }

  // --- Case 2: stickiness should influence symmetric tie-cases. ---
  {
    std::vector<sim::FormationMember> members;
    members.push_back(sim::FormationMember{10ull, math::Vec3d{0, 0, -100}});
    members.push_back(sim::FormationMember{20ull, math::Vec3d{0, 0, -100}});

    // Without prior hints, we expect the deterministic solver ordering:
    // smaller id -> lower slot.
    sim::FormationPlanParams p0 = pp;
    p0.costQuantizeMeters = 1.0;
    p0.stickySlotPenaltyMeters = 0.0;

    const sim::FormationPlan base = sim::planFormation(
      leaderPos,
      leaderVel,
      leaderFwd,
      leaderUp,
      std::span<const sim::FormationMember>(members.data(), members.size()),
      seed,
      p0
    );

    auto slotFor = [&](const sim::FormationPlan& plan, core::u64 id) {
      for (const auto& a : plan.assignments) {
        if (a.memberId == id) return a.slotIndex;
      }
      return -1;
    };

    CHECK(slotFor(base, 10ull) == 0);
    CHECK(slotFor(base, 20ull) == 1);

    // With prior hints + penalty, the solver should prefer to keep the hinted slots.
    sim::FormationPlanParams p1 = pp;
    p1.costQuantizeMeters = 1.0;
    p1.stickySlotPenaltyMeters = 1000.0; // 1 km worth of penalty is enough for tie-breaking.

    const sim::FormationSlotHint hints[] = {
      {10ull, 1},
      {20ull, 0},
    };

    const sim::FormationPlan sticky = sim::planFormation(
      leaderPos,
      leaderVel,
      leaderFwd,
      leaderUp,
      std::span<const sim::FormationMember>(members.data(), members.size()),
      seed,
      p1,
      std::span<const sim::FormationSlotHint>(hints, 2)
    );

    CHECK(slotFor(sticky, 10ull) == 1);
    CHECK(slotFor(sticky, 20ull) == 0);
  }

  return failures;
}
