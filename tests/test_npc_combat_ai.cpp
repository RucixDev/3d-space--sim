#include "stellar/sim/NpcCombatAI.h"

#include "test_harness.h"

#include <iostream>

using namespace stellar;

static void checkPipsInvariant(const sim::Pips& p, int& failures) {
  CHECK(p.eng >= 0);
  CHECK(p.wep >= 0);
  CHECK(p.sys >= 0);
  CHECK(p.eng <= sim::kPipMax);
  CHECK(p.wep <= sim::kPipMax);
  CHECK(p.sys <= sim::kPipMax);
  CHECK((p.eng + p.wep + p.sys) == sim::kPipTotal);
}

int test_npc_combat_ai() {
  int failures = 0;

  const core::u64 seed = 0xC0FFEEu;

  {
    sim::NpcPipContext ctx{};
    ctx.role = sim::NpcCombatRole::Pirate;
    ctx.aiSkill = 0.4;
    auto d = sim::decideNpcPips(ctx, seed);
    checkPipsInvariant(d.pips, failures);
    CHECK(d.orbitStrafe01 == 0.0);
  }

  {
    // Traders who are fleeing should strongly prefer ENG and allow boosting.
    sim::NpcPipContext ctx{};
    ctx.role = sim::NpcCombatRole::Trader;
    ctx.aiSkill = 0.2;
    ctx.fleeing = true;
    ctx.inCombat = true;
    auto d = sim::decideNpcPips(ctx, seed + 1);
    checkPipsInvariant(d.pips, failures);
    CHECK(d.pips.eng >= 3);
    CHECK(d.allowBoost);
    CHECK(d.orbitStrafe01 > 0.0);
  }

  {
    // In combat with low shields, SYS should get attention.
    sim::NpcPipContext ctx{};
    ctx.role = sim::NpcCombatRole::Police;
    ctx.aiSkill = 0.6;
    ctx.inCombat = true;
    ctx.underFire = true;
    ctx.shieldFrac = 0.18;
    auto d = sim::decideNpcPips(ctx, seed + 2);
    checkPipsInvariant(d.pips, failures);
    CHECK(d.pips.sys >= 2);
    CHECK(d.orbitStrafe01 > 0.0);
  }

  {
    // If we want to fire and WEP cap is low, bias toward WEP.
    sim::NpcPipContext ctx{};
    ctx.role = sim::NpcCombatRole::Pirate;
    ctx.aiSkill = 0.8;
    ctx.inCombat = true;
    ctx.wantsToFire = true;
    ctx.wepCapFrac = 0.08;
    auto d = sim::decideNpcPips(ctx, seed + 3);
    checkPipsInvariant(d.pips, failures);
    CHECK(d.pips.wep >= 3);
    CHECK(d.allowBoost);
  }

  {
    // Higher skill should (usually) yield at least as much orbit/jink strafe.
    sim::NpcPipContext lo{};
    lo.role = sim::NpcCombatRole::Pirate;
    lo.aiSkill = 0.1;
    lo.inCombat = true;

    sim::NpcPipContext hi = lo;
    hi.aiSkill = 0.9;

    auto dlo = sim::decideNpcPips(lo, seed + 4);
    auto dhi = sim::decideNpcPips(hi, seed + 4);

    CHECK(dhi.orbitStrafe01 >= dlo.orbitStrafe01);
  }

  return failures;
}
