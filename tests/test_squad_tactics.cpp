#include "stellar/sim/SquadTactics.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

int test_squad_tactics() {
  int failures = 0;

  // Pirates: healthy 3-ship pack in a pirate-friendly system should be confident.
  {
    sim::SquadSnapshot s{};
    s.totalCount = 3;
    s.aliveCount = 3;
    s.underFireCount = 0;
    s.leaderAlive = true;
    s.avgHullFrac = 0.95;
    s.avgShieldFrac = 0.90;

    sim::PirateSquadContext ctx{};
    ctx.security01 = 0.20;
    ctx.piracy01 = 0.85;
    ctx.playerCargoValueCr = 4500.0;
    ctx.committedRaid = false;

    const auto a = sim::evaluatePirateSquad(s, ctx, /*universeSeed*/ 1234u, /*groupId*/ 55u);
    const auto b = sim::evaluatePirateSquad(s, ctx, /*universeSeed*/ 1234u, /*groupId*/ 55u);

    CHECK(std::abs(a.morale01 - b.morale01) < 1e-12);
    CHECK(a.morale01 > 0.55);
    CHECK(!a.retreat);
    CHECK(a.retreatDurationSec == 0.0);
  }

  // Pirates: leader down + casualties => retreat is likely.
  {
    sim::SquadSnapshot s{};
    s.totalCount = 3;
    s.aliveCount = 2;
    s.underFireCount = 2;
    s.leaderAlive = false;
    s.avgHullFrac = 0.55;
    s.avgShieldFrac = 0.25;

    sim::PirateSquadContext ctx{};
    ctx.security01 = 0.70;
    ctx.piracy01 = 0.25;
    ctx.playerCargoValueCr = 800.0;
    ctx.committedRaid = false;

    const auto d = sim::evaluatePirateSquad(s, ctx, /*universeSeed*/ 1234u, /*groupId*/ 99u);

    CHECK(d.morale01 >= 0.0 && d.morale01 <= 1.0);
    CHECK(d.retreat);
    CHECK(d.retreatDurationSec >= 60.0);
  }

  // Police: if wanted + under fire + casualties, they should request reinforcements.
  {
    sim::SquadSnapshot s{};
    s.totalCount = 3;
    s.aliveCount = 2;
    s.underFireCount = 2;
    s.leaderAlive = true;
    s.avgHullFrac = 0.55;
    s.avgShieldFrac = 0.30;

    sim::PoliceSquadContext ctx{};
    ctx.security01 = 0.75;
    ctx.playerWanted = true;
    ctx.playerBountyCr = 9000.0;
    ctx.policeHeat = 1.0;
    ctx.fightingPlayer = true;

    const auto d = sim::evaluatePoliceSquad(s, ctx, /*universeSeed*/ 42u, /*groupId*/ 11u);

    CHECK(d.morale01 >= 0.0 && d.morale01 <= 1.0);
    CHECK(d.callReinforcements);
    CHECK(d.reinforcementUrgency01 > 0.0);
  }

  // Police: very high heat should suppress further escalation.
  {
    sim::SquadSnapshot s{};
    s.totalCount = 3;
    s.aliveCount = 2;
    s.underFireCount = 2;
    s.leaderAlive = true;
    s.avgHullFrac = 0.55;
    s.avgShieldFrac = 0.30;

    sim::PoliceSquadContext ctxLow{};
    ctxLow.security01 = 0.75;
    ctxLow.playerWanted = true;
    ctxLow.playerBountyCr = 9000.0;
    ctxLow.policeHeat = 1.0;
    ctxLow.fightingPlayer = true;

    sim::PoliceSquadContext ctxHi = ctxLow;
    ctxHi.policeHeat = 6.0;

    const auto lo = sim::evaluatePoliceSquad(s, ctxLow, /*universeSeed*/ 42u, /*groupId*/ 22u);
    const auto hi = sim::evaluatePoliceSquad(s, ctxHi, /*universeSeed*/ 42u, /*groupId*/ 22u);

    CHECK(hi.reinforcementUrgency01 <= lo.reinforcementUrgency01 + 1e-12);
  }

  return failures;
}
