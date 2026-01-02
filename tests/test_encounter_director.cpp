#include "stellar/sim/EncounterDirector.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approxEq(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

int test_encounter_director() {
  int fails = 0;

  // Deterministic seeding + basic pirate/trader cadence.
  {
    sim::EncounterDirectorState a = sim::makeEncounterDirector(1337u, 0.0);
    sim::EncounterDirectorState b = sim::makeEncounterDirector(1337u, 0.0);

    sim::EncounterDirectorContext ctx{};
    ctx.timeDays = 0.02;
    ctx.combatSimEnabled = true;
    ctx.stationCount = 2;
    ctx.localFactionId = 42;

    sim::EncounterDirectorCounts counts{};

    const int pirateCapA = sim::planPirateSpawn(a, ctx, counts);
    const int pirateCapB = sim::planPirateSpawn(b, ctx, counts);
    if (pirateCapA != 3 || pirateCapB != 3) {
      std::cerr << "[test_encounter_director] pirate cap mismatch. got="
                << pirateCapA << "/" << pirateCapB << " expected=3\n";
      ++fails;
    }
    if (!approxEq(a.nextPirateSpawnDays, b.nextPirateSpawnDays)) {
      std::cerr << "[test_encounter_director] expected deterministic pirate next time. a="
                << a.nextPirateSpawnDays << " b=" << b.nextPirateSpawnDays << "\n";
      ++fails;
    }

    const int traderA = sim::planTraderSpawn(a, ctx, counts);
    const int traderB = sim::planTraderSpawn(b, ctx, counts);
    if (traderA != 1 || traderB != 1) {
      std::cerr << "[test_encounter_director] expected trader spawn. got="
                << traderA << "/" << traderB << "\n";
      ++fails;
    }
    if (!approxEq(a.nextTraderSpawnDays, b.nextTraderSpawnDays)) {
      std::cerr << "[test_encounter_director] expected deterministic trader next time. a="
                << a.nextTraderSpawnDays << " b=" << b.nextTraderSpawnDays << "\n";
      ++fails;
    }

    const double traderDt = a.nextTraderSpawnDays - ctx.timeDays;
    if (traderDt < (70.0 / 86400.0) - 1e-12 || traderDt > (140.0 / 86400.0) + 1e-12) {
      std::cerr << "[test_encounter_director] trader next time out of bounds. dtDays="
                << traderDt << "\n";
      ++fails;
    }
  }

  // System-level knobs (piracy/security/traffic) influence cadence in the expected direction.
  {
    sim::EncounterDirectorCounts counts{};

    // Pirates: higher piracy + lower security => tighter cadence.
    {
      sim::EncounterDirectorState hi = sim::makeEncounterDirector(4242u, 0.0);
      sim::EncounterDirectorState lo = sim::makeEncounterDirector(4242u, 0.0);

      sim::EncounterDirectorContext ctxHi{};
      ctxHi.timeDays = 0.02;
      ctxHi.combatSimEnabled = true;
      ctxHi.stationCount = 2;
      ctxHi.localFactionId = 42;
      ctxHi.security01 = 0.0;
      ctxHi.piracy01 = 1.0;
      ctxHi.traffic01 = 0.5;

      sim::EncounterDirectorContext ctxLo = ctxHi;
      ctxLo.security01 = 1.0;
      ctxLo.piracy01 = 0.0;

      (void)sim::planPirateSpawn(hi, ctxHi, counts);
      (void)sim::planPirateSpawn(lo, ctxLo, counts);

      const double dtHi = hi.nextPirateSpawnDays - ctxHi.timeDays;
      const double dtLo = lo.nextPirateSpawnDays - ctxLo.timeDays;
      if (!(dtHi < dtLo)) {
        std::cerr << "[test_encounter_director] expected higher piracy to shorten pirate cadence. dtHi="
                  << dtHi << " dtLo=" << dtLo << "\n";
        ++fails;
      }
    }

    // Traders: higher traffic + lower piracy => tighter cadence.
    {
      sim::EncounterDirectorState hi = sim::makeEncounterDirector(1337u, 0.0);
      sim::EncounterDirectorState lo = sim::makeEncounterDirector(1337u, 0.0);

      sim::EncounterDirectorContext ctxHi{};
      ctxHi.timeDays = 0.02;
      ctxHi.combatSimEnabled = true;
      ctxHi.stationCount = 2;
      ctxHi.localFactionId = 42;
      ctxHi.traffic01 = 1.0;
      ctxHi.piracy01 = 0.0;

      sim::EncounterDirectorContext ctxLo = ctxHi;
      ctxLo.traffic01 = 0.0;
      ctxLo.piracy01 = 1.0;

      (void)sim::planTraderSpawn(hi, ctxHi, counts);
      (void)sim::planTraderSpawn(lo, ctxLo, counts);

      const double dtHi = hi.nextTraderSpawnDays - ctxHi.timeDays;
      const double dtLo = lo.nextTraderSpawnDays - ctxLo.timeDays;
      if (!(dtHi < dtLo)) {
        std::cerr << "[test_encounter_director] expected higher traffic to shorten trader cadence. dtHi="
                  << dtHi << " dtLo=" << dtLo << "\n";
        ++fails;
      }
    }

    // Police: higher security => tighter response times + higher desired presence.
    {
      sim::EncounterDirectorState hi = sim::makeEncounterDirector(9001u, 0.0);
      sim::EncounterDirectorState lo = sim::makeEncounterDirector(9001u, 0.0);

      sim::EncounterDirectorContext ctxHi{};
      ctxHi.timeDays = 0.02;
      ctxHi.combatSimEnabled = true;
      ctxHi.localFactionId = 12;
      ctxHi.playerWantedHere = false;
      ctxHi.localBountyCr = 0.0;
      ctxHi.policeHeat = 0.0;
      ctxHi.policeAlert = false;
      ctxHi.security01 = 1.0;
      ctxHi.traffic01 = 1.0;

      sim::EncounterDirectorContext ctxLo = ctxHi;
      ctxLo.security01 = 0.0;
      ctxLo.traffic01 = 0.0;

      sim::EncounterDirectorCounts c{};
      c.aliveTotal = 0;
      c.alivePolice = 0;
      c.alivePirates = 0;

      const auto pHi = sim::planPoliceSpawn(hi, ctxHi, c);
      const auto pLo = sim::planPoliceSpawn(lo, ctxLo, c);

      const double dtHi = hi.nextPoliceSpawnDays - ctxHi.timeDays;
      const double dtLo = lo.nextPoliceSpawnDays - ctxLo.timeDays;
      if (!(dtHi < dtLo)) {
        std::cerr << "[test_encounter_director] expected higher security to shorten police cadence. dtHi="
                  << dtHi << " dtLo=" << dtLo << "\n";
        ++fails;
      }
      if (!(pHi.desiredPolice >= pLo.desiredPolice)) {
        std::cerr << "[test_encounter_director] expected higher security/traffic to not reduce desired police. hi="
                  << pHi.desiredPolice << " lo=" << pLo.desiredPolice << "\n";
        ++fails;
      }
    }
  }

  // Desired police computation saturates/clamps.
  {
    sim::EncounterDirectorContext ctx{};
    ctx.localFactionId = 7;
    ctx.playerWantedHere = true;
    ctx.localBountyCr = 5000.0;
    ctx.policeHeat = 4.0;
    ctx.localRep = -30.0;

    sim::EncounterDirectorCounts counts{};
    counts.alivePirates = 1;

    const int desired = sim::computeDesiredPolice(ctx, counts);
    if (desired != 7) {
      std::cerr << "[test_encounter_director] desired police mismatch. got="
                << desired << " expected=7\n";
      ++fails;
    }
  }

  // Police spawn scheduling uses heat/alert bounds.
  {
    sim::EncounterDirectorState st = sim::makeEncounterDirector(9001u, 0.0);

    sim::EncounterDirectorContext ctx{};
    ctx.timeDays = 0.02;
    ctx.combatSimEnabled = true;
    ctx.localFactionId = 12;
    ctx.playerWantedHere = true;
    ctx.localBountyCr = 2000.0;
    ctx.policeHeat = 0.0;
    ctx.policeAlert = false;

    sim::EncounterDirectorCounts counts{};
    counts.aliveTotal = 0;
    counts.alivePolice = 0;
    counts.alivePirates = 1;

    const auto p = sim::planPoliceSpawn(st, ctx, counts);
    if (p.desiredPolice < 3) {
      std::cerr << "[test_encounter_director] expected desired police >= 3. got="
                << p.desiredPolice << "\n";
      ++fails;
    }
    if (p.spawnMaxCount <= 0) {
      std::cerr << "[test_encounter_director] expected police spawn request.\n";
      ++fails;
    }

    const double dt = st.nextPoliceSpawnDays - ctx.timeDays;
    if (dt < (55.0 / 86400.0) - 1e-12 || dt > (95.0 / 86400.0) + 1e-12) {
      std::cerr << "[test_encounter_director] police next time out of bounds. dtDays="
                << dt << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_encounter_director] PASS\n";
  }
  return fails;
}
