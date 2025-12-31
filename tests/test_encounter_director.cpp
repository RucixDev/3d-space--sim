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
