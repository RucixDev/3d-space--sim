#include "stellar/sim/Derelict.h"

#include <cmath>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

int test_derelict() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  const core::u64 seed = 424242ull;
  const SystemId sys = 12345ull;
  const core::u64 sig = 111222333ull;
  const double t = 1000.25;

  // Determinism.
  {
    const DerelictPlan a = planDerelictEncounter(seed, sys, sig, t,
                                                 /*piracy01=*/0.62,
                                                 /*security01=*/0.41,
                                                 /*contest01=*/0.22,
                                                 /*missionSite=*/false,
                                                 /*includeDayStamp=*/true);
    const DerelictPlan b = planDerelictEncounter(seed, sys, sig, t,
                                                 /*piracy01=*/0.62,
                                                 /*security01=*/0.41,
                                                 /*contest01=*/0.22,
                                                 /*missionSite=*/false,
                                                 /*includeDayStamp=*/true);

    if (a.scenario != b.scenario || a.wreckClass != b.wreckClass || a.ambush != b.ambush || a.pirateCount != b.pirateCount) {
      std::cerr << "[test_derelict] plan not deterministic (flags)\n";
      ++fails;
    }
    if (a.salvageCommodity != b.salvageCommodity) {
      std::cerr << "[test_derelict] plan not deterministic (commodity)\n";
      ++fails;
    }
    if (!nearly(a.salvageUnits, b.salvageUnits) || a.salvagePods != b.salvagePods ||
        a.hasDataCore != b.hasDataCore || !nearly(a.dataUnits, b.dataUnits) ||
        !nearly(a.risk, b.risk)) {
      std::cerr << "[test_derelict] plan not deterministic (numeric)\n";
      ++fails;
    }
  }

  // Basic bounds/invariants.
  {
    const DerelictPlan p = planDerelictEncounter(seed, sys, sig, t,
                                                 /*piracy01=*/0.62,
                                                 /*security01=*/0.41,
                                                 /*contest01=*/0.22,
                                                 /*missionSite=*/false,
                                                 /*includeDayStamp=*/true);

    if (p.risk < 0.0 || p.risk > 1.0) {
      std::cerr << "[test_derelict] risk out of bounds\n";
      ++fails;
    }
    if (p.wreckClass > 2) {
      std::cerr << "[test_derelict] wreckClass out of bounds\n";
      ++fails;
    }
    if (p.hasSalvage) {
      if (p.salvageUnits <= 1e-6) {
        std::cerr << "[test_derelict] hasSalvage requires positive salvageUnits\n";
        ++fails;
      }
      if (p.salvagePods < 1) {
        std::cerr << "[test_derelict] hasSalvage requires salvagePods >= 1\n";
        ++fails;
      }
    } else {
      if (p.salvagePods != 0) {
        std::cerr << "[test_derelict] no-salvage plan should not have pods\n";
        ++fails;
      }
    }
    if (p.ambush) {
      if (p.pirateCount < 2) {
        std::cerr << "[test_derelict] ambush should spawn at least 2 pirates\n";
        ++fails;
      }
    } else {
      if (p.pirateCount != 0) {
        std::cerr << "[test_derelict] non-ambush should not have pirates\n";
        ++fails;
      }
    }
    if (!p.hasDataCore && p.dataUnits > 1e-9) {
      std::cerr << "[test_derelict] no data core should have 0 dataUnits\n";
      ++fails;
    }
  }

  // Mission sites should remain stable across days when includeDayStamp=false.
  {
    const DerelictPlan d0 = planDerelictEncounter(seed, sys, sig, 1000.1,
                                                  /*piracy01=*/0.70,
                                                  /*security01=*/0.30,
                                                  /*contest01=*/0.35,
                                                  /*missionSite=*/true,
                                                  /*includeDayStamp=*/false);
    const DerelictPlan d1 = planDerelictEncounter(seed, sys, sig, 1001.9,
                                                  /*piracy01=*/0.70,
                                                  /*security01=*/0.30,
                                                  /*contest01=*/0.35,
                                                  /*missionSite=*/true,
                                                  /*includeDayStamp=*/false);
    const bool same = (d0.scenario == d1.scenario && d0.wreckClass == d1.wreckClass && d0.ambush == d1.ambush &&
                       d0.pirateCount == d1.pirateCount && d0.salvageCommodity == d1.salvageCommodity &&
                       nearly(d0.salvageUnits, d1.salvageUnits) && d0.salvagePods == d1.salvagePods &&
                       d0.hasDataCore == d1.hasDataCore && nearly(d0.dataUnits, d1.dataUnits));
    if (!same) {
      std::cerr << "[test_derelict] mission derelict should not change across days\n";
      ++fails;
    }
  }

  return fails;
}
