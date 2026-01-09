#include "stellar/sim/TrafficEscort.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

int test_traffic_escort() {
  int failures = 0;

  // Determinism: identical inputs should yield identical outputs.
  {
    const auto a = sim::planTrafficEscortContract(/*universeSeed=*/1234u,
                                                 /*systemId=*/55u,
                                                 /*convoyId=*/0xABCDEFu,
                                                 /*timeDays=*/42.25,
                                                 /*timeToArriveDays=*/0.08,
                                                 /*cargoValueCr=*/18000.0,
                                                 /*piracy01=*/0.65,
                                                 /*security01=*/0.55,
                                                 /*contest01=*/0.25,
                                                 /*piratesPresent=*/false);

    const auto b = sim::planTrafficEscortContract(/*universeSeed=*/1234u,
                                                 /*systemId=*/55u,
                                                 /*convoyId=*/0xABCDEFu,
                                                 /*timeDays=*/42.25,
                                                 /*timeToArriveDays=*/0.08,
                                                 /*cargoValueCr=*/18000.0,
                                                 /*piracy01=*/0.65,
                                                 /*security01=*/0.55,
                                                 /*contest01=*/0.25,
                                                 /*piratesPresent=*/false);

    CHECK(a.offer == b.offer);
    CHECK(std::abs(a.risk01 - b.risk01) < 1e-12);
    CHECK(std::abs(a.maxRangeKm - b.maxRangeKm) < 1e-9);
    CHECK(std::abs(a.durationDays - b.durationDays) < 1e-12);
    CHECK(std::abs(a.rewardCr - b.rewardCr) < 1e-9);
    CHECK(std::abs(a.bonusPerPirateCr - b.bonusPerPirateCr) < 1e-9);
    CHECK(std::abs(a.repReward - b.repReward) < 1e-12);
  }

  // Higher piracy should generally increase risk and reward.
  {
    const auto low = sim::planTrafficEscortContract(999u,
                                                   1u,
                                                   777u,
                                                   10.0,
                                                   0.10,
                                                   22000.0,
                                                   /*piracy01=*/0.10,
                                                   /*security01=*/0.60,
                                                   /*contest01=*/0.10,
                                                   /*piratesPresent=*/false);

    const auto high = sim::planTrafficEscortContract(999u,
                                                    1u,
                                                    777u,
                                                    10.0,
                                                    0.10,
                                                    22000.0,
                                                    /*piracy01=*/0.90,
                                                    /*security01=*/0.60,
                                                    /*contest01=*/0.10,
                                                    /*piratesPresent=*/false);

    CHECK(high.risk01 >= low.risk01 - 1e-12);
    CHECK(high.rewardCr >= low.rewardCr - 1e-9);
    CHECK(high.repReward >= low.repReward - 1e-12);
  }

  // Duration should never exceed the remaining time-to-arrive.
  {
    const double tta = 12.0 / 86400.0; // 12 seconds
    const auto p = sim::planTrafficEscortContract(1u,
                                                 2u,
                                                 3u,
                                                 5.0,
                                                 tta,
                                                 50000.0,
                                                 0.8,
                                                 0.2,
                                                 0.6,
                                                 /*piratesPresent=*/true);

    CHECK(p.offer); // still an offer (tta > 0)
    CHECK(p.durationDays <= tta + 1e-12);
  }

  // No time remaining => no offer.
  {
    const auto p = sim::planTrafficEscortContract(1u,
                                                 2u,
                                                 3u,
                                                 5.0,
                                                 /*timeToArriveDays=*/0.0,
                                                 50000.0,
                                                 0.8,
                                                 0.2,
                                                 0.6,
                                                 /*piratesPresent=*/true);
    CHECK(!p.offer);
    CHECK(p.durationDays == 0.0);
  }

  return failures;
}
