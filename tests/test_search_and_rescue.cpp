#include "test_harness.h"

#include "stellar/sim/SearchAndRescue.h"

#include <cmath>

using namespace stellar;

static double getRep(const sim::SaveGame& s, core::u32 factionId) {
  for (const auto& r : s.reputation) {
    if (r.factionId == factionId) return r.rep;
  }
  return 0.0;
}

int test_search_and_rescue() {
  int failures = 0;

  // ---------------------------------------------------------------------------
  // Seat accounting + addRescuedPod()
  // ---------------------------------------------------------------------------
  {
    sim::SaveGame s{};
    s.passengerSeats = 2;
    s.timeDays = 5.0;

    // One active passenger mission consumes one seat.
    sim::Mission pm{};
    pm.type = sim::MissionType::Passenger;
    pm.units = 1.0;
    pm.completed = false;
    pm.failed = false;
    s.missions.push_back(pm);

    CHECK(sim::passengerMissionSeatUsage(s) == 1);
    CHECK(sim::freePassengerSeats(s) == 1);

    sim::RescuedPod pod{};
    pod.id = 1;
    pod.recoveredSystem = 123;
    pod.registryFactionId = 3;
    pod.recoveredDay = 5.0;
    pod.lifeSupportEndDay = 6.0;
    pod.fromPlayerKill = false;

    const char* reason = nullptr;
    CHECK(sim::addRescuedPod(s, pod, s.timeDays, &reason));
    CHECK(reason == nullptr);
    CHECK(sim::rescuedPodSeatUsage(s) == 1);
    CHECK(sim::occupiedPassengerSeatsTotal(s) == 2);
    CHECK(sim::freePassengerSeats(s) == 0);

    // No free seats left.
    sim::RescuedPod pod2 = pod;
    pod2.id = 2;
    reason = nullptr;
    CHECK(!sim::addRescuedPod(s, pod2, s.timeDays, &reason));
    CHECK(reason != nullptr);

    // Even without authority, we can disembark the pods (no reward).
    const double creditsBefore = s.credits;
    const auto r = sim::applyRescuePodTurnIn(s, s.timeDays, /*stationFactionId=*/3, /*hasAuthority=*/false);
    CHECK(r.ok);
    CHECK(r.creditsPaid == 0.0);
    CHECK(r.repAwarded == 0.0);
    CHECK(s.rescuedPods.empty());
    CHECK(s.credits == creditsBefore);
  }

  // ---------------------------------------------------------------------------
  // Turn-in payouts + reputation application
  // ---------------------------------------------------------------------------
  {
    sim::SaveGame s{};
    s.passengerSeats = 8;
    s.timeDays = 10.0;
    s.credits = 100.0;

    // Alive, rewardable.
    sim::RescuedPod a{};
    a.id = 10;
    a.recoveredSystem = 111;
    a.registryFactionId = 7;
    a.recoveredDay = 10.0;
    a.lifeSupportEndDay = 11.0;
    a.fromPlayerKill = false;

    // Expired, rewardable.
    sim::RescuedPod b{};
    b.id = 11;
    b.recoveredSystem = 111;
    b.registryFactionId = 7;
    b.recoveredDay = 9.0;
    b.lifeSupportEndDay = 9.5;
    b.fromPlayerKill = false;

    // No reward (anti-farm).
    sim::RescuedPod c{};
    c.id = 12;
    c.recoveredSystem = 111;
    c.registryFactionId = 7;
    c.recoveredDay = 10.0;
    c.lifeSupportEndDay = 11.0;
    c.fromPlayerKill = true;

    // Independent pod (registryFactionId==0); should bucket to station faction.
    sim::RescuedPod d{};
    d.id = 13;
    d.recoveredSystem = 222;
    d.registryFactionId = 0;
    d.recoveredDay = 10.0;
    d.lifeSupportEndDay = 10.8;
    d.fromPlayerKill = false;

    s.rescuedPods = {a, b, c, d};

    const auto q = sim::quoteRescuePodTurnIn(s, s.timeDays, /*stationFactionId=*/7, /*hasAuthority=*/true);
    CHECK(q.ok);
    CHECK(q.totalPods == 4);
    CHECK(q.noRewardPods == 1);
    CHECK(q.expiredPods == 1);
    CHECK(q.rewardablePods == 3);
    CHECK(q.creditsTotal > 0.0);
    CHECK(q.repTotal > 0.0);

    const double creditsBefore = s.credits;
    const double repBefore = getRep(s, 7);

    const auto r = sim::applyRescuePodTurnIn(s, s.timeDays, /*stationFactionId=*/7, /*hasAuthority=*/true);
    CHECK(r.ok);
    CHECK(r.turnedIn == 4);
    CHECK(std::abs(r.creditsPaid - q.creditsTotal) < 1e-9);
    CHECK(std::abs(r.repAwarded - q.repTotal) < 1e-9);
    CHECK(s.rescuedPods.empty());
    CHECK(s.credits > creditsBefore);

    const double repAfter = getRep(s, 7);
    CHECK(repAfter > repBefore);
    CHECK(std::abs((repAfter - repBefore) - q.repTotal) < 1e-6);
  }

  return failures;
}
