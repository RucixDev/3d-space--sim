#include "stellar/sim/DockingClearanceService.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approxEq(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

int test_docking_clearance() {
  int fails = 0;

  sim::Station st{};
  st.id = 123;
  st.name = "Test Station";
  st.type = econ::StationType::Outpost;
  st.commsRangeKm = 1000.0;

  const core::u64 seed = 1337u;

  // Out of range: no request attempt should be recorded.
  {
    sim::DockingClearanceState cs{};
    sim::DockingClearanceParams params{};
    params.minGrantProb = 0.0;
    params.maxGrantProb = 1.0;

    const auto d = sim::requestDockingClearance(seed, st, 10.0, 2000.0, cs, 0.0, params);
    if (d.status != sim::DockingClearanceStatus::OutOfRange) {
      std::cerr << "[test_docking_clearance] expected OutOfRange\n";
      ++fails;
    }
    if (cs.requestCount != 0) {
      std::cerr << "[test_docking_clearance] expected requestCount unchanged for OutOfRange. got="
                << cs.requestCount << "\n";
      ++fails;
    }
  }

  // Force grant: clamp the probability to 1.0.
  {
    sim::DockingClearanceState cs{};
    sim::DockingClearanceParams params{};
    params.baseGrantProb = 1.0;
    params.minGrantProb = 1.0;
    params.maxGrantProb = 1.0;
    params.trafficPenalty = 0.0;

    const double t = 1.0;
    const auto d = sim::requestDockingClearance(seed, st, t, 500.0, cs, 0.0, params);
    if (d.status != sim::DockingClearanceStatus::Granted || !d.hasClearance) {
      std::cerr << "[test_docking_clearance] expected Granted. got status="
                << (int)d.status << " hasClearance=" << (int)d.hasClearance << "\n";
      ++fails;
    }
    const double expectedExpires = t + params.grantDurationSec / 86400.0;
    if (!approxEq(cs.expiresDays, expectedExpires)) {
      std::cerr << "[test_docking_clearance] expiresDays mismatch. got="
                << cs.expiresDays << " expected=" << expectedExpires << "\n";
      ++fails;
    }
    if (!sim::dockingClearanceValid(cs, t)) {
      std::cerr << "[test_docking_clearance] expected clearance valid immediately after grant\n";
      ++fails;
    }

    // Calling again while valid should be a no-op unless allowRefresh=true.
    const core::u32 rc = cs.requestCount;
    const auto d2 = sim::requestDockingClearance(seed, st, t + 0.00001, 500.0, cs, 0.0, params);
    if (d2.status != sim::DockingClearanceStatus::Granted || !d2.hasClearance) {
      std::cerr << "[test_docking_clearance] expected Granted (pre-existing clearance)\n";
      ++fails;
    }
    if (cs.requestCount != rc) {
      std::cerr << "[test_docking_clearance] expected requestCount unchanged when already granted. got="
                << cs.requestCount << " expected=" << rc << "\n";
      ++fails;
    }
  }

  // Force deny and verify cooldown/throttle.
  {
    sim::DockingClearanceState cs{};
    sim::DockingClearanceParams params{};
    params.baseGrantProb = 0.0;
    params.minGrantProb = 0.0;
    params.maxGrantProb = 0.0;

    const double t = 2.0;
    const auto d = sim::requestDockingClearance(seed, st, t, 500.0, cs, 0.0, params);
    if (d.status != sim::DockingClearanceStatus::Denied || d.hasClearance) {
      std::cerr << "[test_docking_clearance] expected Denied. got status="
                << (int)d.status << " hasClearance=" << (int)d.hasClearance << "\n";
      ++fails;
    }
    if (cs.cooldownUntilDays <= t) {
      std::cerr << "[test_docking_clearance] expected cooldownUntilDays in the future. got="
                << cs.cooldownUntilDays << " now=" << t << "\n";
      ++fails;
    }

    // Re-request during cooldown should be throttled and should not consume a requestCount.
    const core::u32 rc = cs.requestCount;
    const double mid = t + (params.denyCooldownSec / 86400.0) * 0.5;
    const auto d2 = sim::requestDockingClearance(seed, st, mid, 500.0, cs, 0.0, params);
    if (d2.status != sim::DockingClearanceStatus::Throttled) {
      std::cerr << "[test_docking_clearance] expected Throttled during cooldown. got="
                << (int)d2.status << "\n";
      ++fails;
    }
    if (cs.requestCount != rc) {
      std::cerr << "[test_docking_clearance] expected requestCount unchanged during throttle. got="
                << cs.requestCount << " expected=" << rc << "\n";
      ++fails;
    }
  }

  // Determinism: same seed + station + request order yields the same decision.
  {
    sim::DockingClearanceState a{};
    sim::DockingClearanceState b{};

    sim::DockingClearanceParams params{};
    params.minGrantProb = 0.0;
    params.maxGrantProb = 1.0;

    const double t = 5.0;
    const double traffic = 0.25;

    const auto da = sim::requestDockingClearance(seed, st, t, 500.0, a, traffic, params);
    const auto db = sim::requestDockingClearance(seed, st, t, 500.0, b, traffic, params);

    if (da.status != db.status || a.granted != b.granted || !approxEq(da.pGrant, db.pGrant, 1e-12)) {
      std::cerr << "[test_docking_clearance] expected deterministic decisions. aStatus="
                << (int)da.status << " bStatus=" << (int)db.status << " aGranted="
                << (int)a.granted << " bGranted=" << (int)b.granted << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_docking_clearance] PASS\n";
  }
  return fails;
}
