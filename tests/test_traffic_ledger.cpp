#include "stellar/sim/Traffic.h"
#include "stellar/sim/TrafficLedger.h"
#include "stellar/sim/Units.h"
#include "stellar/sim/WorldIds.h"

#include "stellar/sim/Universe.h"

#include <cmath>
#include <iostream>
#include <unordered_map>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_traffic_ledger() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  const core::u64 seed = 777777ull;

  // Pick a system with stations (and ideally >2) for more interesting traffic.
  Universe u(seed);
  const auto stubs = u.queryNearby({0,0,0}, 60.0, 16);
  if (stubs.empty()) {
    std::cerr << "[test_traffic_ledger] universe returned no stubs\n";
    return 1;
  }

  const StarSystem* sys = nullptr;
  for (const auto& stub : stubs) {
    const auto& s = u.getSystem(stub.id, &stub);
    if (s.stations.size() >= 2) {
      sys = &s;
      if (s.stations.size() >= 3) break;
    }
  }
  if (!sys) {
    std::cerr << "[test_traffic_ledger] expected a system with >=2 stations\n";
    return 1;
  }

  const double t0 = 1200.25;

  std::unordered_map<SystemId, int> stamps;
  TrafficLedger ledger;
  simulateNpcTradeTraffic(u, *sys, t0, stamps, /*kMaxBackfillDays=*/14, &ledger);

  // Repeating the same day should not add extra shipments.
  const std::size_t n0 = ledger.shipments.size();
  simulateNpcTradeTraffic(u, *sys, t0, stamps, /*kMaxBackfillDays=*/14, &ledger);
  if (ledger.shipments.size() != n0) {
    std::cerr << "[test_traffic_ledger] repeated simulateNpcTradeTraffic appended shipments\n";
    ++fails;
  }

  // Basic sanity check of schedule metadata.
  for (const auto& s : ledger.shipments) {
    if (!isDeterministicWorldId(s.id)) {
      std::cerr << "[test_traffic_ledger] shipment id is not tagged deterministic\n";
      ++fails;
      break;
    }
    if (s.systemId != sys->stub.id) {
      std::cerr << "[test_traffic_ledger] shipment systemId mismatch\n";
      ++fails;
      break;
    }
    if (s.fromStation == 0 || s.toStation == 0 || s.fromStation == s.toStation) {
      std::cerr << "[test_traffic_ledger] invalid station ids\n";
      ++fails;
      break;
    }
    if (!(s.departDay >= (double)s.dayStamp && s.departDay < (double)s.dayStamp + 1.0)) {
      std::cerr << "[test_traffic_ledger] departDay not within dayStamp\n";
      ++fails;
      break;
    }
    if (s.arriveDay + 1e-12 < s.departDay) {
      std::cerr << "[test_traffic_ledger] arriveDay < departDay\n";
      ++fails;
      break;
    }
    if (s.distKm < -1e-6) {
      std::cerr << "[test_traffic_ledger] negative distKm\n";
      ++fails;
      break;
    }
    if (s.units < -1e-6) {
      std::cerr << "[test_traffic_ledger] negative units\n";
      ++fails;
      break;
    }

    // Schedule metric consistency: distKm/speedKmS should match endpoints + timeline.
    const Station* fromSt = nullptr;
    const Station* toSt = nullptr;
    for (const auto& st : sys->stations) {
      if (st.id == s.fromStation) fromSt = &st;
      if (st.id == s.toStation) toSt = &st;
    }
    if (fromSt && toSt) {
      const double dist = (stationPosKm(*toSt, s.arriveDay) - stationPosKm(*fromSt, s.departDay)).length();
      const double durSec = std::max(1e-12, (s.arriveDay - s.departDay) * kSecondsPerDay);
      const double spd = dist / durSec;

      const double distErr = std::abs(dist - s.distKm);
      const double distTol = std::max(1e-3, dist * 1e-9);
      if (distErr > distTol) {
        std::cerr << "[test_traffic_ledger] distKm mismatch vs endpoints\n";
        ++fails;
        break;
      }

      const double spdErr = std::abs(spd - s.speedKmS);
      const double spdTol = std::max(1e-3, spd * 1e-9);
      if (spdErr > spdTol) {
        std::cerr << "[test_traffic_ledger] speedKmS mismatch vs endpoints\n";
        ++fails;
        break;
      }
    }
  }

  // Determinism check: same seed/system/time should produce identical shipment logs.
  Universe u2(seed);
  const auto& sys2 = u2.getSystem(sys->stub.id, &sys->stub);

  std::unordered_map<SystemId, int> stamps2;
  TrafficLedger ledger2;
  simulateNpcTradeTraffic(u2, sys2, t0, stamps2, /*kMaxBackfillDays=*/14, &ledger2);

  if (ledger.shipments.size() != ledger2.shipments.size()) {
    std::cerr << "[test_traffic_ledger] shipment count mismatch across identical seeds\n";
    ++fails;
  }

  const std::size_t n = std::min(ledger.shipments.size(), ledger2.shipments.size());
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = ledger.shipments[i];
    const auto& b = ledger2.shipments[i];

    if (a.id != b.id || a.dayStamp != b.dayStamp || a.fromStation != b.fromStation || a.toStation != b.toStation || a.commodity != b.commodity) {
      std::cerr << "[test_traffic_ledger] shipment identity mismatch\n";
      ++fails;
      break;
    }
    if (!nearly(a.units, b.units, 1e-9)) {
      std::cerr << "[test_traffic_ledger] shipment units mismatch\n";
      ++fails;
      break;
    }
    if (!nearly(a.departDay, b.departDay, 1e-12) || !nearly(a.arriveDay, b.arriveDay, 1e-12)) {
      std::cerr << "[test_traffic_ledger] shipment schedule mismatch\n";
      ++fails;
      break;
    }
    if (!nearly(a.distKm, b.distKm, 1e-6) || !nearly(a.speedKmS, b.speedKmS, 1e-6)) {
      std::cerr << "[test_traffic_ledger] shipment schedule metrics mismatch\n";
      ++fails;
      break;
    }
  }

  if (fails == 0) {
    std::cout << "[test_traffic_ledger] PASS\n";
  }
  return fails;
}