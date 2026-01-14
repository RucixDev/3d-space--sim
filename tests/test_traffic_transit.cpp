#include "test_harness.h"

#include "stellar/sim/TrafficTransit.h"
#include "stellar/sim/TrafficLedger.h"

#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>

using namespace stellar;
using namespace stellar::sim;

static double totalInventory(Universe& u, const StarSystem& sys, double timeDays) {
  double sum = 0.0;
  for (const auto& st : sys.stations) {
    const auto& es = u.stationEconomy(st, timeDays);
    for (double v : es.inventory) sum += v;
  }
  return sum;
}

int test_traffic_transit() {
  int failures = 0;

  // Find a system with at least 2 stations (prefer exactly 2 for determinism).
  Universe finder(1337);
  const SystemStub* stub = nullptr;
  const auto stubs = finder.queryNearby({0, 0, 0}, /*radiusLy=*/500.0, /*maxResults=*/256);
  for (const auto& s : stubs) {
    if (s.stationCount == 2) { stub = &s; break; }
  }
  if (!stub) {
    for (const auto& s : stubs) {
      if (s.stationCount >= 2) { stub = &s; break; }
    }
  }

  CHECK(stub != nullptr);

  Universe u(1337);
  const StarSystem& sys = u.getSystem(stub->id, stub);
  CHECK(sys.stations.size() >= 2);

  const double t0 = 10.5;

  // Make inventories extremely lopsided so shipments are guaranteed.
  // Even-indexed stations are full, odd-indexed stations are empty.
  for (std::size_t si = 0; si < sys.stations.size(); ++si) {
    const auto& st = sys.stations[si];
    auto& es = u.stationEconomy(st, t0);
    for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
      es.inventory[i] = (si % 2 == 0) ? st.economyModel.capacity[i] : 0.0;
    }
  }

  std::unordered_map<SystemId, int> stamps;

  TrafficLedger ledger;
  // Force shipments to remain in-flight for a while so we can assert that inventory is not credited immediately.
  ledger.params.minDurationDays = 0.6;
  ledger.params.maxDurationDays = 0.6;
  ledger.params.speedMinKmS = 0.0;
  ledger.params.speedMaxKmS = 1.0;

  const double totalBefore = totalInventory(u, sys, t0);

  simulateNpcTradeTrafficTransit(u, sys, t0, stamps, /*kMaxBackfillDays=*/14, &ledger);

  const double totalAfterPlan = totalInventory(u, sys, t0);

  // We should have generated at least one in-flight shipment.
  CHECK(!ledger.shipments.empty());

  double inTransit = 0.0;
  double maxArrive = t0;
  for (const auto& sh : ledger.shipments) {
    if (sh.systemId != sys.stub.id) continue;
    inTransit += sh.units;
    maxArrive = std::max(maxArrive, sh.arriveDay);
  }
  CHECK(inTransit > 0.0);

  // Transit-mode removes from sources immediately, but does not credit destinations until arrival.
  // Therefore total inventory should drop by exactly the sum of in-flight shipment units (no time advanced).
  const double expectedDrop = inTransit;
  const double actualDrop = totalBefore - totalAfterPlan;
  CHECK(std::abs(actualDrop - expectedDrop) <= 1e-6 * (1.0 + expectedDrop));

  // Idempotence at the same time: calling again should not create new shipments or alter totals.
  const std::size_t ledgerCount = ledger.shipments.size();
  simulateNpcTradeTrafficTransit(u, sys, t0, stamps, /*kMaxBackfillDays=*/14, &ledger);
  CHECK(ledger.shipments.size() == ledgerCount);
  const double totalAfterPlan2 = totalInventory(u, sys, t0);
  CHECK(std::abs(totalAfterPlan2 - totalAfterPlan) <= 1e-6 * (1.0 + totalAfterPlan));

  // Deliver everything: advance time past the latest arrival and block new day simulation.
  const double t1 = maxArrive + 0.25;
  stamps[sys.stub.id] = (int)std::floor(t1);

  const double totalBeforeDelivery = totalInventory(u, sys, t1);
  simulateNpcTradeTrafficTransit(u, sys, t1, stamps, /*kMaxBackfillDays=*/14, &ledger);
  const double totalAfterDelivery = totalInventory(u, sys, t1);

  // After delivery, the in-flight units should be added back into station inventories.
  const double delivered = totalAfterDelivery - totalBeforeDelivery;
  CHECK(std::abs(delivered - inTransit) <= 1e-6 * (1.0 + inTransit));

  // The original in-flight shipments should no longer exist for this system.
  bool anyRemaining = false;
  for (const auto& sh : ledger.shipments) {
    if (sh.systemId == sys.stub.id) { anyRemaining = true; break; }
  }
  CHECK(!anyRemaining);

  return failures;
}
