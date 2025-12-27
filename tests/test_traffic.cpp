#include "stellar/sim/Traffic.h"

#include "stellar/sim/Universe.h"

#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

static bool nearly(double a, double b, double eps = 1e-8) {
  return std::abs(a - b) <= eps;
}

static bool inventoriesEqual(const std::vector<std::array<double, stellar::econ::kCommodityCount>>& a,
                             const std::vector<std::array<double, stellar::econ::kCommodityCount>>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t si = 0; si < a.size(); ++si) {
    for (std::size_t i = 0; i < stellar::econ::kCommodityCount; ++i) {
      if (!nearly(a[si][i], b[si][i])) return false;
    }
  }
  return true;
}

static std::vector<std::array<double, stellar::econ::kCommodityCount>> snapshotInventories(stellar::sim::Universe& u,
                                                                                           const stellar::sim::StarSystem& sys,
                                                                                           double timeDays) {
  std::vector<std::array<double, stellar::econ::kCommodityCount>> out;
  out.reserve(sys.stations.size());
  for (const auto& st : sys.stations) {
    auto& es = u.stationEconomy(st, timeDays);
    out.push_back(es.inventory);
  }
  return out;
}

int test_traffic() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  // Find a system with at least 2 stations (traffic needs at least two endpoints).
  Universe finder(1337);
  const auto stubs = finder.queryNearby({0, 0, 0}, 500.0, 256);
  const SystemStub* chosen = nullptr;
  for (const auto& s : stubs) {
    if (s.stationCount >= 2) { chosen = &s; break; }
  }
  if (!chosen) {
    std::cerr << "[test_traffic] no systems with >=2 stations found in query (seed may have changed?)\n";
    return 1;
  }

  const double t0 = 10.0;
  const int day0 = (int)std::floor(t0);

  // Map-based stamps.
  {
    Universe u(1337);
    const auto& sys = u.getSystem(chosen->id, chosen);
    if (sys.stations.size() < 2) {
      std::cerr << "[test_traffic] stub said >=2 stations but system had <2\n";
      return 1;
    }

    std::unordered_map<SystemId, int> stamps;
    simulateNpcTradeTraffic(u, sys, t0, stamps);

    const auto it = stamps.find(sys.stub.id);
    if (it == stamps.end() || it->second != day0) {
      std::cerr << "[test_traffic] stamps map did not update dayStamp to current day\n";
      ++fails;
    }

    const auto inv1 = snapshotInventories(u, sys, t0);

    // Idempotent when called again with the same day.
    simulateNpcTradeTraffic(u, sys, t0, stamps);
    const auto inv2 = snapshotInventories(u, sys, t0);
    if (!inventoriesEqual(inv1, inv2)) {
      std::cerr << "[test_traffic] repeated simulateNpcTradeTraffic changed inventories for same day\n";
      ++fails;
    }

    // Deterministic across fresh universes (same seed + same stamp rules).
    Universe u2(1337);
    const auto& sys2 = u2.getSystem(chosen->id);
    std::unordered_map<SystemId, int> stamps2;
    simulateNpcTradeTraffic(u2, sys2, t0, stamps2);
    const auto inv3 = snapshotInventories(u2, sys2, t0);
    if (!inventoriesEqual(inv1, inv3)) {
      std::cerr << "[test_traffic] inventories differ across identical runs (map stamps)\n";
      ++fails;
    }
  }

  // SaveGame-friendly vector stamps overload should match the map behavior.
  {
    Universe u(1337);
    const auto& sys = u.getSystem(chosen->id);

    std::vector<SystemTrafficStamp> stamps;
    simulateNpcTradeTraffic(u, sys, t0, stamps);

    if (stamps.size() != 1 || stamps.front().systemId != sys.stub.id || stamps.front().dayStamp != day0) {
      std::cerr << "[test_traffic] stamps vector did not contain expected {systemId, dayStamp}\n";
      ++fails;
    }

    const auto inv1 = snapshotInventories(u, sys, t0);

    // Same-day call should be idempotent.
    simulateNpcTradeTraffic(u, sys, t0, stamps);
    const auto inv2 = snapshotInventories(u, sys, t0);
    if (!inventoriesEqual(inv1, inv2)) {
      std::cerr << "[test_traffic] repeated simulateNpcTradeTraffic changed inventories for same day (vector stamps)\n";
      ++fails;
    }

    // Ensure advancing time updates the stamp.
    const double t1 = 12.0;
    const int day1 = (int)std::floor(t1);
    simulateNpcTradeTraffic(u, sys, t1, stamps);
    if (stamps.empty() || stamps.front().dayStamp != day1) {
      std::cerr << "[test_traffic] stamps vector did not advance to the new day\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_traffic] pass\n";
  return fails;
}
