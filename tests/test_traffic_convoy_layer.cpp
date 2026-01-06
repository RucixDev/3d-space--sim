#include "stellar/sim/Traffic.h"
#include "stellar/sim/TrafficConvoyLayer.h"
#include "stellar/sim/Units.h"
#include "stellar/sim/WorldIds.h"

#include "stellar/sim/Universe.h"

#include <cmath>
#include <iostream>
#include <unordered_map>

static bool nearly(double a, double b, double eps = 1e-7) {
  return std::abs(a - b) <= eps;
}

static bool nearlyVec(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps = 1e-3) {
  return nearly(a.x, b.x, eps) && nearly(a.y, b.y, eps) && nearly(a.z, b.z, eps);
}

static const stellar::sim::Station* findStation(const stellar::sim::StarSystem& sys, stellar::sim::StationId id) {
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

static double perpSpeedKmS(const stellar::sim::TrafficConvoyState& st) {
  // Magnitude of velocity component perpendicular to the lane direction.
  const double vpar = stellar::math::dot(st.velKmS, st.dir);
  const stellar::math::Vec3d vperp = st.velKmS - st.dir * vpar;
  return vperp.length();
}

int test_traffic_convoy_layer() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  const core::u64 seed = 424242ull;

  // Pick a system with at least 2 stations.
  Universe u(seed);
  const auto stubs = u.queryNearby({0, 0, 0}, 80.0, 32);
  if (stubs.empty()) {
    std::cerr << "[test_traffic_convoy_layer] no stubs returned\n";
    return 1;
  }

  const SystemStub* chosen = nullptr;
  for (const auto& s : stubs) {
    if (s.stationCount >= 2) { chosen = &s; break; }
  }
  if (!chosen) {
    std::cerr << "[test_traffic_convoy_layer] expected at least one system with >=2 stations\n";
    return 1;
  }

  const auto& sys = u.getSystem(chosen->id, chosen);
  if (sys.stations.size() < 2) {
    std::cerr << "[test_traffic_convoy_layer] chosen system generated with <2 stations\n";
    return 1;
  }

  const double t = 900.35;

  std::unordered_map<SystemId, int> stamps;
  TrafficLedger ledger;
  simulateNpcTradeTraffic(u, sys, t, stamps, /*kMaxBackfillDays=*/14, &ledger);
  ledger.prune(t);

  TrafficLaneParams lane;
  lane.includeInactive = true;

  const auto c1 = generateTrafficConvoysFromLedger(ledger, sys, t, /*windowDays=*/1, /*includeInactive=*/true, lane);
  const auto c2 = generateTrafficConvoysFromLedger(ledger, sys, t, /*windowDays=*/1, /*includeInactive=*/true, lane);

  if (c1.size() != c2.size()) {
    std::cerr << "[test_traffic_convoy_layer] non-deterministic output within same run\n";
    ++fails;
  }

  const std::size_t n = std::min(c1.size(), c2.size());
  for (std::size_t i = 0; i < n; ++i) {
    if (c1[i].convoy.id != c2[i].convoy.id) {
      std::cerr << "[test_traffic_convoy_layer] id mismatch at index " << i << "\n";
      ++fails;
      break;
    }
    if (!isDeterministicWorldId(c1[i].convoy.id)) {
      std::cerr << "[test_traffic_convoy_layer] convoy id missing deterministic bit\n";
      ++fails;
      break;
    }
  }

  // Endpoint correctness: arc offset should be 0 at depart/arrive.
  if (!c1.empty()) {
    const auto& c = c1.front().convoy;
    const Station* from = findStation(sys, c.fromStation);
    const Station* to = findStation(sys, c.toStation);

    if (from && to) {
      const auto s0 = evaluateTrafficConvoy(c, sys, c.departDay, lane);
      const auto s1 = evaluateTrafficConvoy(c, sys, c.arriveDay, lane);

      const auto p0 = stationPosKm(*from, c.departDay);
      const auto p1 = stationPosKm(*to, c.arriveDay);

      if (!nearlyVec(s0.posKm, p0, 1e-3)) {
        std::cerr << "[test_traffic_convoy_layer] depart endpoint mismatch\n";
        ++fails;
      }
      if (!nearlyVec(s1.posKm, p1, 1e-3)) {
        std::cerr << "[test_traffic_convoy_layer] arrive endpoint mismatch\n";
        ++fails;
      }

      // Arc easing should also yield ~zero lateral velocity at endpoints.
      const double perp0 = perpSpeedKmS(s0);
      const double perp1 = perpSpeedKmS(s1);
      if (perp0 > 1e-6) {
        std::cerr << "[test_traffic_convoy_layer] non-zero lateral speed at depart endpoint\n";
        ++fails;
      }
      if (perp1 > 1e-6) {
        std::cerr << "[test_traffic_convoy_layer] non-zero lateral speed at arrive endpoint\n";
        ++fails;
      }
    }
  }


  // Schedule hydration: convoy replay should still work when a ledger shipment is missing
  // schedule metadata (common when loading older/corrupted saves).
  {
    TrafficLedger ledger2;
    ledger2.params = ledger.params;

    TrafficShipment sh{};
    sh.id = makeDeterministicWorldId(seed, 0xABCDEF123456ull);
    sh.systemId = sys.stub.id;
    sh.dayStamp = (int)std::floor(t);
    sh.fromStation = sys.stations[0].id;
    sh.toStation = sys.stations[1].id;
    sh.factionId = sys.stations[0].factionId;
    sh.commodity = econ::CommodityId::Food;
    sh.units = 10.0;
    // Intentionally leave departDay/arriveDay/distKm/speedKmS at 0.

    ledger2.record(sh);

    const auto v1 = generateTrafficConvoysFromLedger(ledger2, sys, t, /*windowDays=*/0, /*includeInactive=*/true, lane);
    const auto v2 = generateTrafficConvoysFromLedger(ledger2, sys, t, /*windowDays=*/0, /*includeInactive=*/true, lane);

    if (v1.empty()) {
      std::cerr << "[test_traffic_convoy_layer] expected hydrated convoy from missing schedule metadata\n";
      ++fails;
    } else {
      const auto& c = v1.front().convoy;
      if (c.id != sh.id) {
        std::cerr << "[test_traffic_convoy_layer] hydrated convoy id mismatch\n";
        ++fails;
      }
      if (!isDeterministicWorldId(c.id)) {
        std::cerr << "[test_traffic_convoy_layer] hydrated convoy id missing deterministic bit\n";
        ++fails;
      }
      if (!v2.empty()) {
        if (!nearly(c.departDay, v2.front().convoy.departDay, 1e-9) || !nearly(c.arriveDay, v2.front().convoy.arriveDay, 1e-9)) {
          std::cerr << "[test_traffic_convoy_layer] schedule hydration is not deterministic\n";
          ++fails;
        }
      }

      const double day = (double)sh.dayStamp;
      if (c.departDay < day - 1e-9 || c.departDay > day + 1.0 + 1e-9) {
        std::cerr << "[test_traffic_convoy_layer] hydrated departDay not within dayStamp window\n";
        ++fails;
      }
      if (!(c.arriveDay > c.departDay)) {
        std::cerr << "[test_traffic_convoy_layer] hydrated arriveDay not after departDay\n";
        ++fails;
      }

      // Now query at mid-flight with includeInactive=false and ensure we see it as active.
      const double midT = 0.5 * (c.departDay + c.arriveDay);
      const auto mid = generateTrafficConvoysFromLedger(ledger2, sys, midT, /*windowDays=*/1, /*includeInactive=*/false, lane);
      bool found = false;
      for (const auto& v : mid) {
        if (v.convoy.id != c.id) continue;
        found = true;
        if (!v.state.active) {
          std::cerr << "[test_traffic_convoy_layer] hydrated convoy should be active at midT\n";
          ++fails;
        }
        break;
      }
      if (!found) {
        std::cerr << "[test_traffic_convoy_layer] hydrated convoy missing from mid-flight query\n";
        ++fails;
      }
    }
  }

  if (fails == 0) std::cout << "[test_traffic_convoy_layer] PASS\n";
  return fails;
}
