#include "stellar/sim/TrafficLanes.h"

#include "stellar/sim/Units.h"
#include "stellar/sim/Universe.h"
#include "stellar/sim/WorldIds.h"

#include <cmath>
#include <iostream>
#include <unordered_set>

static bool nearly(double a, double b, double eps = 1e-7) {
  return std::abs(a - b) <= eps;
}

static bool nearlyVec(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps = 1e-6) {
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

int test_traffic_lanes() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  // Find a system with at least 2 stations.
  Universe finder(1337);
  const auto stubs = finder.queryNearby({0, 0, 0}, 500.0, 256);
  const SystemStub* chosen = nullptr;
  for (const auto& s : stubs) {
    if (s.stationCount >= 2) { chosen = &s; break; }
  }
  if (!chosen) {
    std::cerr << "[test_traffic_lanes] no systems with >=2 stations found in query\n";
    return 1;
  }

  const double t = 10.5;
  TrafficLaneParams p;
  p.includeInactive = true;
  p.genWindowDays = 0; // keep comparisons simple

  Universe u(1337);
  const auto& sys = u.getSystem(chosen->id, chosen);

  const auto convoys1 = generateTrafficConvoys(u.seed(), sys, t, p);
  const auto convoys2 = generateTrafficConvoys(u.seed(), sys, t, p);

  if (convoys1.size() != convoys2.size()) {
    std::cerr << "[test_traffic_lanes] generateTrafficConvoys not deterministic within same run\n";
    ++fails;
  }

  const std::size_t n = std::min(convoys1.size(), convoys2.size());
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = convoys1[i];
    const auto& b = convoys2[i];

    if (a.convoy.id != b.convoy.id) {
      std::cerr << "[test_traffic_lanes] id mismatch at index " << i << "\n";
      ++fails;
      break;
    }

    if (!isDeterministicWorldId(a.convoy.id)) {
      std::cerr << "[test_traffic_lanes] convoy id missing deterministic bit\n";
      ++fails;
    }

    if (a.convoy.fromStation == a.convoy.toStation) {
      std::cerr << "[test_traffic_lanes] fromStation == toStation\n";
      ++fails;
    }
    if (!(a.convoy.departDay < a.convoy.arriveDay)) {
      std::cerr << "[test_traffic_lanes] departDay not < arriveDay\n";
      ++fails;
    }
    if (!(a.convoy.units > 0.0)) {
      std::cerr << "[test_traffic_lanes] non-positive units\n";
      ++fails;
    }

    if (a.convoy.fromStation != b.convoy.fromStation || a.convoy.toStation != b.convoy.toStation ||
        a.convoy.commodity != b.convoy.commodity || !nearly(a.convoy.units, b.convoy.units) ||
        !nearly(a.convoy.departDay, b.convoy.departDay) || !nearly(a.convoy.arriveDay, b.convoy.arriveDay) ||
        a.state.active != b.state.active || !nearly(a.state.progress01, b.state.progress01) ||
        !nearlyVec(a.state.posKm, b.state.posKm) || !nearlyVec(a.state.velKmS, b.state.velKmS)) {
      std::cerr << "[test_traffic_lanes] schedule/state mismatch at index " << i << "\n";
      ++fails;
      break;
    }
  }

  // Deterministic across identical universes.
  {
    Universe u2(1337);
    const auto& sys2 = u2.getSystem(chosen->id);
    const auto convoys3 = generateTrafficConvoys(u2.seed(), sys2, t, p);
    if (convoys1.size() != convoys3.size()) {
      std::cerr << "[test_traffic_lanes] convoy count differs across identical universes\n";
      ++fails;
    }
    const std::size_t m = std::min(convoys1.size(), convoys3.size());
    for (std::size_t i = 0; i < m; ++i) {
      if (convoys1[i].convoy.id != convoys3[i].convoy.id) {
        std::cerr << "[test_traffic_lanes] id differs across identical universes\n";
        ++fails;
        break;
      }
      if (!nearlyVec(convoys1[i].state.posKm, convoys3[i].state.posKm) ||
          !nearlyVec(convoys1[i].state.velKmS, convoys3[i].state.velKmS)) {
        std::cerr << "[test_traffic_lanes] state differs across identical universes\n";
        ++fails;
        break;
      }
    }
  }

  // Endpoint correctness: the lane arc should be 0 at t=0 and t=1.
  if (!convoys1.empty()) {
    const auto& c = convoys1.front().convoy;
    const Station* from = findStation(sys, c.fromStation);
    const Station* to = findStation(sys, c.toStation);
    if (from && to) {
      const auto s0 = evaluateTrafficConvoy(c, sys, c.departDay, p);
      const auto s1 = evaluateTrafficConvoy(c, sys, c.arriveDay, p);
      const auto p0 = stationPosKm(*from, c.departDay);
      const auto p1 = stationPosKm(*to, c.arriveDay);
      if (!nearlyVec(s0.posKm, p0, 1e-3)) {
        std::cerr << "[test_traffic_lanes] depart endpoint not at origin station\n";
        ++fails;
      }
      if (!nearlyVec(s1.posKm, p1, 1e-3)) {
        std::cerr << "[test_traffic_lanes] arrive endpoint not at destination station\n";
        ++fails;
      }

      // Arc easing should also yield ~zero lateral velocity at endpoints.
      const double perp0 = perpSpeedKmS(s0);
      const double perp1 = perpSpeedKmS(s1);
      if (perp0 > 1e-6) {
        std::cerr << "[test_traffic_lanes] non-zero lateral speed at depart endpoint\n";
        ++fails;
      }
      if (perp1 > 1e-6) {
        std::cerr << "[test_traffic_lanes] non-zero lateral speed at arrive endpoint\n";
        ++fails;
      }
    }
  }

  if (fails == 0) std::cout << "[test_traffic_lanes] pass\n";
  return fails;
}
