#include "stellar/sim/NavRoute.h"

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

int test_nav() {
  int fails = 0;

  using namespace stellar;

  auto makeStub = [](sim::SystemId id, double x) {
    sim::SystemStub s{};
    s.id = id;
    s.seed = id * 1337ULL;
    s.name = "S" + std::to_string((unsigned long long)id);
    s.posLy = math::Vec3d{x, 0.0, 0.0};
    s.primaryClass = sim::StarClass::G;
    s.planetCount = 1;
    s.stationCount = 1;
    s.factionId = 0;
    return s;
  };

  // A simple linear chain of systems 8 ly apart.
  std::vector<sim::SystemStub> nodes;
  nodes.push_back(makeStub(1, 0.0));
  nodes.push_back(makeStub(2, 8.0));
  nodes.push_back(makeStub(3, 16.0));
  nodes.push_back(makeStub(4, 24.0));
  nodes.push_back(makeStub(5, 32.0));
  nodes.push_back(makeStub(6, 40.0));
  nodes.push_back(makeStub(99, 1000.0)); // distractor far away

  sim::RoutePlanStats stats{};
  const auto route = sim::plotRouteAStarHops(nodes, 1, 6, 10.0, &stats);

  if (route.empty()) {
    std::cerr << "[test_nav] expected a route, got empty\n";
    ++fails;
  } else {
    if (route.front() != 1 || route.back() != 6) {
      std::cerr << "[test_nav] route endpoints mismatch\n";
      ++fails;
    }
    if (route.size() != 6) {
      std::cerr << "[test_nav] expected 6 nodes in route, got " << route.size() << "\n";
      ++fails;
    }
    std::string err;
    if (!sim::validateRoute(nodes, route, 10.0, &err)) {
      std::cerr << "[test_nav] route validation failed: " << err << "\n";
      ++fails;
    }
    if (!stats.reached || stats.hops != 5) {
      std::cerr << "[test_nav] stats mismatch reached=" << stats.reached << " hops=" << stats.hops << "\n";
      ++fails;
    }

    // Determinism check (same input => same output).
    sim::RoutePlanStats stats2{};
    const auto route2 = sim::plotRouteAStarHops(nodes, 1, 6, 10.0, &stats2);
    if (route2 != route) {
      std::cerr << "[test_nav] non-deterministic result\n";
      ++fails;
    }
  }

  // Unreachable: gap larger than jump range.
  std::vector<sim::SystemStub> gap;
  gap.push_back(makeStub(1, 0.0));
  gap.push_back(makeStub(2, 8.0));
  gap.push_back(makeStub(4, 24.0)); // gap (8 -> 24) = 16 > 10
  gap.push_back(makeStub(6, 40.0));

  const auto noRoute = sim::plotRouteAStarHops(gap, 1, 6, 10.0, nullptr);
  if (!noRoute.empty()) {
    std::cerr << "[test_nav] expected no route in gapped graph, got size=" << noRoute.size() << "\n";
    ++fails;
  }

  // start==goal: should trivially succeed (if present).
  const auto self = sim::plotRouteAStarHops(nodes, 3, 3, 10.0, nullptr);
  if (self.size() != 1 || self.front() != 3) {
    std::cerr << "[test_nav] expected self route [3], got size=" << self.size() << "\n";
    ++fails;
  }


  // Cost models: hops vs distance vs fuel.
  // Build a tiny graph where the minimum-hop path is a detour, but the minimum-distance
  // path uses more (shorter) hops.
  auto makeStubXY = [](sim::SystemId id, double x, double y) {
    sim::SystemStub s{};
    s.id = id;
    s.seed = id * 1337ULL;
    s.name = "S" + std::to_string((unsigned long long)id);
    s.posLy = math::Vec3d{x, y, 0.0};
    s.primaryClass = sim::StarClass::G;
    s.planetCount = 1;
    s.stationCount = 1;
    s.factionId = 0;
    return s;
  };

  std::vector<sim::SystemStub> costNodes;
  costNodes.push_back(makeStubXY(1, 0.0, 0.0));    // start
  costNodes.push_back(makeStubXY(2, 5.0, 0.0));    // waypoint B (near start)
  costNodes.push_back(makeStubXY(3, 21.0, 0.0));   // waypoint C (near goal, but out of start range)
  costNodes.push_back(makeStubXY(4, 30.0, 0.0));   // goal
  costNodes.push_back(makeStubXY(5, 15.0, 13.0));  // detour node A (bridges start->goal in 2 hops)

  const double jr = 20.0;

  // Hop-minimizing route should take the detour: 1 -> 5 -> 4 (2 hops).
  sim::RoutePlanStats hs{};
  const auto hopRoute = sim::plotRouteAStarHops(costNodes, 1, 4, jr, &hs);
  if (hopRoute != std::vector<sim::SystemId>({1, 5, 4})) {
    std::cerr << "[test_nav] expected hop route [1,5,4], got size=" << hopRoute.size() << "\n";
    ++fails;
  } else {
    std::string err;
    if (!sim::validateRoute(costNodes, hopRoute, jr, &err)) {
      std::cerr << "[test_nav] hop route validation failed: " << err << "\n";
      ++fails;
    }
  }

  // Distance-minimizing route should take the 3-hop near-straight path: 1 -> 2 -> 3 -> 4.
  sim::RoutePlanStats ds{};
  const auto distRoute = sim::plotRouteAStarCost(costNodes, 1, 4, jr, /*costPerJump=*/0.0, /*costPerLy=*/1.0, &ds);
  if (distRoute != std::vector<sim::SystemId>({1, 2, 3, 4})) {
    std::cerr << "[test_nav] expected distance route [1,2,3,4], got size=" << distRoute.size() << "\n";
    ++fails;
  } else {
    std::string err;
    if (!sim::validateRoute(costNodes, distRoute, jr, &err)) {
      std::cerr << "[test_nav] distance route validation failed: " << err << "\n";
      ++fails;
    }
    const double expectCost = sim::routeCost(costNodes, distRoute, 0.0, 1.0);
    if (std::abs(ds.cost - expectCost) > 1e-6) {
      std::cerr << "[test_nav] distance cost mismatch stats=" << ds.cost << " expect=" << expectCost << "\n";
      ++fails;
    }
  }

  // Fuel-like cost model should also prefer the 3-hop path in this setup.
  sim::RoutePlanStats fs{};
  const auto fuelRoute = sim::plotRouteAStarCost(costNodes, 1, 4, jr, /*costPerJump=*/2.0, /*costPerLy=*/0.5, &fs);
  if (fuelRoute != std::vector<sim::SystemId>({1, 2, 3, 4})) {
    std::cerr << "[test_nav] expected fuel route [1,2,3,4], got size=" << fuelRoute.size() << "\n";
    ++fails;
  } else {
    const double expectFuel = sim::routeCost(costNodes, fuelRoute, 2.0, 0.5);
    if (std::abs(fs.cost - expectFuel) > 1e-6) {
      std::cerr << "[test_nav] fuel cost mismatch stats=" << fs.cost << " expect=" << expectFuel << "\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_nav] pass\n";
  return fails;
}
