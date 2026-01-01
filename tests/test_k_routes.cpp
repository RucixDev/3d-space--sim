#include "stellar/sim/NavRoute.h"

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

int test_k_routes() {
  int fails = 0;

  using namespace stellar;

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

  // Grid-like graph with exactly 3 loopless routes from 1 -> 4.
  //  (0,1) 3 ---- 4 (1,1) ---- 6 (2,1)
  //          \     |            |
  //           \    |            |
  //            \   |            |
  //  (0,0) 1 ---- 2 (1,0) ---- 5 (2,0)
  std::vector<sim::SystemStub> nodes;
  nodes.push_back(makeStubXY(1, 0.0, 0.0));
  nodes.push_back(makeStubXY(2, 1.0, 0.0));
  nodes.push_back(makeStubXY(3, 0.0, 1.0));
  nodes.push_back(makeStubXY(4, 1.0, 1.0));
  nodes.push_back(makeStubXY(5, 2.0, 0.0));
  nodes.push_back(makeStubXY(6, 2.0, 1.0));

  const double jr = 1.05;

  const auto routes = sim::plotKRoutesAStarHops(nodes, 1, 4, jr, 3);
  if (routes.size() != 3) {
    std::cerr << "[test_k_routes] expected 3 routes, got " << routes.size() << "\n";
    ++fails;
  }

  const std::vector<std::vector<sim::SystemId>> expect = {
    {1, 2, 4},
    {1, 3, 4},
    {1, 2, 5, 6, 4},
  };

  const std::size_t take = std::min<std::size_t>(routes.size(), expect.size());
  for (std::size_t i = 0; i < take; ++i) {
    if (routes[i].path != expect[i]) {
      std::cerr << "[test_k_routes] route[" << i << "] mismatch\n";
      ++fails;
    }

    std::string err;
    if (!sim::validateRoute(nodes, routes[i].path, jr, &err)) {
      std::cerr << "[test_k_routes] validateRoute failed for route[" << i << "]: " << err << "\n";
      ++fails;
    }

    const int hops = (int)routes[i].path.size() - 1;
    if (routes[i].hops != hops) {
      std::cerr << "[test_k_routes] hops mismatch for route[" << i << "] got=" << routes[i].hops
                << " expect=" << hops << "\n";
      ++fails;
    }

    const double dist = sim::routeDistanceLy(nodes, routes[i].path);
    if (std::abs(routes[i].distanceLy - dist) > 1e-9) {
      std::cerr << "[test_k_routes] distance mismatch for route[" << i << "] got=" << routes[i].distanceLy
                << " expect=" << dist << "\n";
      ++fails;
    }

    const double cost = sim::routeCost(nodes, routes[i].path, /*costPerJump=*/1.0, /*costPerLy=*/0.0);
    if (std::abs(routes[i].cost - cost) > 1e-9) {
      std::cerr << "[test_k_routes] cost mismatch for route[" << i << "] got=" << routes[i].cost
                << " expect=" << cost << "\n";
      ++fails;
    }
  }

  // Determinism: run again and compare.
  const auto routes2 = sim::plotKRoutesAStarHops(nodes, 1, 4, jr, 3);
  if (routes2.size() != routes.size()) {
    std::cerr << "[test_k_routes] determinism size mismatch\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < routes.size(); ++i) {
      if (routes2[i].path != routes[i].path) {
        std::cerr << "[test_k_routes] non-deterministic route ordering/content at i=" << i << "\n";
        ++fails;
        break;
      }
    }
  }

  // Non-decreasing cost.
  for (std::size_t i = 1; i < routes.size(); ++i) {
    if (routes[i].cost + 1e-12 < routes[i - 1].cost) {
      std::cerr << "[test_k_routes] cost not non-decreasing\n";
      ++fails;
      break;
    }
  }

  if (fails == 0) std::cout << "[test_k_routes] pass\n";
  return fails;
}
