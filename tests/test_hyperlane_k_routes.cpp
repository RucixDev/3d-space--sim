#include "stellar/proc/HyperlaneKRoutes.h"

#include "stellar/proc/Hyperlanes.h"
#include "stellar/sim/Celestial.h"

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

using namespace stellar;

int test_hyperlane_k_routes() {
  using sim::SystemId;
  using sim::SystemStub;

  // Build a tiny deterministic graph with 3 simple routes between 1 and 4.
  // Layout (all distances = 1):
  //   1-2-4
  //   1-3-4
  //   1-2-5-6-4
  std::vector<SystemStub> nodes;
  auto mk = [](SystemId id, double x, double y) {
    SystemStub s{};
    s.id = id;
    s.name = "S" + std::to_string((unsigned long long)id);
    s.posLy = math::Vec3d{x, y, 0.0};
    return s;
  };

  nodes.push_back(mk(1, 0.0, 0.0));
  nodes.push_back(mk(2, 1.0, 0.0));
  nodes.push_back(mk(3, 0.0, 1.0));
  nodes.push_back(mk(4, 1.0, 1.0));
  nodes.push_back(mk(5, 2.0, 0.0));
  nodes.push_back(mk(6, 2.0, 1.0));

  proc::HyperlaneNetwork net{};
  auto addEdge = [&](SystemId a, SystemId b) {
    proc::HyperlaneEdge e{};
    if ((core::u64)a < (core::u64)b) {
      e.a = a;
      e.b = b;
    } else {
      e.a = b;
      e.b = a;
    }
    e.distanceLy = 1.0;
    e.bandwidth01 = 1.0;
    e.risk01 = 0.0;
    net.edges.push_back(e);
  };

  addEdge(1, 2);
  addEdge(2, 4);
  addEdge(1, 3);
  addEdge(3, 4);
  addEdge(2, 5);
  addEdge(5, 6);
  addEdge(6, 4);

  proc::HyperlaneTravelParams travel{};
  travel.riskWeight = 0.0;
  travel.bandwidthBias = 0.0;
  travel.minBandwidthFactor = 1.0;

  const auto routes = proc::plotKHyperlaneRoutesAStarCost(nodes, net, /*start=*/1, /*goal=*/4, travel, /*k=*/3, 20000);
  assert(routes.size() == 3);

  const std::vector<SystemId> r0 = {1, 2, 4};
  const std::vector<SystemId> r1 = {1, 3, 4};
  const std::vector<SystemId> r2 = {1, 2, 5, 6, 4};

  assert(routes[0].path == r0);
  assert(routes[1].path == r1);
  assert(routes[2].path == r2);

  assert(routes[0].metrics.reachable);
  assert(routes[1].metrics.reachable);
  assert(routes[2].metrics.reachable);

  assert(routes[0].metrics.hops == 2);
  assert(routes[1].metrics.hops == 2);
  assert(routes[2].metrics.hops == 4);

  // With risk/bandwidth neutralized, cost == distance.
  assert(std::abs(routes[0].metrics.costLy - 2.0) < 1e-9);
  assert(std::abs(routes[1].metrics.costLy - 2.0) < 1e-9);
  assert(std::abs(routes[2].metrics.costLy - 4.0) < 1e-9);

  // Determinism: repeated solves should return identical ordered paths.
  const auto routes2 = proc::plotKHyperlaneRoutesAStarCost(nodes, net, 1, 4, travel, 3, 20000);
  assert(routes2.size() == routes.size());
  for (std::size_t i = 0; i < routes.size(); ++i) {
    assert(routes2[i].path == routes[i].path);
    assert(std::abs(routes2[i].metrics.costLy - routes[i].metrics.costLy) < 1e-12);
  }

  return 0;
}
