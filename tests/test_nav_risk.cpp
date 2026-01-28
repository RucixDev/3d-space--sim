#include "test_harness.h"

#include "stellar/sim/NavRoute.h"
#include "stellar/sim/NavRouteBatch.h"

#include <vector>

using namespace stellar;

int test_nav_risk() {
  int failures = 0;

  // Small graph with two possible routes:
  //   1 -> 2 -> 3 (shorter but dangerous)
  //   1 -> 4 -> 3 (longer but safe)
  std::vector<sim::SystemStub> nodes;
  nodes.resize(4);

  nodes[0].id = 1;
  nodes[0].posLy = {0.0, 0.0, 0.0};
  nodes[1].id = 2;
  nodes[1].posLy = {5.0, 0.0, 0.0};
  nodes[2].id = 3;
  nodes[2].posLy = {10.0, 0.0, 0.0};
  nodes[3].id = 4;
  nodes[3].posLy = {5.0, 5.0, 0.0};

  // Risk aligned with the node order above.
  // system 2 is a pirate hotspot; 4 is quiet.
  const std::vector<double> risk = {0.1, 0.9, 0.1, 0.1};

  constexpr double kMaxJump = 8.0;

  // Without risk-awareness, shortest-distance routing wins.
  const auto routeShort = sim::plotRouteAStarCost(nodes, 1, 3, kMaxJump, 0.0, 1.0);
  CHECK(routeShort.size() == 3);
  CHECK(routeShort.front() == 1);
  CHECK(routeShort.back() == 3);
  CHECK(routeShort[1] == 2);

  // With strong risk penalty, the planner should go through system 4 instead.
  const auto routeSafe = sim::plotRouteAStarCostRisk(nodes, 1, 3, kMaxJump, 0.0, 1.0, 2.0, risk);
  CHECK(routeSafe.size() == 3);
  CHECK(routeSafe.front() == 1);
  CHECK(routeSafe.back() == 3);
  CHECK(routeSafe[1] == 4);

  // Batch solver should match the same preference.
  const auto batch = sim::computeNavRouteBatchCostRisk(nodes, 1, kMaxJump, 0.0, 1.0, 2.0, risk);
  const auto routeBatch = sim::routeFromBatch(nodes, batch, 3);
  CHECK(routeBatch.size() == 3);
  CHECK(routeBatch.front() == 1);
  CHECK(routeBatch.back() == 3);
  CHECK(routeBatch[1] == 4);

  return failures;
}
