#include "test_harness.h"

#include "stellar/sim/NavRouteBatch.h"
#include "stellar/sim/NavRoute.h"

#include <iostream>
#include <string>
#include <vector>

int test_nav_route_batch() {
  int failures = 0;

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

  const double jr = 10.0;

  // Batch hops solve should match the A* hop route.
  const auto batchH = sim::computeNavRouteBatchHops(nodes, /*startId=*/1, jr);
  const auto routeH = sim::routeFromBatch(nodes, batchH, /*goalId=*/6);
  CHECK(routeH == std::vector<sim::SystemId>({1, 2, 3, 4, 5, 6}));

  // start==goal
  const auto self = sim::routeFromBatch(nodes, batchH, /*goalId=*/1);
  CHECK(self == std::vector<sim::SystemId>({1}));

  // Unreachable: gap larger than jump range.
  std::vector<sim::SystemStub> gap;
  gap.push_back(makeStub(1, 0.0));
  gap.push_back(makeStub(2, 8.0));
  gap.push_back(makeStub(4, 24.0)); // (8 -> 24) = 16 > 10
  gap.push_back(makeStub(6, 40.0));

  const auto batchGap = sim::computeNavRouteBatchHops(gap, /*startId=*/1, jr);
  const auto noRoute = sim::routeFromBatch(gap, batchGap, /*goalId=*/6);
  CHECK(noRoute.empty());

  // Cost models: mirror the existing A* tests to ensure path selection matches.
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

  const double jr2 = 20.0;

  // Hops: 1 -> 5 -> 4
  {
    const auto b = sim::computeNavRouteBatchHops(costNodes, 1, jr2);
    const auto r = sim::routeFromBatch(costNodes, b, 4);
    CHECK(r == std::vector<sim::SystemId>({1, 5, 4}));
  }

  // Distance: 1 -> 2 -> 3 -> 4
  {
    const auto b = sim::computeNavRouteBatchDistance(costNodes, 1, jr2);
    const auto r = sim::routeFromBatch(costNodes, b, 4);
    CHECK(r == std::vector<sim::SystemId>({1, 2, 3, 4}));
  }

  // Fuel-like: also 1 -> 2 -> 3 -> 4
  {
    const auto b = sim::computeNavRouteBatchCost(costNodes, 1, jr2, /*costPerJump=*/2.0, /*costPerLy=*/0.5);
    const auto r = sim::routeFromBatch(costNodes, b, 4);
    CHECK(r == std::vector<sim::SystemId>({1, 2, 3, 4}));
  }

  if (failures == 0) {
    std::cout << "[test_nav_route_batch] pass\n";
  }
  return failures;
}
