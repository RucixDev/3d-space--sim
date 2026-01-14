#include "stellar/proc/TradeFlow.h"

#include "test_harness.h"

#include <cmath>
#include <vector>

using namespace stellar;

int test_trade_flow() {
  int failures = 0;
  const core::u64 seed = 1234567;

  // 4-system toy network:
  //   1 <-> 3 <-> 2 (trade endpoints)
  //   3 <-> 4 (dead-end)
  // Only system 1/2 are compatible, so all flow should route via node 3.
  sim::SystemStub s1{};
  s1.id = 1;
  s1.posLy = {0.0, 0.0, 0.0};
  s1.seed = 111;

  sim::SystemStub s2{};
  s2.id = 2;
  s2.posLy = {20.0, 0.0, 0.0};
  s2.seed = 222;

  sim::SystemStub s3{};
  s3.id = 3;
  s3.posLy = {10.0, 0.0, 0.0};
  s3.seed = 333;

  sim::SystemStub s4{};
  s4.id = 4;
  s4.posLy = {10.0, 10.0, 0.0};
  s4.seed = 444;

  std::vector<sim::SystemStub> stubs{s1, s2, s3, s4};

  proc::TradeProfile p1{};
  proc::TradeProfile p2{};
  proc::TradeProfile p3{};
  proc::TradeProfile p4{};

  // Make only (1,2) have non-zero compatibility by requiring symmetric export/import overlap.
  p1.exportScore.fill(1.0);
  p1.importScore.fill(1.0);
  p2.exportScore.fill(1.0);
  p2.importScore.fill(1.0);
  p3.exportScore.fill(0.0);
  p3.importScore.fill(0.0);
  p4.exportScore.fill(0.0);
  p4.importScore.fill(0.0);

  // Mass doesn't matter for (3,4) since compatibility is forced to 0, but keep values sensible.
  p1.hub = p1.population = p1.wealth = 1.0;
  p2.hub = p2.population = p2.wealth = 1.0;
  p3.hub = p3.population = p3.wealth = 0.3;
  p4.hub = p4.population = p4.wealth = 0.1;

  std::vector<proc::TradeProfile> profiles{p1, p2, p3, p4};

  proc::HyperlaneNetwork net{};
  net.edges.push_back(proc::HyperlaneEdge{1, 3, 10.0, 1.0, 0.0});
  net.edges.push_back(proc::HyperlaneEdge{2, 3, 10.0, 1.0, 0.0});
  net.edges.push_back(proc::HyperlaneEdge{3, 4, 10.0, 1.0, 0.0});

  proc::TradeFlowParams fp{};
  fp.commodityCompatibilityWeight = 1.0; // requires compatibility to be >0
  fp.gravityExponent = 1.0;
  fp.minCostLy = 1.0;
  fp.travelParams.riskWeight = 0.0;
  fp.travelParams.bandwidthBias = 0.0;
  fp.travelParams.minBandwidthFactor = 1.0;

  const auto flow = proc::computeTradeFlow(seed, stubs, profiles, net, fp);
  CHECK(flow.edges.size() == 2);
  CHECK(flow.nodes.size() == 4);
  CHECK(flow.totalFlow > 0.0);

  CHECK(flow.edges[0].a == 1 && flow.edges[0].b == 3);
  CHECK(flow.edges[1].a == 2 && flow.edges[1].b == 3);

  // Both used corridors should have max normalized flow.
  CHECK(flow.edges[0].flow > 0.0);
  CHECK(flow.edges[1].flow > 0.0);
  CHECK(std::abs(flow.edges[0].flow - flow.edges[1].flow) < 1e-9);
  CHECK(flow.edges[0].flow01 > 0.999);
  CHECK(flow.edges[1].flow01 > 0.999);

  // Node 3 should be the transit hub: it sees both incident flows.
  CHECK(flow.nodes[0].id == 3);
  CHECK(flow.nodes[0].traffic01 > 0.999);

  // Nodes 1 and 2 should be equal secondaries.
  CHECK(flow.nodes[1].id == 1);
  CHECK(flow.nodes[2].id == 2);
  CHECK(std::abs(flow.nodes[1].traffic01 - 0.5) < 1e-9);
  CHECK(std::abs(flow.nodes[2].traffic01 - 0.5) < 1e-9);

  // Determinism.
  const auto flow2 = proc::computeTradeFlow(seed, stubs, profiles, net, fp);
  CHECK(flow2.edges.size() == flow.edges.size());
  CHECK(flow2.nodes.size() == flow.nodes.size());
  CHECK(std::abs(flow2.totalFlow - flow.totalFlow) < 1e-12);
  CHECK(flow2.edges[0].a == flow.edges[0].a && flow2.edges[0].b == flow.edges[0].b);
  CHECK(std::abs(flow2.edges[0].flow - flow.edges[0].flow) < 1e-12);

  return failures;
}
