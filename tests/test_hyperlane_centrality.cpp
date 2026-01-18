#include "test_harness.h"

#include "stellar/proc/HyperlaneCentrality.h"
#include "stellar/proc/Hyperlanes.h"
#include "stellar/sim/Celestial.h"

#include <vector>

using namespace stellar;

static sim::SystemStub mkStub(sim::SystemId id, double x, double y) {
  sim::SystemStub s{};
  s.id = id;
  s.posLy = math::Vec3d{x, y, 0.0};
  s.primaryClass = sim::StarClass::G;
  return s;
}

int test_hyperlane_centrality() {
  int failures = 0;

  // A simple chain graph:
  // 1 -- 2 -- 3 -- 4
  // The middle edge should have the highest betweenness (traffic).
  std::vector<sim::SystemStub> nodes;
  nodes.push_back(mkStub((sim::SystemId)1, 0.0, 0.0));
  nodes.push_back(mkStub((sim::SystemId)2, 1.0, 0.0));
  nodes.push_back(mkStub((sim::SystemId)3, 2.0, 0.0));
  nodes.push_back(mkStub((sim::SystemId)4, 3.0, 0.0));

  proc::HyperlaneNetwork net;
  net.edges.push_back(proc::HyperlaneEdge{(sim::SystemId)1, (sim::SystemId)2, 1.0, 1.0, 0.0});
  net.edges.push_back(proc::HyperlaneEdge{(sim::SystemId)2, (sim::SystemId)3, 1.0, 1.0, 0.0});
  net.edges.push_back(proc::HyperlaneEdge{(sim::SystemId)3, (sim::SystemId)4, 1.0, 1.0, 0.0});

  proc::HyperlaneTravelParams travel;
  travel.riskWeight = 0.0;
  travel.bandwidthBias = 0.0;

  proc::HyperlaneBetweennessParams params;
  params.travel = travel;
  params.sampleSources = 0; // exact (all sources) for small graphs

  const auto res = proc::estimateHyperlaneBetweennessCentrality(nodes, net, params);
  CHECK((int)res.edgeBetweenness.size() == (int)net.edges.size());
  CHECK(res.maxEdge >= 0.0);

  if (res.edgeBetweenness.size() == 3) {
    const double e0 = res.edgeBetweenness[0];
    const double e1 = res.edgeBetweenness[1];
    const double e2 = res.edgeBetweenness[2];
    CHECK(e1 >= e0);
    CHECK(e1 >= e2);
  }

  return failures;
}
