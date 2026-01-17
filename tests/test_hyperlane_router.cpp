#include "test_harness.h"

#include "stellar/proc/HyperlaneRouter.h"

#include <cmath>
#include <string>
#include <vector>

using namespace stellar;

namespace {

sim::SystemStub mkStub(sim::SystemId id, double xLy, double yLy) {
  sim::SystemStub s{};
  s.id = id;
  s.name = "S" + std::to_string((unsigned long long)id);
  s.posLy = math::Vec3d{xLy, yLy, 0.0};
  s.primaryClass = sim::StarClass::G;
  return s;
}

bool approx(double a, double b, double eps = 1e-6) {
  return std::abs(a - b) <= eps;
}

} // namespace

int test_hyperlane_router() {
  int failures = 0;

  {
    // Graph:
    //   A(1) --10--> B(2) --10--> C(3)
    //    \                         /
    //     \--15 (low bw, risky) --/
    //
    // Expected: A->B->C wins under default travel params.

    std::vector<sim::SystemStub> stubs;
    stubs.push_back(mkStub(1, 0.0, 0.0));
    stubs.push_back(mkStub(2, 10.0, 0.0));
    stubs.push_back(mkStub(3, 20.0, 0.0));

    proc::HyperlaneNetwork net;
    net.edges = {
        proc::HyperlaneEdge{1, 2, 10.0, 1.0, 0.0},
        proc::HyperlaneEdge{2, 3, 10.0, 1.0, 0.0},
        proc::HyperlaneEdge{1, 3, 15.0, 0.2, 0.6},
    };

    proc::HyperlaneRouter r(net, stubs);
    proc::HyperlaneTravelParams tp{};
    tp.riskWeight = 0.60;
    tp.bandwidthBias = 0.65;
    tp.minBandwidthFactor = 0.35;

    CHECK(r.compute(1, tp));

    const auto m = r.metricsTo(3);
    CHECK(m.reachable);
    CHECK(m.hops == 2);
    CHECK(approx(m.distanceLy, 20.0));
    CHECK(approx(m.risk01, 0.0));
    CHECK(approx(m.bottleneckBandwidth01, 1.0));

    std::vector<sim::SystemId> path;
    CHECK(r.buildPathTo(3, path));
    CHECK(path.size() == 3);
    CHECK(path[0] == 1);
    CHECK(path[1] == 2);
    CHECK(path[2] == 3);
  }

  {
    // A(1) -- B(2) -- C(3) -- D(4)
    // with a slightly risky / lower-bw edge on C-D.

    std::vector<sim::SystemStub> stubs;
    stubs.push_back(mkStub(1, 0.0, 0.0));
    stubs.push_back(mkStub(2, 10.0, 0.0));
    stubs.push_back(mkStub(3, 20.0, 0.0));
    stubs.push_back(mkStub(4, 30.0, 0.0));

    proc::HyperlaneNetwork net;
    net.edges = {
        proc::HyperlaneEdge{1, 2, 10.0, 1.0, 0.0},
        proc::HyperlaneEdge{2, 3, 10.0, 1.0, 0.0},
        proc::HyperlaneEdge{3, 4, 10.0, 0.5, 0.1},
    };

    proc::HyperlaneRouter r(net, stubs);
    proc::HyperlaneTravelParams tp{};
    tp.riskWeight = 0.60;
    tp.bandwidthBias = 0.65;
    tp.minBandwidthFactor = 0.35;

    CHECK(r.compute(1, tp));

    const auto m = r.metricsTo(4);
    CHECK(m.reachable);
    CHECK(m.hops == 3);
    CHECK(approx(m.distanceLy, 30.0));

    // Compound risk: 1 - Π(1-risk). Only last edge has risk=0.1.
    CHECK(approx(m.risk01, 0.1));

    // Bottleneck should be the min bandwidth along the path.
    CHECK(approx(m.bottleneckBandwidth01, 0.5));
  }

  {
    // Unreachable destination => metrics.reachable=false.
    std::vector<sim::SystemStub> stubs;
    stubs.push_back(mkStub(1, 0.0, 0.0));
    stubs.push_back(mkStub(2, 10.0, 0.0));
    stubs.push_back(mkStub(3, 20.0, 0.0));

    proc::HyperlaneNetwork net;
    net.edges = {
        proc::HyperlaneEdge{1, 2, 10.0, 1.0, 0.0},
    };

    proc::HyperlaneRouter r(net, stubs);
    CHECK(r.compute(1));

    const auto m = r.metricsTo(3);
    CHECK(!m.reachable);
    CHECK(m.hops == 0);
    CHECK(approx(m.costLy, 0.0));
  }

  return failures;
}
