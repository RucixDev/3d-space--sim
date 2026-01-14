#include "stellar/proc/TradeRoutes.h"

#include "stellar/proc/TradeProfile.h"

#include "test_harness.h"

#include <cmath>
#include <vector>

using namespace stellar;

int test_trade_routes_hyperlanes() {
  int failures = 0;
  // 3-system toy: origin wants to export ORE.
  sim::SystemStub origin{};
  origin.id = 1;
  origin.posLy = {0.0, 0.0, 0.0};
  origin.seed = 111;

  sim::SystemStub a{};
  a.id = 2;
  a.posLy = {5.0, 0.0, 0.0};
  a.seed = 222;

  sim::SystemStub b{};
  b.id = 3;
  b.posLy = {20.0, 0.0, 0.0};
  b.seed = 333;

  // Profiles: origin exports Ore strongly. A imports Ore strongly, B imports Ore slightly less.
  proc::TradeProfile op{};
  op.exportScore.fill(0.0);
  op.importScore.fill(0.0);
  op.exportScore[(std::size_t)econ::CommodityId::Ore] = 1.0;
  op.lawlessness = 0.2;

  proc::TradeProfile ap{};
  ap.exportScore.fill(0.0);
  ap.importScore.fill(0.0);
  ap.importScore[(std::size_t)econ::CommodityId::Ore] = 1.0;
  ap.lawlessness = 0.2;

  proc::TradeProfile bp{};
  bp.exportScore.fill(0.0);
  bp.importScore.fill(0.0);
  bp.importScore[(std::size_t)econ::CommodityId::Ore] = 0.8;
  bp.lawlessness = 0.2;

  std::vector<sim::SystemStub> cands{origin, a, b};
  std::vector<proc::TradeProfile> cps{op, ap, bp};

  proc::TradeRouteParams rp{};
  rp.maxRoutes = 4;
  rp.dropWeakRoutes = false;
  rp.maxDistanceLy = 0.0;
  rp.distanceExponent = 1.0;

  // Direct (no lanes): best export should be A (closer + same affinity).
  {
    const auto direct = proc::suggestTradeRoutes(origin, op, cands, cps, rp);
    CHECK(!direct.exports.empty());
    CHECK(direct.exports[0].otherSystem == a.id);
  }

  // Hyperlane graph: origin <-> B only. A is unreachable.
  proc::HyperlaneNetwork net{};
  proc::HyperlaneEdge e{};
  e.a = 1;
  e.b = 3;
  e.distanceLy = 20.0;
  e.bandwidth01 = 1.0;
  e.risk01 = 0.1;
  net.edges.push_back(e);

  // Lane routing should drop A and pick B.
  {
    proc::HyperlaneTravelParams tp{};
    tp.riskWeight = 0.0;
    tp.bandwidthBias = 0.0;
    tp.minBandwidthFactor = 1.0;

    const auto lanes = proc::suggestTradeRoutes(origin, op, cands, cps, net, rp, tp);
    CHECK(!lanes.exports.empty());
    CHECK(lanes.exports[0].otherSystem == b.id);
    CHECK(lanes.exports[0].laneHops == 1);
    CHECK(std::abs(lanes.exports[0].laneDistanceLy - 20.0) < 1e-9);
    CHECK(std::abs(lanes.exports[0].directDistanceLy - 20.0) < 1e-9);
    CHECK(lanes.exports[0].laneBottleneckBandwidth01 > 0.9);
  }

  return failures;
}
