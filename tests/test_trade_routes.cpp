#include "test_harness.h"

#include "stellar/proc/TradeRoutes.h"

#include <array>
#include <vector>

int test_trade_routes() {
  int failures = 0;

  using namespace stellar;
  using namespace stellar::econ;
  using namespace stellar::proc;
  using namespace stellar::sim;

  SystemStub origin;
  origin.id = 1;
  origin.name = "Origin";
  origin.posLy = math::Vec3d{0.0, 0.0, 0.0};
  origin.factionId = 42;

  TradeProfile op{};
  op.hub = 0.65;
  op.population = 0.70;
  op.wealth = 0.45;
  op.lawlessness = 0.20;
  op.exportScore.fill(0.0);
  op.importScore.fill(0.0);
  op.exportScore[(std::size_t)CommodityId::Ore] = 1.0;
  op.importScore[(std::size_t)CommodityId::Medicine] = 1.0;

  SystemStub a;
  a.id = 2;
  a.name = "A";
  a.posLy = math::Vec3d{10.0, 0.0, 0.0};
  a.factionId = 42; // same faction (bonus)

  TradeProfile ap{};
  ap.hub = 0.50;
  ap.population = 0.60;
  ap.wealth = 0.40;
  ap.lawlessness = 0.15;
  ap.exportScore.fill(0.0);
  ap.importScore.fill(0.0);
  ap.importScore[(std::size_t)CommodityId::Ore] = 1.0;
  ap.exportScore[(std::size_t)CommodityId::Medicine] = 0.20;

  SystemStub b;
  b.id = 3;
  b.name = "B";
  b.posLy = math::Vec3d{12.0, 0.0, 0.0};
  b.factionId = 7;

  TradeProfile bp{};
  bp.hub = 0.50;
  bp.population = 0.60;
  bp.wealth = 0.40;
  bp.lawlessness = 0.15;
  bp.exportScore.fill(0.0);
  bp.importScore.fill(0.0);
  bp.importScore[(std::size_t)CommodityId::Ore] = 0.20;
  bp.exportScore[(std::size_t)CommodityId::Medicine] = 0.05;

  SystemStub c;
  c.id = 4;
  c.name = "C";
  c.posLy = math::Vec3d{11.0, 0.0, 0.0};
  c.factionId = 7;

  TradeProfile cp{};
  cp.hub = 0.55;
  cp.population = 0.62;
  cp.wealth = 0.42;
  cp.lawlessness = 0.15;
  cp.exportScore.fill(0.0);
  cp.importScore.fill(0.0);
  cp.exportScore[(std::size_t)CommodityId::Medicine] = 1.0;
  cp.importScore[(std::size_t)CommodityId::Ore] = 0.05;

  std::vector<SystemStub> cand{origin, a, b, c};
  std::vector<TradeProfile> prof{op, ap, bp, cp};

  TradeRouteParams params{};
  params.maxRoutes = 4;
  params.maxDistanceLy = 0.0;
  params.distanceExponent = 1.35;
  params.distanceSofteningLy = 1.0;
  params.dropWeakRoutes = true;
  params.minAffinity = 0.01;

  const auto routes = suggestTradeRoutes(origin, op, cand, prof, params);

  // Export side: ore should want to go to A (highest import affinity, plus same-faction bonus).
  CHECK(!routes.exports.empty());
  CHECK(routes.exports[0].otherSystem == 2);
  CHECK(routes.exports[0].commodity == CommodityId::Ore);

  // Import side: medicine should come from C (strongest export affinity).
  CHECK(!routes.imports.empty());
  CHECK(routes.imports[0].otherSystem == 4);
  CHECK(routes.imports[0].commodity == CommodityId::Medicine);

  // Max-routes clamp.
  TradeRouteParams p2 = params;
  p2.maxRoutes = 1;
  const auto routes2 = suggestTradeRoutes(origin, op, cand, prof, p2);
  CHECK(routes2.exports.size() <= 1);
  CHECK(routes2.imports.size() <= 1);

  return failures;
}
