#include "stellar/proc/TradeIntel.h"
#include "stellar/proc/TradePrices.h"

#include "stellar/sim/Celestial.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

static constexpr std::size_t idx(econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

int test_trade_prices() {
  int failures = 0;

  // ---- Multiplier directionality -----------------------------------------
  proc::TradeProfile exporter{};
  exporter.exportScore.fill(0.0);
  exporter.importScore.fill(0.0);
  exporter.exportScore[idx(econ::CommodityId::Food)] = 1.0;
  exporter.hub = 0.0;
  exporter.population = 0.0;

  proc::TradeProfile importer{};
  importer.exportScore.fill(0.0);
  importer.importScore.fill(0.0);
  importer.importScore[idx(econ::CommodityId::Food)] = 1.0;
  importer.hub = 0.0;
  importer.population = 0.0;

  const core::u64 seed = 123456789ull;
  const double mExp = proc::estimateTradePriceMultiplier(seed, exporter, econ::CommodityId::Food);
  const double mImp = proc::estimateTradePriceMultiplier(seed, importer, econ::CommodityId::Food);

  CHECK(mExp < 1.0);
  CHECK(mImp > 1.0);

  // Determinism: repeated calls must match bit-for-bit (within FP noise bounds).
  const double mImp2 = proc::estimateTradePriceMultiplier(seed, importer, econ::CommodityId::Food);
  CHECK(std::abs(mImp - mImp2) < 1e-12);

  // Hub/pop dampening: big markets should be closer to neutral (1.0).
  proc::TradeProfile importerStable = importer;
  importerStable.hub = 1.0;
  importerStable.population = 1.0;
  const double mImpStable = proc::estimateTradePriceMultiplier(seed, importerStable, econ::CommodityId::Food);
  CHECK(mImpStable > 1.0);
  CHECK(mImpStable < mImp);

  // ---- Basic intel loop profitability ------------------------------------
  sim::SystemStub a{};
  a.id = 1;
  a.seed = 111;
  a.name = "A";
  a.posLy = {0.0, 0.0, 0.0};
  a.factionId = 1;

  sim::SystemStub b{};
  b.id = 2;
  b.seed = 222;
  b.name = "B";
  b.posLy = {12.0, 0.0, 0.0};
  b.factionId = 1;

  proc::TradeProfile pA{};
  pA.exportScore.fill(0.0);
  pA.importScore.fill(0.0);
  pA.exportScore[idx(econ::CommodityId::Food)] = 1.0;
  pA.importScore[idx(econ::CommodityId::Luxury)] = 1.0;
  pA.hub = 0.25;
  pA.population = 0.25;

  proc::TradeProfile pB{};
  pB.exportScore.fill(0.0);
  pB.importScore.fill(0.0);
  pB.importScore[idx(econ::CommodityId::Food)] = 1.0;
  pB.exportScore[idx(econ::CommodityId::Luxury)] = 1.0;
  pB.hub = 0.25;
  pB.population = 0.25;

  std::vector<sim::SystemStub> cands{b};
  std::vector<proc::TradeProfile> cps{pB};

  proc::TradeRouteParams rp{};
  rp.maxRoutes = 8;

  proc::TradeIntelParams ip{};
  ip.bidAskSpread = 0.10;
  ip.buyFeeRate = 0.0;
  ip.sellFeeRate = 0.0;
  ip.cargoCapacityKg = 200.0;
  ip.maxLoops = 8;

  const auto rep = proc::buildTradeIntel(a, pA, cands, cps, rp, ip);

  CHECK(!rep.exports.empty());
  CHECK(!rep.imports.empty());
  CHECK(!rep.loops.empty());

  // We engineered a simple complementary pair:
  //  A exports Food to B, B exports Luxury to A.
  CHECK(rep.exports[0].commodity == econ::CommodityId::Food);
  CHECK(rep.imports[0].commodity == econ::CommodityId::Luxury);

  // Both legs should be plausibly profitable in the macro model.
  CHECK(rep.exports[0].profitPerUnitCr > 0.0);
  CHECK(rep.imports[0].profitPerUnitCr > 0.0);
  CHECK(rep.loops[0].roundTripProfitCr > 0.0);

  return failures;
}
