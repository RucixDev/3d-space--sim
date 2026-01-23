#include "test_harness.h"

#include "stellar/sim/BlackMarket.h"
#include "stellar/sim/SecurityModel.h"
#include "stellar/sim/Universe.h"

#include <array>
#include <cmath>

using namespace stellar;

int test_black_market() {
  int failures = 0;

  const core::u64 seed = 1337ull;
  sim::Universe u(seed);
  const auto stubs = u.queryNearby({0, 0, 0}, 60.0, 8);
  CHECK(!stubs.empty());
  if (stubs.empty()) return failures;

  const auto& sys = u.getSystem(stubs[0].id, &stubs[0]);
  CHECK(!sys.stations.empty());
  if (sys.stations.empty()) return failures;

  const auto& st = sys.stations[0];
  const auto sec = sim::systemSecurityProfile(seed, sys);
  const auto law = sim::lawProfile(seed, st.factionId);

  // Determinism: same inputs => same outputs.
  const auto bmA = sim::blackMarketProfile(seed, sys.stub.id, st, sec, law, /*timeDays=*/5.0, /*rep=*/10.0);
  const auto bmB = sim::blackMarketProfile(seed, sys.stub.id, st, sec, law, /*timeDays=*/5.0, /*rep=*/10.0);
  CHECK(bmA.available == bmB.available);
  CHECK(std::abs(bmA.access01 - bmB.access01) < 1e-12);
  CHECK(std::abs(bmA.risk01 - bmB.risk01) < 1e-12);
  CHECK(std::abs(bmA.fenceCut - bmB.fenceCut) < 1e-12);
  CHECK(std::abs(bmA.bidMul - bmB.bidMul) < 1e-12);
  CHECK(std::abs(bmA.askMul - bmB.askMul) < 1e-12);
  CHECK(std::abs(bmA.stingChance - bmB.stingChance) < 1e-12);

  // Bounds.
  CHECK(bmA.access01 >= 0.0 && bmA.access01 <= 1.0);
  CHECK(bmA.risk01 >= 0.0 && bmA.risk01 <= 1.0);
  CHECK(bmA.fenceCut >= 0.0 && bmA.fenceCut <= 1.0);
  CHECK(bmA.bidMul > 0.0);
  CHECK(bmA.askMul > 0.0);
  CHECK(bmA.stingChance >= 0.0 && bmA.stingChance <= 1.0);

  // Find an illegal commodity at this station.
  econ::CommodityId illegalCid = econ::CommodityId::Food;
  bool foundIllegal = false;
  for (std::size_t i = 0; i < std::min<std::size_t>(econ::kCommodityCount, 32); ++i) {
    const auto cid = (econ::CommodityId)i;
    if (sim::blackMarketEligibleCommodity(seed, st, cid)) {
      illegalCid = cid;
      foundIllegal = true;
      break;
    }
  }
  CHECK(foundIllegal);
  if (!foundIllegal) return failures;

  // Find a day where the fence is available (within a reasonable window).
  sim::BlackMarketProfile bmAvail{};
  bool foundAvail = false;
  for (int day = 0; day < 64; ++day) {
    const auto bm = sim::blackMarketProfile(seed, sys.stub.id, st, sec, law, (double)day, /*rep=*/25.0);
    if (bm.available) {
      bmAvail = bm;
      foundAvail = true;
      break;
    }
  }
  CHECK(foundAvail);
  if (!foundAvail) return failures;

  // Prepare cargo with illegal goods.
  std::array<double, econ::kCommodityCount> cargo{};
  cargo.fill(0.0);
  cargo[(std::size_t)illegalCid] = 12.0;
  double credits = 1000.0;

  // Success path: pick an eventSeed that doesn't sting.
  core::u64 eventSeedNoSting = 1;
  for (core::u64 s = 1; s < 256; ++s) {
    if (!sim::rollBlackMarketSting(s, bmAvail, /*heat=*/0.0)) {
      eventSeedNoSting = s;
      break;
    }
  }

  const double creditsBefore = credits;
  const double cargoBefore = cargo[(std::size_t)illegalCid];

  const auto r = sim::sellToBlackMarket(seed,
                                       eventSeedNoSting,
                                       st,
                                       bmAvail,
                                       law,
                                       /*heat=*/0.0,
                                       illegalCid,
                                       /*units=*/5.0,
                                       /*bidAskSpread=*/0.10,
                                       /*midOverride=*/nullptr,
                                       credits,
                                       cargo);
  CHECK(r.ok);
  CHECK(!r.stung);
  CHECK(r.unitsSold > 0.0);
  CHECK(r.payoutCr >= 0.0);
  CHECK(credits >= creditsBefore);
  CHECK(cargo[(std::size_t)illegalCid] <= cargoBefore);

  // Sting path: pick an eventSeed that *does* sting.
  sim::BlackMarketProfile bmRisky = bmAvail;
  bmRisky.stingChance = std::max(0.25, bmRisky.stingChance); // ensure reasonable chance.

  core::u64 eventSeedSting = 1;
  bool foundStingSeed = false;
  for (core::u64 s = 1; s < 2048; ++s) {
    if (sim::rollBlackMarketSting(s, bmRisky, /*heat=*/50.0)) {
      eventSeedSting = s;
      foundStingSeed = true;
      break;
    }
  }
  CHECK(foundStingSeed);
  if (!foundStingSeed) return failures;

  // Add more illegal cargo to make confiscation observable.
  cargo[(std::size_t)illegalCid] += 7.0;
  const double cargoBeforeSting = cargo[(std::size_t)illegalCid];
  const double creditsBeforeSting = credits;

  const auto r2 = sim::sellToBlackMarket(seed,
                                        eventSeedSting,
                                        st,
                                        bmRisky,
                                        law,
                                        /*heat=*/50.0,
                                        illegalCid,
                                        /*units=*/3.0,
                                        /*bidAskSpread=*/0.10,
                                        /*midOverride=*/nullptr,
                                        credits,
                                        cargo);
  CHECK(r2.ok);
  CHECK(r2.stung);
  CHECK(r2.scan.illegalValueCr >= 0.0);
  CHECK(r2.enforcement.fineCr >= 0.0);
  CHECK(r2.creditsDelta <= 1e-6);
  CHECK(credits <= creditsBeforeSting + 1e-6);
  // Illegal cargo should not increase after confiscation.
  CHECK(cargo[(std::size_t)illegalCid] <= cargoBeforeSting + 1e-6);



  // ---------------------------------------------------------------------------
  // Buy path (BM ask side)
  // ---------------------------------------------------------------------------

  // Successful purchase: force sting chance to 0 to make the outcome deterministic.
  {
    sim::BlackMarketProfile bmSafe = bmAvail;
    bmSafe.stingChance = 0.0;
    bmSafe.available = true;

    std::array<double, econ::kCommodityCount> buyCargo{};
    buyCargo.fill(0.0);
    double buyCredits = 8000.0;

    const double creditsBeforeBuy = buyCredits;
    const double cargoBeforeBuy = buyCargo[(std::size_t)illegalCid];

    const auto rb = sim::buyFromBlackMarket(seed,
                                           /*eventSeed=*/123,
                                           st,
                                           bmSafe,
                                           law,
                                           /*heat=*/0.0,
                                           illegalCid,
                                           /*units=*/4.0,
                                           /*bidAskSpread=*/0.10,
                                           /*midOverride=*/nullptr,
                                           buyCredits,
                                           buyCargo);

    CHECK(rb.ok);
    CHECK(!rb.stung);
    CHECK(rb.unitsBought > 0.0);
    CHECK(rb.costCr >= 0.0);
    CHECK(buyCredits <= creditsBeforeBuy + 1e-6);
    CHECK(buyCargo[(std::size_t)illegalCid] >= cargoBeforeBuy + 1e-6);
  }

  // Sting purchase: set stingChance to 1 so it always stings, then ensure enforcement triggers.
  {
    sim::BlackMarketProfile bmAlwaysSting = bmAvail;
    bmAlwaysSting.stingChance = 1.0;
    bmAlwaysSting.available = true;

    std::array<double, econ::kCommodityCount> buyCargo{};
    buyCargo.fill(0.0);
    buyCargo[(std::size_t)illegalCid] = 2.0; // already carrying something illegal
    double buyCredits = 8000.0;

    const double creditsBeforeBuy = buyCredits;
    const double cargoBeforeBuy = buyCargo[(std::size_t)illegalCid];

    const auto rb2 = sim::buyFromBlackMarket(seed,
                                            /*eventSeed=*/777,
                                            st,
                                            bmAlwaysSting,
                                            law,
                                            /*heat=*/25.0,
                                            illegalCid,
                                            /*units=*/3.0,
                                            /*bidAskSpread=*/0.10,
                                            /*midOverride=*/nullptr,
                                            buyCredits,
                                            buyCargo);

    CHECK(rb2.ok);
    CHECK(rb2.stung);
    CHECK(rb2.scan.illegalValueCr >= 0.0);
    CHECK(rb2.enforcement.fineCr >= 0.0);
    CHECK(buyCredits <= creditsBeforeBuy + 1e-6);
    // After enforcement, illegal cargo should not increase.
    CHECK(buyCargo[(std::size_t)illegalCid] <= cargoBeforeBuy + 1e-6);
  }

  (void)creditsBefore;
  (void)creditsBeforeSting;
  (void)cargoBeforeSting;

  return failures;
}
