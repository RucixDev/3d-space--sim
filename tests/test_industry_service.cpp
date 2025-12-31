#include "stellar/sim/IndustryService.h"

#include "stellar/econ/Commodity.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {

static bool nearly(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

} // namespace

int test_industry_service() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;
  using namespace stellar::econ;

  Station st{};
  st.id = 9001;
  st.name = "Refinery Services";
  st.type = StationType::Refinery;
  st.factionId = 7;
  st.feeRate = 0.08;

  const auto* recipe = findIndustryRecipe(IndustryRecipeId::SmeltOre);
  if (!recipe) {
    std::cerr << "[industry_service] missing SmeltOre recipe\n";
    return 1;
  }

  std::vector<IndustryOrder> orders;
  core::u64 nextId = 1;

  std::array<double, kCommodityCount> cargo{};
  cargo.fill(0.0);
  cargo[(std::size_t)CommodityId::Ore] = 500.0;

  double credits = 100000.0;
  const double rep = 10.0;
  const double feeEff = st.feeRate;

  // Submit an order.
  IndustrySubmitResult s1 = submitIndustryOrder(nextId, orders, cargo, credits,
                                                st, /*timeDays=*/10.0, rep,
                                                *recipe, /*batches=*/10, feeEff);
  if (!s1.ok) {
    std::cerr << "[industry_service] submit #1 failed (" << (s1.reason ? s1.reason : "?") << ")\n";
    ++fails;
  }
  if (orders.size() != 1) {
    std::cerr << "[industry_service] expected 1 order after submit\n";
    ++fails;
  }
  if (nextId != 2) {
    std::cerr << "[industry_service] nextId did not increment\n";
    ++fails;
  }
  if (!(s1.order.readyDay > s1.order.submittedDay)) {
    std::cerr << "[industry_service] expected readyDay > submittedDay\n";
    ++fails;
  }
  if (!nearly(cargo[(std::size_t)CommodityId::Ore], 500.0 - s1.quote.inputAUnits, 1e-4)) {
    std::cerr << "[industry_service] cargo input deduction mismatch\n";
    ++fails;
  }

  // Queueing: a second submit at the same station/time should start after the first finishes.
  cargo[(std::size_t)CommodityId::Ore] += 1000.0;
  IndustrySubmitResult s2 = submitIndustryOrder(nextId, orders, cargo, credits,
                                                st, /*timeDays=*/10.0, rep,
                                                *recipe, /*batches=*/10, feeEff);
  if (!s2.ok) {
    std::cerr << "[industry_service] submit #2 failed (" << (s2.reason ? s2.reason : "?") << ")\n";
    ++fails;
  }
  if (orders.size() != 2) {
    std::cerr << "[industry_service] expected 2 orders after second submit\n";
    ++fails;
  }
  if (!nearly(s2.order.readyDay - s1.order.readyDay, s2.quote.timeDays, 1e-6)) {
    std::cerr << "[industry_service] expected queued readyDay offset by quote.timeDays\n";
    ++fails;
  }

  // Claim gating: cannot claim before ready.
  {
    std::array<double, kCommodityCount> shipCargo{};
    shipCargo.fill(0.0);
    const double capKg = 3.0; // Metals are 3kg/unit -> can hold ~1 unit.

    const auto r = claimIndustryOrderToCargo(orders[0], st, /*timeDays=*/10.0, shipCargo, capKg);
    if (r.ok || std::string(r.reason ? r.reason : "") != "not_ready") {
      std::cerr << "[industry_service] expected claim to fail with not_ready\n";
      ++fails;
    }
  }

  // Partial claim into cargo (limited by mass), then store remainder to warehouse.
  {
    std::array<double, kCommodityCount> shipCargo{};
    shipCargo.fill(0.0);
    std::vector<StationStorage> storage;

    const double capKg = 3.0; // 1 unit
    const double tReady = orders[0].readyDay + 1e-3;

    const auto r1 = claimIndustryOrderToCargo(orders[0], st, tReady, shipCargo, capKg);
    if (!r1.ok) {
      std::cerr << "[industry_service] claim to cargo failed (" << (r1.reason ? r1.reason : "?") << ")\n";
      ++fails;
    }
    if (!(r1.unitsMoved > 0.0)) {
      std::cerr << "[industry_service] expected some units moved to cargo\n";
      ++fails;
    }
    if (!(shipCargo[(std::size_t)orders[0].output] >= r1.unitsMoved - 1e-6)) {
      std::cerr << "[industry_service] shipCargo not increased by claimed units\n";
      ++fails;
    }
    if (orders[0].claimed) {
      // With a tiny cargo cap, the order should not be fully claimed yet.
      std::cerr << "[industry_service] expected partial claim to leave order unclaimed\n";
      ++fails;
    }

    const double remainingBeforeStore = orders[0].outputUnits;
    const auto r2 = moveIndustryOrderOutputToWarehouse(orders[0], st, tReady, rep, storage);
    if (!r2.ok) {
      std::cerr << "[industry_service] store to warehouse failed (" << (r2.reason ? r2.reason : "?") << ")\n";
      ++fails;
    }
    if (!orders[0].claimed) {
      std::cerr << "[industry_service] expected order to be claimed after warehouse move\n";
      ++fails;
    }
    if (!nearly(r2.unitsMoved, remainingBeforeStore, 1e-6)) {
      std::cerr << "[industry_service] warehouse moved units mismatch\n";
      ++fails;
    }

    const auto* e = findStorage(storage, st.id);
    if (!e) {
      std::cerr << "[industry_service] expected storage entry after warehouse move\n";
      ++fails;
    } else if (!(e->cargo[(std::size_t)orders[0].output] >= r2.unitsMoved - 1e-6)) {
      std::cerr << "[industry_service] warehouse cargo not increased\n";
      ++fails;
    }
  }

  // Prune should remove claimed orders.
  pruneClaimedIndustryOrders(orders);
  if (orders.size() != 1) {
    std::cerr << "[industry_service] expected prune to remove claimed orders\n";
    ++fails;
  }

  if (fails == 0) {
    std::cout << "[test_industry_service] PASS\n";
  }
  return fails;
}
