#include "stellar/sim/Warehouse.h"

#include "stellar/econ/Commodity.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {

static bool nearly(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

} // namespace

int test_warehouse() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;
  using namespace stellar::econ;

  Station st{};
  st.id = 4242;
  st.name = "Warehouse Test";
  st.type = StationType::Industrial;
  st.factionId = 7;
  st.feeRate = 0.08;

  std::vector<StationStorage> storage;
  std::array<double, kCommodityCount> shipCargo{};
  shipCargo.fill(0.0);

  shipCargo[(std::size_t)CommodityId::Ore] = 100.0;

  const double rep = 0.0;

  // Deposit some ore.
  {
    const double t0 = 10.0;
    const auto r = depositToStorage(storage, st, t0, rep, shipCargo, CommodityId::Ore, 50.0);
    if (!r.ok) {
      std::cerr << "[warehouse] deposit failed\n";
      ++fails;
    }
    if (!nearly(r.unitsMoved, 50.0)) {
      std::cerr << "[warehouse] deposit moved units mismatch\n";
      ++fails;
    }

    if (!nearly(shipCargo[(std::size_t)CommodityId::Ore], 50.0)) {
      std::cerr << "[warehouse] shipCargo not reduced\n";
      ++fails;
    }

    const auto* e = findStorage(storage, st.id);
    if (!e) {
      std::cerr << "[warehouse] missing storage entry after deposit\n";
      ++fails;
    } else if (!nearly(e->cargo[(std::size_t)CommodityId::Ore], 50.0)) {
      std::cerr << "[warehouse] stored cargo mismatch\n";
      ++fails;
    }
  }

  // Fee accrual should be positive and proportional to time.
  {
    auto* e = findStorage(storage, st.id);
    if (!e) {
      std::cerr << "[warehouse] missing entry for fee test\n";
      ++fails;
    } else {
      StationStorage tmp = *e;
      const double daily = estimateStorageDailyFeeCr(tmp, rep);
      if (!(daily > 0.0)) {
        std::cerr << "[warehouse] expected positive daily fee\n";
        ++fails;
      }

      const double t1 = tmp.lastFeeDay + 10.0;
      accrueStorageFees(tmp, t1, rep);
      if (!nearly(tmp.feesDueCr, daily * 10.0, 1e-3)) {
        std::cerr << "[warehouse] fee accrual not proportional (got=" << tmp.feesDueCr << " expected=" << daily * 10.0 << ")\n";
        ++fails;
      }
    }
  }

  // Withdrawal should be blocked if fees are due and credits are insufficient.
  {
    auto* e = findStorage(storage, st.id);
    if (!e) {
      std::cerr << "[warehouse] missing entry for withdraw test\n";
      ++fails;
    } else {
      const double t2 = e->lastFeeDay + 5.0;
      accrueStorageFees(*e, t2, rep);

      double credits = 0.0;
      const auto r = withdrawFromStorage(storage, st, t2, rep, shipCargo, credits, 240.0, CommodityId::Ore, 10.0);
      if (r.ok || std::string(r.reason ? r.reason : "") != "fees_due") {
        std::cerr << "[warehouse] expected withdraw to fail with fees_due\n";
        ++fails;
      }
    }
  }

  // With enough credits, withdrawal should auto-pay fees.
  {
    auto* e = findStorage(storage, st.id);
    if (!e) {
      std::cerr << "[warehouse] missing entry for auto-pay test\n";
      ++fails;
    } else {
      const double t3 = e->lastFeeDay + 1.0;
      accrueStorageFees(*e, t3, rep);

      const double due = e->feesDueCr;
      if (!(due > 0.0)) {
        std::cerr << "[warehouse] expected fees due before auto-pay withdraw\n";
        ++fails;
      }

      double credits = due + 100.0;
      const double shipOreBefore = shipCargo[(std::size_t)CommodityId::Ore];
      const auto r = withdrawFromStorage(storage, st, t3, rep, shipCargo, credits, 240.0, CommodityId::Ore, 10.0);
      if (!r.ok) {
        std::cerr << "[warehouse] withdraw failed unexpectedly\n";
        ++fails;
      }
      if (!(r.creditsPaid > 0.0)) {
        std::cerr << "[warehouse] expected auto-pay to spend credits\n";
        ++fails;
      }
      if (!nearly(shipCargo[(std::size_t)CommodityId::Ore], shipOreBefore + r.unitsMoved)) {
        std::cerr << "[warehouse] ship cargo not increased by moved units\n";
        ++fails;
      }

      // After auto-pay, fees should be zero.
      if (e->feesDueCr > 1e-3) {
        std::cerr << "[warehouse] expected fees to be cleared by auto-pay withdraw\n";
        ++fails;
      }
    }
  }

  // Capacity enforcement: massive deposit should not exceed capacity.
  {
    shipCargo[(std::size_t)CommodityId::Ore] = 100000.0;
    const double t4 = 123.0;
    const auto r = depositToStorage(storage, st, t4, rep, shipCargo, CommodityId::Ore, 100000.0);
    if (!r.ok) {
      std::cerr << "[warehouse] large deposit failed\n";
      ++fails;
    }

    const auto* e = findStorage(storage, st.id);
    if (e) {
      const double mass = storageMassKg(*e);
      const double cap = storageCapacityKg(*e);
      if (mass > cap + 1e-3) {
        std::cerr << "[warehouse] storage exceeded capacity\n";
        ++fails;
      }
    }
  }

  // Prune should remove empty/no-fee entries.
  {
    // Empty the station storage manually.
    if (auto* e = findStorage(storage, st.id)) {
      e->cargo.fill(0.0);
      e->feesDueCr = 0.0;
    }
    pruneEmptyStorage(storage);
    if (findStorage(storage, st.id) != nullptr) {
      std::cerr << "[warehouse] pruneEmptyStorage did not remove empty entry\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_warehouse] pass\n";
  return fails;
}
