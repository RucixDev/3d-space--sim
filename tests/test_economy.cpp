#include "stellar/econ/Market.h"
#include "stellar/sim/Universe.h"

#include <iostream>

int test_economy() {
  int fails = 0;

  stellar::sim::Universe u(123);
  auto stubs = u.queryNearby({0,0,0}, 60.0, 8);
  if (stubs.empty()) {
    std::cerr << "[test_economy] no systems found\n";
    return 1;
  }

  const auto& sys = u.getSystem(stubs.front().id, &stubs.front());
  if (sys.stations.empty()) {
    std::cerr << "[test_economy] system has no stations\n";
    return 1;
  }

  const auto& station = sys.stations.front();

  // Determinism check: two fresh universes with the same seed should
  // produce identical station economy states at the same timestamp.
  {
    stellar::sim::Universe uA(123);
    stellar::sim::Universe uB(123);

    auto stubsA = uA.queryNearby({0,0,0}, 60.0, 8);
    auto stubsB = uB.queryNearby({0,0,0}, 60.0, 8);

    const auto& sysA = uA.getSystem(stubsA.front().id, &stubsA.front());
    const auto& sysB = uB.getSystem(stubsB.front().id, &stubsB.front());

    const auto& stA = sysA.stations.front();
    const auto& stB = sysB.stations.front();

    auto& econA = uA.stationEconomy(stA, 10.0);
    auto& econB = uB.stationEconomy(stB, 10.0);

    const auto cid = static_cast<std::size_t>(stellar::econ::CommodityId::Food);
    const double invA = econA.inventory[cid];
    const double invB = econB.inventory[cid];
    if (std::abs(invA - invB) > 1e-8) {
      std::cerr << "[test_economy] determinism mismatch (Food inventory) "
                << invA << " vs " << invB << "\n";
      ++fails;
    }
  }

  auto& econ0 = u.stationEconomy(station, 0.0);
  const auto q0 = stellar::econ::quote(econ0, station.economyModel, stellar::econ::CommodityId::Food);

  auto& econ10 = u.stationEconomy(station, 10.0);
  const auto q10 = stellar::econ::quote(econ10, station.economyModel, stellar::econ::CommodityId::Food);

  if (econ10.history[static_cast<std::size_t>(stellar::econ::CommodityId::Food)].empty()) {
    std::cerr << "[test_economy] expected price history samples\n";
    ++fails;
  }

  // Basic trade check
  double credits = 1000.0;
  const double invBefore = econ10.inventory[static_cast<std::size_t>(stellar::econ::CommodityId::Food)];
  auto tr = stellar::econ::buy(econ10, station.economyModel, stellar::econ::CommodityId::Food, 1.0, credits, 0.10, station.feeRate);
  if (!tr.ok) {
    std::cerr << "[test_economy] buy failed: " << (tr.reason ? tr.reason : "") << "\n";
    ++fails;
  } else {
    const double invAfter = econ10.inventory[static_cast<std::size_t>(stellar::econ::CommodityId::Food)];
    if (!(invAfter < invBefore)) {
      std::cerr << "[test_economy] inventory didn't decrease after buy\n";
      ++fails;
    }
    if (!(credits < 1000.0)) {
      std::cerr << "[test_economy] credits didn't decrease after buy\n";
      ++fails;
    }
  }

  // Inventory helper clamps
  {
    using stellar::econ::CommodityId;
    const std::size_t iFood = static_cast<std::size_t>(CommodityId::Food);
    const double before = econ10.inventory[iFood];

    const double take = stellar::econ::takeInventory(econ10, station.economyModel, CommodityId::Food, before + 9999.0);
    if (std::abs(take - before) > 1e-9) {
      std::cerr << "[test_economy] takeInventory didn't clamp to available inventory\n";
      ++fails;
    }
    if (econ10.inventory[iFood] < -1e-9) {
      std::cerr << "[test_economy] takeInventory produced negative inventory\n";
      ++fails;
    }

    // Fill to capacity and ensure addInventory clamps.
    econ10.inventory[iFood] = 0.0;
    const double cap = station.economyModel.capacity[iFood];
    const double add = stellar::econ::addInventory(econ10, station.economyModel, CommodityId::Food, cap + 9999.0);
    if (std::abs(add - cap) > 1e-9) {
      std::cerr << "[test_economy] addInventory didn't clamp to capacity\n";
      ++fails;
    }
    if (econ10.inventory[iFood] > cap + 1e-9) {
      std::cerr << "[test_economy] addInventory exceeded capacity\n";
      ++fails;
    }
  }

  // Robustness: clamp fee/spread so pathological inputs cannot produce
  // negative payouts or NaNs.
  {
    using namespace stellar::econ;

    StationEconomyModel m{};
    StationEconomyState s{};
    for (std::size_t i = 0; i < kCommodityCount; ++i) {
      m.capacity[i] = 100.0;
      m.desiredStock[i] = 50.0;
      s.inventory[i] = 10.0;
    }

    CommodityId cid = CommodityId::Food;
    const std::size_t iFood2 = static_cast<std::size_t>(cid);

    // Out-of-range spread should not produce negative bids.
    {
      const auto q = quote(s, m, cid, /*bidAskSpread=*/10.0);
      if (q.bid < -1e-9) {
        std::cerr << "[test_economy] quote() produced negative bid with huge spread\n";
        ++fails;
      }
      if (!(q.bid <= q.mid + 1e-9 && q.mid <= q.ask + 1e-9)) {
        std::cerr << "[test_economy] quote() invariant bid<=mid<=ask violated\n";
        ++fails;
      }
    }

    // Fee > 1 should be clamped so selling never *reduces* player credits.
    {
      double credits2 = 0.0;
      const double before = credits2;
      const auto tr = sell(s, m, cid, 1.0, credits2, 0.10, /*fee=*/2.0);
      if (!tr.ok) {
        std::cerr << "[test_economy] sell() failed under fee clamp\n";
        ++fails;
      } else if (credits2 + 1e-9 < before) {
        std::cerr << "[test_economy] sell() decreased credits with fee>1 (should clamp)\n";
        ++fails;
      } else if (tr.creditsDelta < -1e-9) {
        std::cerr << "[test_economy] sell() returned negative creditsDelta with fee>1\n";
        ++fails;
      }
    }

    // Corrupted station inventory (>capacity) should be handled safely by helpers.
    {
      s.inventory[iFood2] = 1e9; // corrupted
      const double taken = takeInventory(s, m, cid, 1.0);
      (void)taken;
      if (s.inventory[iFood2] > m.capacity[iFood2] + 1e-6) {
        std::cerr << "[test_economy] takeInventory() did not clamp oversized inventory\n";
        ++fails;
      }
    }
  }

  (void)q0; (void)q10;

  if (fails == 0) std::cout << "[test_economy] pass\n";
  return fails;
}
