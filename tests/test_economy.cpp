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

  (void)q0; (void)q10;

  if (fails == 0) std::cout << "[test_economy] pass\n";
  return fails;
}
