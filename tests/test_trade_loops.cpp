#include "stellar/sim/TradeLoopScanner.h"

#include "stellar/econ/Commodity.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int test_trade_loops() {
  int fails = 0;

  stellar::sim::Universe u(42);
  auto stubs = u.queryNearby({0,0,0}, 120.0, 64);
  if (stubs.size() < 2) {
    std::cerr << "[test_trade_loops] expected at least two systems\n";
    return 1;
  }

  // Find two systems with at least one station each.
  int idxA = -1;
  int idxB = -1;
  for (int i = 0; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[i].id, &stubs[i]);
    if (!sys.stations.empty()) { idxA = i; break; }
  }
  for (int i = idxA + 1; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[i].id, &stubs[i]);
    if (!sys.stations.empty()) { idxB = i; break; }
  }

  if (idxA < 0 || idxB < 0) {
    std::cerr << "[test_trade_loops] could not find two systems with stations\n";
    return 1;
  }

  const auto& fromStub = stubs[(std::size_t)idxA];
  const auto& toStub = stubs[(std::size_t)idxB];

  const auto& fromSys = u.getSystem(fromStub.id, &fromStub);
  const auto& toSys = u.getSystem(toStub.id, &toStub);

  const auto& fromSt = fromSys.stations.front();
  const auto& toSt = toSys.stations.front();

  // Force deterministic conditions for a profitable round-trip:
  //  - Food:  from -> to
  //  - Ore:   to -> from
  auto& econFrom = u.stationEconomy(fromSt, 0.0);
  auto& econTo = u.stationEconomy(toSt, 0.0);

  const std::size_t food = (std::size_t)stellar::econ::CommodityId::Food;
  const std::size_t ore  = (std::size_t)stellar::econ::CommodityId::Ore;

  econFrom.inventory[food] = std::max(econFrom.inventory[food], fromSt.economyModel.capacity[food] * 0.95);
  econTo.inventory[food] = 0.0;

  econTo.inventory[ore] = std::max(econTo.inventory[ore], toSt.economyModel.capacity[ore] * 0.95);
  econFrom.inventory[ore] = 0.0;

  stellar::sim::TradeLoopScanParams p{};
  p.manifest.cargoCapacityKg = 50.0;
  p.manifest.bidAskSpread = 0.10;
  p.manifest.stepKg = 1.0;
  p.manifest.simulatePriceImpact = false; // stable, easier to reason about

  p.legs = 2;
  p.maxLegCandidates = 8;
  p.maxResults = 8;
  p.includeSameSystem = true;

  std::vector<stellar::sim::SystemStub> candidates;
  candidates.push_back(fromStub);
  candidates.push_back(toStub);

  const auto loops = stellar::sim::scanTradeLoops(u, fromStub, fromSt, 0.0, candidates, p);
  if (loops.empty()) {
    std::cerr << "[test_trade_loops] expected at least one loop\n";
    return 1;
  }

  const auto& best = loops.front();
  if (best.legs.size() != 2) {
    std::cerr << "[test_trade_loops] expected 2-leg loop, got " << best.legs.size() << "\n";
    ++fails;
  }

  // First leg should go to the selected destination station.
  if (!best.legs.empty() && best.legs[0].toStation != toSt.id) {
    std::cerr << "[test_trade_loops] expected first leg toStation=" << (std::uint64_t)toSt.id
              << " got " << (std::uint64_t)best.legs[0].toStation << "\n";
    ++fails;
  }

  // Second leg should return to origin.
  if (best.legs.size() >= 2 && best.legs[1].toStation != fromSt.id) {
    std::cerr << "[test_trade_loops] expected second leg toStation=" << (std::uint64_t)fromSt.id
              << " got " << (std::uint64_t)best.legs[1].toStation << "\n";
    ++fails;
  }

  if (!(best.totalProfitCr > 0.0)) {
    std::cerr << "[test_trade_loops] expected positive total profit, got " << best.totalProfitCr << "\n";
    ++fails;
  }

  // Determinism: running again should match the same best leg target + profit.
  const auto loops2 = stellar::sim::scanTradeLoops(u, fromStub, fromSt, 0.0, candidates, p);
  if (loops2.empty()) {
    std::cerr << "[test_trade_loops] second run returned no loops\n";
    ++fails;
  } else {
    const auto& best2 = loops2.front();
    if (!best2.legs.empty() && best2.legs[0].toStation != best.legs[0].toStation) {
      std::cerr << "[test_trade_loops] determinism mismatch: first leg toStation changed\n";
      ++fails;
    }
    if (std::fabs(best2.totalProfitCr - best.totalProfitCr) > 1e-6) {
      std::cerr << "[test_trade_loops] determinism mismatch: total profit changed\n";
      ++fails;
    }
  }

  return fails;
}
