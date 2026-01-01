#include "stellar/sim/TradeRunPlanner.h"

#include "stellar/econ/Commodity.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int test_trade_runs() {
  int fails = 0;

  stellar::sim::Universe u(42);
  auto stubs = u.queryNearby({0,0,0}, 200.0, 96);
  if (stubs.size() < 3) {
    std::cerr << "[test_trade_runs] expected at least three systems\n";
    return 1;
  }

  // Find three systems with at least one station each.
  int idxA = -1;
  int idxB = -1;
  int idxC = -1;

  for (int i = 0; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[i].id, &stubs[i]);
    if (!sys.stations.empty()) { idxA = i; break; }
  }
  for (int i = idxA + 1; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[i].id, &stubs[i]);
    if (!sys.stations.empty()) { idxB = i; break; }
  }
  for (int i = idxB + 1; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[i].id, &stubs[i]);
    if (!sys.stations.empty()) { idxC = i; break; }
  }

  if (idxA < 0 || idxB < 0 || idxC < 0) {
    std::cerr << "[test_trade_runs] could not find three systems with stations\n";
    return 1;
  }

  const auto& stubA = stubs[(std::size_t)idxA];
  const auto& stubB = stubs[(std::size_t)idxB];
  const auto& stubC = stubs[(std::size_t)idxC];

  const auto& sysA = u.getSystem(stubA.id, &stubA);
  const auto& sysB = u.getSystem(stubB.id, &stubB);
  const auto& sysC = u.getSystem(stubC.id, &stubC);

  const auto& stA = sysA.stations.front();
  const auto& stB = sysB.stations.front();
  const auto& stC = sysC.stations.front();

  // Force deterministic conditions for a profitable 2-leg chain:
  //   A -> B: Food
  //   B -> C: Ore
  auto& econA = u.stationEconomy(stA, 0.0);
  auto& econB = u.stationEconomy(stB, 0.0);
  auto& econC = u.stationEconomy(stC, 0.0);

  const std::size_t food = (std::size_t)stellar::econ::CommodityId::Food;
  const std::size_t ore  = (std::size_t)stellar::econ::CommodityId::Ore;

  // A has abundant Food, B has none -> A sells Food profitably to B.
  econA.inventory[food] = std::max(econA.inventory[food], stA.economyModel.capacity[food] * 0.95);
  econB.inventory[food] = 0.0;

  // B has abundant Ore, C has none -> B sells Ore profitably to C.
  econB.inventory[ore] = std::max(econB.inventory[ore], stB.economyModel.capacity[ore] * 0.95);
  econC.inventory[ore] = 0.0;

  // Avoid making B->A attractive via Ore.
  econA.inventory[ore] = std::max(econA.inventory[ore], stA.economyModel.capacity[ore] * 0.95);

  stellar::sim::TradeRunScanParams p{};
  p.manifest.cargoCapacityKg = 50.0;
  p.manifest.bidAskSpread = 0.10;
  p.manifest.stepKg = 1.0;
  p.manifest.simulatePriceImpact = false;

  p.legs = 2;
  p.beamWidth = 24;
  p.maxLegCandidates = 10;
  p.maxResults = 8;
  p.maxStations = 128;
  p.includeSameSystem = true;
  p.loopless = true;

  // Use a huge jump range so reachability is not the limiting factor in this unit test.
  p.jumpRangeLy = 10000.0;
  p.routeCostPerJump = 1.0;
  p.routeCostPerLy = 0.0;
  p.scoreMode = stellar::sim::TradeRunScoreMode::TotalProfit;

  std::vector<stellar::sim::SystemStub> candidates;
  candidates.push_back(stubA);
  candidates.push_back(stubB);
  candidates.push_back(stubC);

  const auto runs = stellar::sim::planTradeRuns(u, stubA, stA, 0.0, candidates, p);
  if (runs.empty()) {
    std::cerr << "[test_trade_runs] expected at least one run\n";
    return 1;
  }

  const auto& best = runs.front();
  if (best.legs.size() != 2) {
    std::cerr << "[test_trade_runs] expected 2-leg run, got " << best.legs.size() << "\n";
    ++fails;
  }

  if (!best.legs.empty() && best.legs[0].toStation != stB.id) {
    std::cerr << "[test_trade_runs] expected first leg toStation=" << (std::uint64_t)stB.id
              << " got " << (std::uint64_t)best.legs[0].toStation << "\n";
    ++fails;
  }

  if (best.legs.size() >= 2 && best.legs[1].toStation != stC.id) {
    std::cerr << "[test_trade_runs] expected second leg toStation=" << (std::uint64_t)stC.id
              << " got " << (std::uint64_t)best.legs[1].toStation << "\n";
    ++fails;
  }

  if (!(best.totalProfitCr > 0.0)) {
    std::cerr << "[test_trade_runs] expected positive total profit, got " << best.totalProfitCr << "\n";
    ++fails;
  }

  // Determinism: running again should match the same best run (station sequence + profit).
  const auto runs2 = stellar::sim::planTradeRuns(u, stubA, stA, 0.0, candidates, p);
  if (runs2.empty()) {
    std::cerr << "[test_trade_runs] second run returned no runs\n";
    ++fails;
  } else {
    const auto& best2 = runs2.front();
    if (best2.legs.size() != best.legs.size()) {
      std::cerr << "[test_trade_runs] determinism mismatch: legs count changed\n";
      ++fails;
    } else if (best2.legs.size() >= 2) {
      if (best2.legs[0].toStation != best.legs[0].toStation ||
          best2.legs[1].toStation != best.legs[1].toStation) {
        std::cerr << "[test_trade_runs] determinism mismatch: station sequence changed\n";
        ++fails;
      }
    }
    if (std::fabs(best2.totalProfitCr - best.totalProfitCr) > 1e-6) {
      std::cerr << "[test_trade_runs] determinism mismatch: total profit changed\n";
      ++fails;
    }
  }

  return fails;
}
