#include "stellar/sim/TradeLoopScanner.h"

#include "stellar/core/JobSystem.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

static bool nearlyEqual(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

static bool samePlan(const stellar::econ::CargoManifestPlan& a,
                     const stellar::econ::CargoManifestPlan& b) {
  if (!nearlyEqual(a.cargoFilledKg, b.cargoFilledKg, 1e-6)) return false;
  if (!nearlyEqual(a.netBuyCr, b.netBuyCr, 1e-4)) return false;
  if (!nearlyEqual(a.netSellCr, b.netSellCr, 1e-4)) return false;
  if (!nearlyEqual(a.netProfitCr, b.netProfitCr, 1e-4)) return false;
  if (a.lines.size() != b.lines.size()) return false;
  for (std::size_t i = 0; i < a.lines.size(); ++i) {
    const auto& la = a.lines[i];
    const auto& lb = b.lines[i];
    if (la.commodity != lb.commodity) return false;
    if (!nearlyEqual(la.units, lb.units, 1e-6)) return false;
    if (!nearlyEqual(la.massKg, lb.massKg, 1e-6)) return false;
    if (!nearlyEqual(la.netProfitCr, lb.netProfitCr, 1e-4)) return false;
    if (!nearlyEqual(la.netProfitPerKg, lb.netProfitPerKg, 1e-6)) return false;
  }
  return true;
}

static bool sameLeg(const stellar::sim::TradeLoopLeg& a,
                    const stellar::sim::TradeLoopLeg& b) {
  if (a.fromSystem != b.fromSystem) return false;
  if (a.fromStation != b.fromStation) return false;
  if (a.toSystem != b.toSystem) return false;
  if (a.toStation != b.toStation) return false;
  if (!nearlyEqual(a.distanceLy, b.distanceLy, 1e-6)) return false;
  if (!nearlyEqual(a.feeFrom, b.feeFrom, 1e-9)) return false;
  if (!nearlyEqual(a.feeTo, b.feeTo, 1e-9)) return false;
  return samePlan(a.manifest, b.manifest);
}

static bool sameLoop(const stellar::sim::TradeLoop& a,
                     const stellar::sim::TradeLoop& b) {
  if (a.legs.size() != b.legs.size()) return false;
  if (!nearlyEqual(a.totalProfitCr, b.totalProfitCr, 1e-4)) return false;
  if (!nearlyEqual(a.totalDistanceLy, b.totalDistanceLy, 1e-6)) return false;
  if (!nearlyEqual(a.profitPerLy, b.profitPerLy, 1e-6)) return false;
  for (std::size_t i = 0; i < a.legs.size(); ++i) {
    if (!sameLeg(a.legs[i], b.legs[i])) return false;
  }
  return true;
}

int test_trade_loops_parallel() {
  int fails = 0;

  stellar::sim::Universe u(42);
  auto stubs = u.queryNearby({0, 0, 0}, 120.0, 64);
  if (stubs.size() < 2) {
    std::cerr << "[test_trade_loops_parallel] expected at least two systems\n";
    return 1;
  }

  // Find two systems with at least one station each.
  int idxA = -1;
  int idxB = -1;
  for (int i = 0; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[(std::size_t)i].id, &stubs[(std::size_t)i]);
    if (!sys.stations.empty()) {
      idxA = i;
      break;
    }
  }
  for (int i = idxA + 1; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[(std::size_t)i].id, &stubs[(std::size_t)i]);
    if (!sys.stations.empty()) {
      idxB = i;
      break;
    }
  }
  if (idxA < 0 || idxB < 0) {
    std::cerr << "[test_trade_loops_parallel] could not find two systems with stations\n";
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
  const std::size_t ore = (std::size_t)stellar::econ::CommodityId::Ore;

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

  const auto serial = stellar::sim::scanTradeLoops(u, fromStub, fromSt, 0.0, candidates, p);
  if (serial.empty()) {
    std::cerr << "[test_trade_loops_parallel] expected at least one serial loop\n";
    return 1;
  }

  stellar::core::JobSystem jobs(4);
  const auto parallel = stellar::sim::scanTradeLoopsParallel(jobs, u, fromStub, fromSt, 0.0, candidates, p);
  if (parallel.empty()) {
    std::cerr << "[test_trade_loops_parallel] expected at least one parallel loop\n";
    return 1;
  }

  // Basic shape checks (same as the serial test).
  const auto& best = parallel.front();
  if (best.legs.size() != 2) {
    std::cerr << "[test_trade_loops_parallel] expected 2-leg loop, got " << best.legs.size() << "\n";
    ++fails;
  }

  // First leg should go to the selected destination station.
  if (!best.legs.empty() && best.legs[0].toStation != toSt.id) {
    std::cerr << "[test_trade_loops_parallel] expected first leg toStation=" << (std::uint64_t)toSt.id
              << " got " << (std::uint64_t)best.legs[0].toStation << "\n";
    ++fails;
  }

  // Second leg should return to origin.
  if (best.legs.size() >= 2 && best.legs[1].toStation != fromSt.id) {
    std::cerr << "[test_trade_loops_parallel] expected second leg toStation=" << (std::uint64_t)fromSt.id
              << " got " << (std::uint64_t)best.legs[1].toStation << "\n";
    ++fails;
  }

  if (!(best.totalProfitCr > 0.0)) {
    std::cerr << "[test_trade_loops_parallel] expected positive total profit, got " << best.totalProfitCr << "\n";
    ++fails;
  }

  // Strong equivalence: serial and parallel should match exactly after sorting.
  if (serial.size() != parallel.size()) {
    std::cerr << "[test_trade_loops_parallel] size mismatch: serial=" << serial.size()
              << " parallel=" << parallel.size() << "\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < serial.size(); ++i) {
      if (!sameLoop(serial[i], parallel[i])) {
        std::cerr << "[test_trade_loops_parallel] loop mismatch at i=" << i << "\n";
        ++fails;
        break;
      }
    }
  }

  // Determinism: running again should match exactly.
  const auto parallel2 = stellar::sim::scanTradeLoopsParallel(jobs, u, fromStub, fromSt, 0.0, candidates, p);
  if (parallel2.size() != parallel.size()) {
    std::cerr << "[test_trade_loops_parallel] determinism mismatch: size changed\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < parallel.size(); ++i) {
      if (!sameLoop(parallel[i], parallel2[i])) {
        std::cerr << "[test_trade_loops_parallel] determinism mismatch at i=" << i << "\n";
        ++fails;
        break;
      }
    }
  }

  return fails;
}
