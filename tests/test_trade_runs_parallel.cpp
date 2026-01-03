#include "stellar/sim/TradeRunPlanner.h"

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

static bool sameLeg(const stellar::sim::TradeRunLeg& a,
                    const stellar::sim::TradeRunLeg& b) {
  if (a.fromSystem != b.fromSystem) return false;
  if (a.fromStation != b.fromStation) return false;
  if (a.toSystem != b.toSystem) return false;
  if (a.toStation != b.toStation) return false;

  if (a.route != b.route) return false;
  if (a.routeHops != b.routeHops) return false;
  if (!nearlyEqual(a.routeDistanceLy, b.routeDistanceLy, 1e-6)) return false;
  if (!nearlyEqual(a.routeCost, b.routeCost, 1e-6)) return false;

  if (!nearlyEqual(a.feeFrom, b.feeFrom, 1e-9)) return false;
  if (!nearlyEqual(a.feeTo, b.feeTo, 1e-9)) return false;

  return samePlan(a.manifest, b.manifest);
}

static bool sameRun(const stellar::sim::TradeRun& a,
                    const stellar::sim::TradeRun& b) {
  if (a.legs.size() != b.legs.size()) return false;

  if (!nearlyEqual(a.totalProfitCr, b.totalProfitCr, 1e-4)) return false;
  if (!nearlyEqual(a.totalRouteDistanceLy, b.totalRouteDistanceLy, 1e-6)) return false;
  if (!nearlyEqual(a.totalRouteCost, b.totalRouteCost, 1e-6)) return false;
  if (a.totalHops != b.totalHops) return false;

  if (!nearlyEqual(a.profitPerLy, b.profitPerLy, 1e-6)) return false;
  if (!nearlyEqual(a.profitPerHop, b.profitPerHop, 1e-6)) return false;
  if (!nearlyEqual(a.profitPerCost, b.profitPerCost, 1e-6)) return false;

  for (std::size_t i = 0; i < a.legs.size(); ++i) {
    if (!sameLeg(a.legs[i], b.legs[i])) return false;
  }
  return true;
}

int test_trade_runs_parallel() {
  int fails = 0;

  stellar::sim::Universe u(42);
  auto stubs = u.queryNearby({0,0,0}, 200.0, 96);
  if (stubs.size() < 3) {
    std::cerr << "[test_trade_runs_parallel] expected at least three systems\n";
    return 1;
  }

  // Find three systems with at least one station each.
  int idxA = -1;
  int idxB = -1;
  int idxC = -1;

  for (int i = 0; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[(std::size_t)i].id, &stubs[(std::size_t)i]);
    if (!sys.stations.empty()) { idxA = i; break; }
  }
  for (int i = idxA + 1; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[(std::size_t)i].id, &stubs[(std::size_t)i]);
    if (!sys.stations.empty()) { idxB = i; break; }
  }
  for (int i = idxB + 1; i < (int)stubs.size(); ++i) {
    const auto& sys = u.getSystem(stubs[(std::size_t)i].id, &stubs[(std::size_t)i]);
    if (!sys.stations.empty()) { idxC = i; break; }
  }

  if (idxA < 0 || idxB < 0 || idxC < 0) {
    std::cerr << "[test_trade_runs_parallel] could not find three systems with stations\n";
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

  econA.inventory[food] = std::max(econA.inventory[food], stA.economyModel.capacity[food] * 0.95);
  econB.inventory[food] = 0.0;

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

  const auto serial = stellar::sim::planTradeRuns(u, stubA, stA, 0.0, candidates, p);
  if (serial.empty()) {
    std::cerr << "[test_trade_runs_parallel] expected at least one serial run\n";
    return 1;
  }

  stellar::core::JobSystem jobs(4);
  const auto parallel = stellar::sim::planTradeRunsParallel(jobs, u, stubA, stA, 0.0, candidates, p);
  if (parallel.empty()) {
    std::cerr << "[test_trade_runs_parallel] expected at least one parallel run\n";
    return 1;
  }

  // Strong equivalence: serial and parallel should match exactly after sorting.
  if (serial.size() != parallel.size()) {
    std::cerr << "[test_trade_runs_parallel] size mismatch: serial=" << serial.size()
              << " parallel=" << parallel.size() << "\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < serial.size(); ++i) {
      if (!sameRun(serial[i], parallel[i])) {
        std::cerr << "[test_trade_runs_parallel] run mismatch at i=" << i << "\n";
        ++fails;
        break;
      }
    }
  }

  // Determinism: running again should match exactly.
  const auto parallel2 = stellar::sim::planTradeRunsParallel(jobs, u, stubA, stA, 0.0, candidates, p);
  if (parallel2.size() != parallel.size()) {
    std::cerr << "[test_trade_runs_parallel] determinism mismatch: size changed\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < parallel.size(); ++i) {
      if (!sameRun(parallel[i], parallel2[i])) {
        std::cerr << "[test_trade_runs_parallel] determinism mismatch at i=" << i << "\n";
        ++fails;
        break;
      }
    }
  }

  return fails;
}
