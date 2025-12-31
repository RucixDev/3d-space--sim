#include "stellar/sim/IndustryScanner.h"
#include "stellar/sim/Universe.h"
#include "stellar/econ/Market.h"

#include <algorithm>
#include <cmath>
#include <iostream>

int test_industry_scanner() {
  int fails = 0;

  stellar::sim::Universe u(42);
  auto stubs = u.queryNearby({0,0,0}, 140.0, 64);
  if (stubs.size() < 2) {
    std::cerr << "[test_industry_scanner] expected at least two systems\n";
    return 1;
  }

  const auto ore = stellar::econ::CommodityId::Ore;
  const auto met = stellar::econ::CommodityId::Metals;
  const int oreIdx = (int)ore;
  const int metIdx = (int)met;

  stellar::sim::SystemStub fromStub{};
  stellar::sim::SystemStub toStub{};
  stellar::sim::SystemId fromSysId = 0;
  stellar::sim::SystemId toSysId = 0;
  stellar::sim::StationId fromStationId = 0;
  stellar::sim::StationId toStationId = 0;

  auto pickOrigin = [&]() {
    for (const auto& stub : stubs) {
      const auto& sys = u.getSystem(stub.id, &stub);
      for (const auto& st : sys.stations) {
        if (st.type == stellar::econ::StationType::Refinery && st.economyModel.capacity[oreIdx] > 10.0) {
          fromStub = stub;
          fromSysId = stub.id;
          fromStationId = st.id;
          return true;
        }
      }
    }
    return false;
  };

  auto pickDest = [&]() {
    for (const auto& stub : stubs) {
      if (stub.id == fromSysId) continue;
      const auto& sys = u.getSystem(stub.id, &stub);
      for (const auto& st : sys.stations) {
        if (st.economyModel.capacity[metIdx] > 50.0) {
          toStub = stub;
          toSysId = stub.id;
          toStationId = st.id;
          return true;
        }
      }
    }
    return false;
  };

  if (!pickOrigin() || fromStationId == 0) {
    std::cerr << "[test_industry_scanner] couldn't find refinery origin station\n";
    return 1;
  }
  if (!pickDest() || toStationId == 0) {
    std::cerr << "[test_industry_scanner] couldn't find destination station\n";
    return 1;
  }

  const auto& fromSys = u.getSystem(fromSysId, &fromStub);
  const auto& toSys = u.getSystem(toSysId, &toStub);

  const stellar::sim::Station* fromStPtr = nullptr;
  const stellar::sim::Station* toStPtr = nullptr;
  for (const auto& st : fromSys.stations) {
    if (st.id == fromStationId) { fromStPtr = &st; break; }
  }
  for (const auto& st : toSys.stations) {
    if (st.id == toStationId) { toStPtr = &st; break; }
  }
  if (!fromStPtr || !toStPtr) {
    std::cerr << "[test_industry_scanner] station lookup failed\n";
    return 1;
  }
  const auto& fromSt = *fromStPtr;
  const auto& toSt = *toStPtr;

  // Force deterministic conditions for a profitable ore->metals processing route.
  auto& fromEcon = u.stationEconomy(fromSt, 0.0);
  auto& toEcon = u.stationEconomy(toSt, 0.0);

  // Make ore abundant at origin => cheap. Make metals scarce at destination => expensive.
  fromEcon.inventory[oreIdx] = std::max(fromEcon.inventory[oreIdx], fromSt.economyModel.capacity[oreIdx] * 0.95);
  toEcon.inventory[metIdx] = 0.0;

  const double cargoKg = 60.0;
  const double feeFrom = 0.0;
  const double feeTo = 0.0;

  auto feeFn = [&](const stellar::sim::Station&) { return 0.0; };

  // Compute expected profit using the same math as the scanner.
  const auto* recipe = stellar::sim::findIndustryRecipe(stellar::sim::IndustryRecipeId::SmeltOre);
  if (!recipe) {
    std::cerr << "[test_industry_scanner] missing SmeltOre recipe\n";
    return 1;
  }

  const auto q1 = stellar::sim::quoteIndustryOrder(*recipe, fromSt.id, fromSt.type, 1.0, feeFrom, 0.0);
  if (q1.output != met) {
    std::cerr << "[test_industry_scanner] expected SmeltOre output MET\n";
    return 1;
  }
  if (q1.inputA != ore) {
    std::cerr << "[test_industry_scanner] expected SmeltOre input ORE\n";
    return 1;
  }

  const double outPerBatch = q1.outputUnits;
  const double outMassPerBatch = outPerBatch * stellar::econ::commodityDef(met).massKg;
  const int maxByCargo = (int)std::floor(cargoKg / std::max(1e-9, outMassPerBatch) + 1e-9);

  const auto qOre = stellar::econ::quote(fromEcon, fromSt.economyModel, ore, 0.10);
  const int maxByInputs = (int)std::floor(std::max(0.0, qOre.inventory) / std::max(1e-9, q1.inputAUnits) + 1e-9);

  const auto qMet = stellar::econ::quote(toEcon, toSt.economyModel, met, 0.10);
  const double capMet = std::max(0.0, toSt.economyModel.capacity[metIdx]);
  const double spaceMet = std::max(0.0, capMet - qMet.inventory);
  const int maxBySpace = (int)std::floor(spaceMet / std::max(1e-9, outPerBatch) + 1e-9);

  const int expBatches = std::max(0, std::min({maxByCargo, maxByInputs, maxBySpace}));
  if (expBatches <= 0) {
    std::cerr << "[test_industry_scanner] expected at least one feasible batch\n";
    return 1;
  }

  const auto q = stellar::sim::quoteIndustryOrder(*recipe, fromSt.id, fromSt.type, (double)expBatches, feeFrom, 0.0);
  const auto qOre2 = stellar::econ::quote(fromEcon, fromSt.economyModel, ore, 0.10);
  const auto qMet2 = stellar::econ::quote(toEcon, toSt.economyModel, met, 0.10);

  const double expInputCost = qOre2.ask * q.inputAUnits * (1.0 + feeFrom);
  const double expRevenue = qMet2.bid * q.outputUnits * (1.0 - feeTo);
  const double expProfit = expRevenue - expInputCost - q.serviceFeeCr;

  stellar::sim::IndustryTradeScanParams scan;
  scan.maxResults = 8;
  scan.perStationLimit = 1;
  scan.cargoCapacityKg = cargoKg;
  scan.cargoUsedKg = 0.0;
  scan.useFreeHold = true;
  scan.bidAskSpread = 0.10;
  scan.processingRep = 0.0;
  scan.minNetProfit = 0.0;
  scan.includeSameSystem = false;

  const std::vector<stellar::sim::SystemStub> candidates{fromStub, toStub};
  const auto got = stellar::sim::scanIndustryTradeOpportunities(u,
                                                                fromStub,
                                                                fromSt,
                                                                0.0,
                                                                candidates,
                                                                scan,
                                                                feeFn);

  if (got.empty()) {
    std::cerr << "[test_industry_scanner] scanner returned no results\n";
    return 1;
  }

  const auto& best = got.front();
  if (best.toSystem != toSysId) {
    std::cerr << "[test_industry_scanner] expected toSystem=" << (std::uint64_t)toSysId
              << " got " << (std::uint64_t)best.toSystem << "\n";
    ++fails;
  }
  if (best.toStation != toStationId) {
    std::cerr << "[test_industry_scanner] expected toStation=" << (std::uint64_t)toStationId
              << " got " << (std::uint64_t)best.toStation << "\n";
    ++fails;
  }
  if (best.recipe != stellar::sim::IndustryRecipeId::SmeltOre) {
    std::cerr << "[test_industry_scanner] expected recipe SmeltOre\n";
    ++fails;
  }
  if (best.output != met) {
    std::cerr << "[test_industry_scanner] expected output MET\n";
    ++fails;
  }
  if (best.netProfitCr <= 0.0) {
    std::cerr << "[test_industry_scanner] expected positive profit\n";
    ++fails;
  }

  const double eps = 1e-6;
  if (std::abs(best.netProfitCr - expProfit) > eps) {
    std::cerr << "[test_industry_scanner] netProfit mismatch: expected " << expProfit
              << " got " << best.netProfitCr << "\n";
    ++fails;
  }
  if ((int)std::llround(best.batches) != expBatches) {
    std::cerr << "[test_industry_scanner] batches mismatch: expected " << expBatches
              << " got " << best.batches << "\n";
    ++fails;
  }

  // Min-profit threshold should filter it out.
  scan.minNetProfit = best.netProfitCr + 1.0;
  const auto got2 = stellar::sim::scanIndustryTradeOpportunities(u,
                                                                 fromStub,
                                                                 fromSt,
                                                                 0.0,
                                                                 candidates,
                                                                 scan,
                                                                 feeFn);
  if (!got2.empty()) {
    std::cerr << "[test_industry_scanner] expected no results with high minNetProfit\n";
    ++fails;
  }

  return fails;
}
