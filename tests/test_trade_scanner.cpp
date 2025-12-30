#include "stellar/econ/RoutePlanner.h"
#include "stellar/sim/TradeScanner.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <iostream>

int test_trade_scanner() {
  int fails = 0;

  stellar::sim::Universe u(42);
  auto stubs = u.queryNearby({0,0,0}, 120.0, 64);
  if (stubs.size() < 2) {
    std::cerr << "[test_trade_scanner] expected at least two systems\n";
    return 1;
  }

  const auto cid = stellar::econ::CommodityId::Food;
  const int cidx = (int)cid;

  // Pick two different stations with non-trivial capacity for Food.
  stellar::sim::SystemStub fromStub{};
  stellar::sim::SystemStub toStub{};
  stellar::sim::SystemId fromSysId = 0;
  stellar::sim::SystemId toSysId = 0;
  stellar::sim::StationId fromStationId = 0;
  stellar::sim::StationId toStationId = 0;

  auto pickStation = [&](bool wantDifferentSystem) {
    for (const auto& stub : stubs) {
      if (wantDifferentSystem && stub.id == fromSysId) continue;

      const auto& sys = u.getSystem(stub.id, &stub);
      for (const auto& st : sys.stations) {
        if (st.economyModel.capacity[cidx] > 50.0) {
          if (!wantDifferentSystem) {
            fromStub = stub;
            fromSysId = stub.id;
            fromStationId = st.id;
          } else {
            toStub = stub;
            toSysId = stub.id;
            toStationId = st.id;
          }
          return true;
        }
      }
    }
    return false;
  };

  if (!pickStation(false) || fromStationId == 0) {
    std::cerr << "[test_trade_scanner] couldn't find origin station\n";
    return 1;
  }
  if (!pickStation(true) || toStationId == 0) {
    std::cerr << "[test_trade_scanner] couldn't find destination station\n";
    return 1;
  }

  const auto& fromSys = u.getSystem(fromSysId, &fromStub);
  const auto& toSys = u.getSystem(toSysId, &toStub);

  const stellar::sim::Station* fromStPtr = nullptr;
  const stellar::sim::Station* toStPtr = nullptr;
  for (const auto& st : fromSys.stations) {
    if (st.id == fromStationId) {
      fromStPtr = &st;
      break;
    }
  }
  for (const auto& st : toSys.stations) {
    if (st.id == toStationId) {
      toStPtr = &st;
      break;
    }
  }
  if (!fromStPtr || !toStPtr) {
    std::cerr << "[test_trade_scanner] station lookup failed\n";
    return 1;
  }
  const auto& fromSt = *fromStPtr;
  const auto& toSt = *toStPtr;

  // Force deterministic conditions for a profitable Food route.
  auto& fromEcon = u.stationEconomy(fromSt, 0.0);
  auto& toEcon = u.stationEconomy(toSt, 0.0);

  fromEcon.inventory[cidx] = std::max(fromEcon.inventory[cidx], fromSt.economyModel.capacity[cidx] * 0.95);
  toEcon.inventory[cidx] = 0.0;

  const double feeFrom = 0.10;
  const double feeTo = 0.05;
  const double cargoKg = 50.0;

  auto feeFn = [&](const stellar::sim::Station& st) {
    return (st.id == fromStationId) ? feeFrom : feeTo;
  };

  const auto routes = stellar::econ::bestRoutesForCargo(fromEcon,
                                                        fromSt.economyModel,
                                                        toEcon,
                                                        toSt.economyModel,
                                                        cargoKg,
                                                        feeFrom,
                                                        feeTo,
                                                        0.10,
                                                        stellar::econ::kCommodityCount);

  auto expectedIt = std::find_if(routes.begin(), routes.end(), [&](const auto& r) { return r.commodity == cid; });
  if (expectedIt == routes.end()) {
    std::cerr << "[test_trade_scanner] expected a Food opportunity, but none was found\n";
    return 1;
  }
  const auto expected = *expectedIt;

  stellar::sim::TradeScanParams scan;
  scan.maxResults = 8;
  scan.perStationLimit = 1;
  scan.cargoCapacityKg = cargoKg;
  scan.cargoUsedKg = 0.0;
  scan.useFreeHold = true;
  scan.bidAskSpread = 0.10;
  scan.minNetProfit = 0.0;
  scan.includeSameSystem = false;
  scan.commodityFilterEnabled = true;
  scan.commodityFilter = cid;

  const std::vector<stellar::sim::SystemStub> candidates{fromStub, toStub};
  const auto got = stellar::sim::scanTradeOpportunities(u,
                                                        fromStub,
                                                        fromSt,
                                                        0.0,
                                                        candidates,
                                                        scan,
                                                        feeFn);

  if (got.empty()) {
    std::cerr << "[test_trade_scanner] scanner returned no results\n";
    return 1;
  }

  const auto& best = got.front();
  if (best.toSystem != toSysId) {
    std::cerr << "[test_trade_scanner] expected toSystem=" << (std::uint64_t)toSysId
              << " got " << (std::uint64_t)best.toSystem << "\n";
    ++fails;
  }
  if (best.toStation != toStationId) {
    std::cerr << "[test_trade_scanner] expected toStation=" << (std::uint64_t)toStationId
              << " got " << (std::uint64_t)best.toStation << "\n";
    ++fails;
  }
  if (best.commodity != cid) {
    std::cerr << "[test_trade_scanner] expected commodity=Food\n";
    ++fails;
  }

  const double eps = 1e-6;
  if (std::abs(best.netProfitTotal - expected.netProfitTotal) > eps) {
    std::cerr << "[test_trade_scanner] netProfitTotal mismatch: expected " << expected.netProfitTotal
              << " got " << best.netProfitTotal << "\n";
    ++fails;
  }
  if (std::abs(best.unitsPossible - expected.unitsPossible) > eps) {
    std::cerr << "[test_trade_scanner] unitsPossible mismatch: expected " << expected.unitsPossible
              << " got " << best.unitsPossible << "\n";
    ++fails;
  }

  // Min-profit threshold should filter it out.
  scan.minNetProfit = expected.netProfitTotal + 1.0;
  const auto got2 = stellar::sim::scanTradeOpportunities(u,
                                                         fromStub,
                                                         fromSt,
                                                         0.0,
                                                         candidates,
                                                         scan,
                                                         feeFn);
  if (!got2.empty()) {
    std::cerr << "[test_trade_scanner] expected no results with high minNetProfit\n";
    ++fails;
  }

  return fails;
}
