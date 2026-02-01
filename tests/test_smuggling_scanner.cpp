#include "test_harness.h"

#include "stellar/sim/SmugglingScanner.h"

#include "stellar/econ/Commodity.h"
#include "stellar/sim/Contraband.h"
#include "stellar/sim/Universe.h"

#include <cmath>
#include <vector>

using namespace stellar;

static bool approxEq(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

int test_smuggling_scanner() {
  int failures = 0;

  sim::Universe u(42ull);

  // Pick a deterministic origin.
  const auto stubs = u.queryNearby(math::Vec3d{0.0, 0.0, 0.0}, 350.0, 180);
  CHECK(!stubs.empty());
  if (stubs.empty()) return failures;

  const sim::SystemStub originStub = stubs.front();
  const auto& originSys = u.getSystem(originStub.id, &originStub);
  CHECK(!originSys.stations.empty());
  if (originSys.stations.empty()) return failures;

  const sim::Station originStation = originSys.stations.front();
  const core::u64 seed = u.seed();

  // Find a destination station with at least one illegal commodity that is legal at origin.
  bool found = false;
  sim::SystemStub destStub{};
  sim::Station destStation{};
  econ::CommodityId chosenCid = econ::CommodityId::Food;

  for (int pass = 0; pass < 2 && !found; ++pass) {
    for (std::size_t si = 0; si < stubs.size() && !found; ++si) {
      const auto& stub = stubs[si];
      if (stub.id == 0 || stub.id == originStub.id) continue;

      const auto& sys = u.getSystem(stub.id, &stub);
      if (sys.stations.empty()) continue;
      if (pass == 0 && sys.stations.size() != 1) continue; // prefer single-station systems for stability

      for (const auto& st : sys.stations) {
        if (st.id == 0) continue;

        const core::u32 illegalMask = sim::illegalCommodityMaskForStation(seed, st.factionId, st.id, st.type);
        if (illegalMask == 0u) continue;

        const std::size_t maxBits = std::min<std::size_t>(econ::kCommodityCount, 32);
        for (std::size_t i = 0; i < maxBits; ++i) {
          if ((illegalMask & ((core::u32)1u << (core::u32)i)) == 0u) continue;
          const auto cid = static_cast<econ::CommodityId>(i);

          // Must be legal at origin so we can buy it on the official market.
          if (sim::isIllegalCommodityAtStation(seed,
                                              originStation.factionId,
                                              originStation.id,
                                              originStation.type,
                                              cid)) {
            continue;
          }

          // Must have real capacity on both ends so the scanner can size the trade.
          if (!(originStation.economyModel.capacity[i] > 1e-6)) continue;
          if (!(st.economyModel.capacity[i] > 1e-6)) continue;

          found = true;
          destStub = stub;
          destStation = st;
          chosenCid = cid;
          break;
        }

        if (found) break;
      }
    }
  }

  CHECK(found);
  if (!found) return failures;

  const double timeDays = 0.0;
  const std::size_t idx = static_cast<std::size_t>(chosenCid);

  // Force a strongly profitable, single-commodity smuggling setup:
  //  - Origin: chosen commodity fully stocked (cheap), everything else empty (expensive).
  //  - Destination: chosen commodity empty (high demand), everything else full (no demand).
  auto& fromEcon = u.stationEconomy(originStation, timeDays);
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double cap = std::max(0.0, originStation.economyModel.capacity[i]);
    fromEcon.inventory[i] = (i == idx) ? cap : 0.0;
  }

  auto& toEcon = u.stationEconomy(destStation, timeDays);
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double cap = std::max(0.0, destStation.economyModel.capacity[i]);
    toEcon.inventory[i] = (i == idx) ? 0.0 : cap;
  }

  sim::SmugglingScanParams p;
  p.cargoCapacityKg = 250.0;
  p.cargoUsedKg = 0.0;
  p.useFreeHold = true;
  p.perStationLimit = 4;
  p.maxResults = 16;
  p.availability = sim::SmugglingAvailabilityMode::Ignore; // keep test independent of daily availability rolls
  p.scoreMode = sim::SmugglingScoreMode::RiskAdjusted;
  p.riskLambda = 0.25;
  p.minScoreCr = 1.0;

  const std::vector<sim::SystemStub> candidates = {destStub};

  const auto r1 = sim::scanSmugglingOpportunities(u, originStub, originStation, timeDays, candidates, p);
  const auto r2 = sim::scanSmugglingOpportunities(u, originStub, originStation, timeDays, candidates, p);

  CHECK(!r1.empty());
  CHECK(r1.size() == r2.size());
  if (r1.empty() || r1.size() != r2.size()) return failures;

  // Ensure the engineered opportunity appears and has sane economics.
  bool sawChosen = false;
  for (std::size_t i = 0; i < r1.size(); ++i) {
    const auto& a = r1[i];
    const auto& b = r2[i];

    CHECK(a.toSystem == b.toSystem);
    CHECK(a.toStation == b.toStation);
    CHECK(a.commodity == b.commodity);
    CHECK(approxEq(a.scoreCr, b.scoreCr));
    CHECK(approxEq(a.expectedProfitCr, b.expectedProfitCr));
    CHECK(approxEq(a.cleanProfitCr, b.cleanProfitCr));

    if (a.toSystem == destStub.id && a.toStation == destStation.id && a.commodity == chosenCid) {
      sawChosen = true;
      CHECK(a.unitsPossible > 0.0);
      CHECK(a.buyAskNetCr > 0.0);
      CHECK(a.sellBidCr > a.buyAskNetCr);
      CHECK(a.cleanProfitCr > 0.0);
      CHECK(a.scoreCr > 0.0);
    }
  }
  CHECK(sawChosen);

  return failures;
}
