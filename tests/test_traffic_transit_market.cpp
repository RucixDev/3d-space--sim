#include "test_harness.h"

#include "stellar/econ/Market.h"

#include "stellar/sim/TrafficTransit.h"
#include "stellar/sim/TrafficLedger.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_map>

using namespace stellar;
using namespace stellar::sim;

namespace {

struct StationSnap {
  const Station* station{};
  econ::StationEconomyState state{};
};

static double profitPerUnitCr(const Station& src,
                             const econ::StationEconomyState& srcState,
                             const Station& dst,
                             const econ::StationEconomyState& dstState,
                             econ::CommodityId cid,
                             double bidAskSpread) {
  const auto qBuy = econ::quote(srcState, src.economyModel, cid, bidAskSpread);
  const auto qSell = econ::quote(dstState, dst.economyModel, cid, bidAskSpread);

  const double buyCost = qBuy.ask * (1.0 + src.feeRate);
  const double sellRev = qSell.bid * (1.0 - dst.feeRate);
  return sellRev - buyCost;
}

} // namespace

int test_traffic_transit_market() {
  int failures = 0;

  // Find a system with >=2 stations.
  Universe finder(1337);
  const SystemStub* stub = nullptr;
  const auto stubs = finder.queryNearby({0, 0, 0}, /*radiusLy=*/500.0, /*maxResults=*/256);
  for (const auto& s : stubs) {
    if (s.stationCount == 2) { stub = &s; break; }
  }
  if (!stub) {
    for (const auto& s : stubs) {
      if (s.stationCount >= 2) { stub = &s; break; }
    }
  }

  CHECK(stub != nullptr);

  Universe u(1337);
  const StarSystem& sys = u.getSystem(stub->id, stub);
  CHECK(sys.stations.size() >= 2);

  const double t0 = 10.5;

  // Make inventories extremely lopsided so arbitrage opportunities are obvious.
  // Even-indexed stations are full, odd-indexed stations are empty.
  for (std::size_t si = 0; si < sys.stations.size(); ++si) {
    const auto& st = sys.stations[si];
    auto& es = u.stationEconomy(st, t0);
    for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
      if (si % 2 == 0) {
        // Ensure the "source" stations are clearly surplus even if desiredStock can
        // exceed capacity for some commodity in some generated systems.
        es.inventory[i] = std::max(st.economyModel.capacity[i], st.economyModel.desiredStock[i] * 1.2);
      } else {
        es.inventory[i] = 0.0;
      }
    }
  }

  // Snapshot pre-simulation states so we can compute "would this have been profitable"
  // without needing to replay the exact per-run inventory mutations.
  std::unordered_map<StationId, StationSnap> before;
  before.reserve(sys.stations.size());
  for (const auto& st : sys.stations) {
    StationSnap snap;
    snap.station = &st;
    snap.state = u.stationEconomy(st, t0); // copy
    before.emplace(st.id, std::move(snap));
  }

  TrafficLedger ledger;
  // Force shipments to remain in-flight for a while so we can assert that inventory is not credited immediately.
  ledger.params.minDurationDays = 0.6;
  ledger.params.maxDurationDays = 0.6;
  ledger.params.speedMinKmS = 0.0;
  ledger.params.speedMaxKmS = 1.0;

  std::unordered_map<SystemId, int> stamps;

  TrafficTransitTradeParams params{};
  params.model = TrafficTradeModel::MarketArbitrage;
  params.bidAskSpread = 0.12;
  params.minProfitPerUnitCr = 0.0;
  params.distancePenaltyPerKm = 0.0;
  params.profitWeight = 1.0;
  params.flowWeight = 0.0;

  simulateNpcTradeTrafficTransitEx(u, sys, t0, stamps, params, /*kMaxBackfillDays=*/14, &ledger);

  CHECK(!ledger.shipments.empty());

  // Market-arbitrage model should not choose trades that were unprofitable under
  // the *pre*-simulation inventory snapshot.
  //
  // Note: Profit can only decrease over the day (sources are drained, sinks fill),
  // so a trade profitable at decision-time must also be profitable against the
  // initial snapshot.
  for (const auto& sh : ledger.shipments) {
    if (sh.systemId != sys.stub.id) continue;

    const auto itSrc = before.find(sh.fromStation);
    const auto itDst = before.find(sh.toStation);
    CHECK(itSrc != before.end());
    CHECK(itDst != before.end());

    const Station& src = *itSrc->second.station;
    const Station& dst = *itDst->second.station;

    const double pp = profitPerUnitCr(src, itSrc->second.state, dst, itDst->second.state, sh.commodity, params.bidAskSpread);
    CHECK(pp > 0.0);
  }

  return failures;
}
