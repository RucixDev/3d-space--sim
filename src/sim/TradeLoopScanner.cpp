#include "stellar/sim/TradeLoopScanner.h"

#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace stellar::sim {

static double clampFee(double feeRate) {
  if (!std::isfinite(feeRate)) return 0.0;
  // Fees in this prototype are expected to be small (0..25%), but be generous.
  return std::clamp(feeRate, 0.0, 0.95);
}

static TradeLoopFeeRateFn defaultFeeFn() {
  return [](const Station& st) { return st.feeRate; };
}

static double distLy(const math::Vec3d& a, const math::Vec3d& b) {
  const auto d = a - b;
  return std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
}

struct StopHandle {
  const SystemStub* stub{nullptr}; // for names + positions + getSystem hint
  SystemId systemId{0};
  std::size_t stationIdx{0};

  // Stable tie-break key.
  core::u64 stationId{0};
};

static const Station& getStation(Universe& u, const StopHandle& h) {
  const auto& sys = u.getSystem(h.systemId, h.stub);
  return sys.stations[h.stationIdx];
}

struct LegCache {
  bool done{false};
  econ::CargoManifestPlan plan{};
  double feeFrom{0.0};
  double feeTo{0.0};
  double distLy{0.0};
};

static TradeLoopLeg buildLeg(Universe& u,
                             const StopHandle& fromH,
                             const StopHandle& toH,
                             const LegCache& e) {
  const auto& fromSt = getStation(u, fromH);
  const auto& toSt = getStation(u, toH);

  TradeLoopLeg leg{};
  leg.fromSystem = fromH.systemId;
  leg.fromStation = fromSt.id;
  leg.toSystem = toH.systemId;
  leg.toStation = toSt.id;

  if (fromH.stub) leg.fromSystemName = fromH.stub->name;
  leg.fromStationName = fromSt.name;
  if (toH.stub) leg.toSystemName = toH.stub->name;
  leg.toStationName = toSt.name;

  leg.distanceLy = e.distLy;
  leg.feeFrom = e.feeFrom;
  leg.feeTo = e.feeTo;
  leg.manifest = e.plan;
  return leg;
}

static TradeLoopScanParams sanitize(const TradeLoopScanParams& in) {
  TradeLoopScanParams p = in;

  if (!std::isfinite(p.manifest.cargoCapacityKg)) p.manifest.cargoCapacityKg = 0.0;
  p.manifest.cargoCapacityKg = std::max(0.0, p.manifest.cargoCapacityKg);

  if (!std::isfinite(p.manifest.bidAskSpread)) p.manifest.bidAskSpread = 0.10;
  p.manifest.bidAskSpread = std::clamp(p.manifest.bidAskSpread, 0.0, 1.0);

  if (!std::isfinite(p.manifest.stepKg)) p.manifest.stepKg = 1.0;
  p.manifest.stepKg = std::clamp(p.manifest.stepKg, 0.1, 100000.0);

  if (!std::isfinite(p.manifest.maxBuyCreditsCr)) p.manifest.maxBuyCreditsCr = 0.0;

  if (!std::isfinite(p.minLegProfitCr)) p.minLegProfitCr = 0.0;
  if (!std::isfinite(p.minLoopProfitCr)) p.minLoopProfitCr = 0.0;
  p.minLegProfitCr = std::max(0.0, p.minLegProfitCr);
  p.minLoopProfitCr = std::max(0.0, p.minLoopProfitCr);

  p.legs = std::clamp<std::size_t>(p.legs, 2, 3);
  p.maxLegCandidates = std::max<std::size_t>(1, p.maxLegCandidates);
  p.maxResults = std::max<std::size_t>(1, p.maxResults);
  p.maxStations = std::max<std::size_t>(2, p.maxStations);

  return p;
}

std::vector<TradeLoop> scanTradeLoops(Universe& u,
                                      const SystemStub& originStub,
                                      const Station& originStation,
                                      double timeDays,
                                      const std::vector<SystemStub>& candidates,
                                      const TradeLoopScanParams& params,
                                      TradeLoopFeeRateFn feeRate) {
  std::vector<TradeLoop> loops;

  const auto p = sanitize(params);
  if (p.manifest.cargoCapacityKg <= 1e-9) return loops;
  if (candidates.empty()) return loops;

  if (!feeRate) feeRate = defaultFeeFn();

  // Resolve origin station index.
  std::size_t originStationIdx = 0;
  {
    const auto& sys = u.getSystem(originStub.id, &originStub);
    bool found = false;
    for (std::size_t i = 0; i < sys.stations.size(); ++i) {
      if (sys.stations[i].id == originStation.id) {
        originStationIdx = i;
        found = true;
        break;
      }
    }
    if (!found) return loops;
  }

  std::vector<StopHandle> stops;
  stops.reserve(p.maxStations);

  // Origin stop is always index 0.
  stops.push_back(StopHandle{ &originStub, originStub.id, originStationIdx, (core::u64)originStation.id });

  // Candidate station nodes.
  for (const auto& stub : candidates) {
    const bool isOriginSystem = (stub.id == originStub.id);
    if (!p.includeSameSystem && isOriginSystem) continue;

    const auto& sys = u.getSystem(stub.id, &stub);
    for (std::size_t si = 0; si < sys.stations.size(); ++si) {
      const auto& st = sys.stations[si];
      if (isOriginSystem && st.id == originStation.id) continue; // already added as origin

      if (stops.size() >= p.maxStations) break;
      stops.push_back(StopHandle{ &stub, stub.id, si, (core::u64)st.id });
    }
    if (stops.size() >= p.maxStations) break;
  }

  if (stops.size() < 2) return loops;
  const std::size_t n = stops.size();

  std::vector<LegCache> cache(n * n);

  auto leg = [&](std::size_t i, std::size_t j) -> const LegCache& {
    auto& e = cache[i * n + j];
    if (e.done) return e;
    e.done = true;

    if (i == j) return e;

    const auto& fromH = stops[i];
    const auto& toH = stops[j];

    const auto& fromSt = getStation(u, fromH);
    const auto& toSt = getStation(u, toH);

    auto& fromEcon = u.stationEconomy(fromSt, timeDays);
    auto& toEcon = u.stationEconomy(toSt, timeDays);

    econ::CargoManifestParams mp = p.manifest;
    const double feeFrom = clampFee(feeRate(fromSt));
    const double feeTo = clampFee(feeRate(toSt));
    mp.fromFeeRate = feeFrom;
    mp.toFeeRate = feeTo;

    e.plan = econ::bestManifestForCargo(fromEcon, fromSt.economyModel, toEcon, toSt.economyModel, mp);
    e.feeFrom = feeFrom;
    e.feeTo = feeTo;

    const math::Vec3d a = fromH.stub ? fromH.stub->posLy : math::Vec3d{0,0,0};
    const math::Vec3d b = toH.stub ? toH.stub->posLy : math::Vec3d{0,0,0};
    e.distLy = distLy(a, b);
    return e;
  };

  struct CandidateLeg {
    std::size_t to{0};
    double profit{0.0};
  };

  auto bestOutgoing = [&](std::size_t fromIdx) {
    std::vector<CandidateLeg> out;
    out.reserve(n);

    for (std::size_t j = 0; j < n; ++j) {
      if (j == fromIdx) continue;
      const auto& e = leg(fromIdx, j);
      const double prof = e.plan.netProfitCr;
      if (!(prof + 1e-9 >= p.minLegProfitCr)) continue;
      out.push_back(CandidateLeg{ j, prof });
    }

    std::sort(out.begin(), out.end(), [&](const CandidateLeg& a, const CandidateLeg& b) {
      // Primary: profit desc.
      if (std::fabs(a.profit - b.profit) > 1e-6) return a.profit > b.profit;

      // Tie-break: deterministic by (systemId, stationId).
      const auto& ha = stops[a.to];
      const auto& hb = stops[b.to];
      if (ha.systemId != hb.systemId) return ha.systemId < hb.systemId;
      return ha.stationId < hb.stationId;
    });

    if (out.size() > p.maxLegCandidates) out.resize(p.maxLegCandidates);
    return out;
  };

  // First expansion: origin -> B candidates.
  const auto originOutgoing = bestOutgoing(0);

  // 2-leg loops: origin -> B -> origin
  for (const auto& cB : originOutgoing) {
    const std::size_t b = cB.to;
    if (b == 0) continue;

    const auto& leg0b = leg(0, b);
    const auto& legb0 = leg(b, 0);

    if (leg0b.plan.netProfitCr <= 1e-9) continue;
    if (legb0.plan.netProfitCr <= 1e-9) continue;

    TradeLoop loop{};
    loop.legs.reserve(2);
    loop.legs.push_back(buildLeg(u, stops[0], stops[b], leg0b));
    loop.legs.push_back(buildLeg(u, stops[b], stops[0], legb0));

    loop.totalProfitCr = leg0b.plan.netProfitCr + legb0.plan.netProfitCr;
    loop.totalDistanceLy = leg0b.distLy + legb0.distLy;
    loop.profitPerLy = loop.totalProfitCr / std::max(1e-6, loop.totalDistanceLy);

    if (loop.totalProfitCr + 1e-9 >= p.minLoopProfitCr) {
      loops.push_back(std::move(loop));
    }
  }

  // 3-leg loops: origin -> B -> C -> origin
  if (p.legs >= 3) {
    for (const auto& cB : originOutgoing) {
      const std::size_t b = cB.to;
      if (b == 0) continue;

      const auto outgoingB = bestOutgoing(b);
      for (const auto& cC : outgoingB) {
        const std::size_t c = cC.to;
        if (c == 0 || c == b) continue;

        const auto& leg0b = leg(0, b);
        const auto& legbc = leg(b, c);
        const auto& legc0 = leg(c, 0);

        if (leg0b.plan.netProfitCr <= 1e-9) continue;
        if (legbc.plan.netProfitCr <= 1e-9) continue;
        if (legc0.plan.netProfitCr <= 1e-9) continue;

        TradeLoop loop{};
        loop.legs.reserve(3);
        loop.legs.push_back(buildLeg(u, stops[0], stops[b], leg0b));
        loop.legs.push_back(buildLeg(u, stops[b], stops[c], legbc));
        loop.legs.push_back(buildLeg(u, stops[c], stops[0], legc0));

        loop.totalProfitCr = leg0b.plan.netProfitCr + legbc.plan.netProfitCr + legc0.plan.netProfitCr;
        loop.totalDistanceLy = leg0b.distLy + legbc.distLy + legc0.distLy;
        loop.profitPerLy = loop.totalProfitCr / std::max(1e-6, loop.totalDistanceLy);

        if (loop.totalProfitCr + 1e-9 >= p.minLoopProfitCr) {
          loops.push_back(std::move(loop));
        }
      }
    }
  }

  // Sort loops.
  std::sort(loops.begin(), loops.end(), [](const TradeLoop& a, const TradeLoop& b) {
    if (std::fabs(a.totalProfitCr - b.totalProfitCr) > 1e-6) return a.totalProfitCr > b.totalProfitCr;
    if (std::fabs(a.totalDistanceLy - b.totalDistanceLy) > 1e-6) return a.totalDistanceLy < b.totalDistanceLy;

    // Deterministic tie-break: compare station ids sequence.
    const auto nA = a.legs.size();
    const auto nB = b.legs.size();
    const auto n = std::min(nA, nB);
    for (std::size_t i = 0; i < n; ++i) {
      if (a.legs[i].toSystem != b.legs[i].toSystem) return a.legs[i].toSystem < b.legs[i].toSystem;
      if (a.legs[i].toStation != b.legs[i].toStation) return a.legs[i].toStation < b.legs[i].toStation;
    }
    return nA < nB;
  });

  if (loops.size() > p.maxResults) loops.resize(p.maxResults);
  return loops;
}

std::vector<TradeLoop> scanTradeLoops(Universe& u,
                                      const SystemStub& originStub,
                                      const Station& originStation,
                                      double timeDays,
                                      double radiusLy,
                                      std::size_t maxSystems,
                                      const TradeLoopScanParams& params,
                                      TradeLoopFeeRateFn feeRate) {
  auto stubs = u.queryNearby(originStub.posLy, radiusLy, maxSystems);
  return scanTradeLoops(u, originStub, originStation, timeDays, stubs, params, std::move(feeRate));
}

} // namespace stellar::sim
