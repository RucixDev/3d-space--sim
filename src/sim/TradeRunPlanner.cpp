#include "stellar/sim/TradeRunPlanner.h"

#include "stellar/core/JobSystem.h"
#include "stellar/sim/NavRoute.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>

namespace stellar::sim {

static double clampFee(double feeRate) {
  if (!std::isfinite(feeRate)) return 0.0;
  // Fees in this prototype are expected to be small (0..25%), but be generous.
  return std::clamp(feeRate, 0.0, 0.95);
}

static TradeRunFeeRateFn defaultFeeFn() {
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

struct RouteCache {
  bool done{false};
  bool reachable{false};
  std::vector<SystemId> route{};
  int hops{0};
  double distLy{0.0};
  double cost{0.0};
};

static TradeRunScanParams sanitize(const TradeRunScanParams& in) {
  TradeRunScanParams p = in;

  if (!std::isfinite(p.manifest.cargoCapacityKg)) p.manifest.cargoCapacityKg = 0.0;
  p.manifest.cargoCapacityKg = std::max(0.0, p.manifest.cargoCapacityKg);

  if (!std::isfinite(p.manifest.bidAskSpread)) p.manifest.bidAskSpread = 0.10;
  p.manifest.bidAskSpread = std::clamp(p.manifest.bidAskSpread, 0.0, 1.0);

  if (!std::isfinite(p.manifest.stepKg)) p.manifest.stepKg = 1.0;
  p.manifest.stepKg = std::clamp(p.manifest.stepKg, 0.1, 100000.0);

  if (!std::isfinite(p.manifest.maxBuyCreditsCr)) p.manifest.maxBuyCreditsCr = 0.0;

  if (!std::isfinite(p.minLegProfitCr)) p.minLegProfitCr = 0.0;
  if (!std::isfinite(p.minRunProfitCr)) p.minRunProfitCr = 0.0;
  p.minLegProfitCr = std::max(0.0, p.minLegProfitCr);
  p.minRunProfitCr = std::max(0.0, p.minRunProfitCr);

  p.legs = std::max<std::size_t>(1, p.legs);
  p.beamWidth = std::max<std::size_t>(1, p.beamWidth);
  p.maxLegCandidates = std::max<std::size_t>(1, p.maxLegCandidates);
  p.maxResults = std::max<std::size_t>(1, p.maxResults);
  p.maxStations = std::max<std::size_t>(2, p.maxStations);

  if (!std::isfinite(p.jumpRangeLy)) p.jumpRangeLy = 0.0;
  p.jumpRangeLy = std::max(0.0, p.jumpRangeLy);

  if (!std::isfinite(p.routeCostPerJump)) p.routeCostPerJump = 1.0;
  if (!std::isfinite(p.routeCostPerLy)) p.routeCostPerLy = 0.0;
  p.routeCostPerJump = std::max(0.0, p.routeCostPerJump);
  p.routeCostPerLy = std::max(0.0, p.routeCostPerLy);

  return p;
}

static double scoreTotals(double profitCr, double distLy, int hops, double cost, TradeRunScoreMode mode) {
  profitCr = std::max(0.0, profitCr);
  distLy = std::max(0.0, distLy);
  cost = std::max(0.0, cost);
  hops = std::max(0, hops);

  switch (mode) {
    case TradeRunScoreMode::TotalProfit:
      return profitCr;
    case TradeRunScoreMode::ProfitPerLy:
      return profitCr / std::max(1e-6, distLy);
    case TradeRunScoreMode::ProfitPerHop:
      return profitCr / (double)std::max(1, hops);
    case TradeRunScoreMode::ProfitPerCost:
      return profitCr / std::max(1e-6, cost);
    default:
      return profitCr;
  }
}

static TradeRunLeg buildLeg(Universe& u,
                            const StopHandle& fromH,
                            const StopHandle& toH,
                            const LegCache& e,
                            const RouteCache& r) {
  const auto& fromSt = getStation(u, fromH);
  const auto& toSt = getStation(u, toH);

  TradeRunLeg leg{};
  leg.fromSystem = fromH.systemId;
  leg.fromStation = fromSt.id;
  leg.toSystem = toH.systemId;
  leg.toStation = toSt.id;

  if (fromH.stub) leg.fromSystemName = fromH.stub->name;
  leg.fromStationName = fromSt.name;
  if (toH.stub) leg.toSystemName = toH.stub->name;
  leg.toStationName = toSt.name;

  leg.feeFrom = e.feeFrom;
  leg.feeTo = e.feeTo;
  leg.manifest = e.plan;

  leg.route = r.route;
  leg.routeHops = r.hops;
  leg.routeDistanceLy = r.distLy;
  leg.routeCost = r.cost;
  return leg;
}

std::vector<TradeRun> planTradeRuns(Universe& u,
                                    const SystemStub& originStub,
                                    const Station& originStation,
                                    double timeDays,
                                    const std::vector<SystemStub>& candidates,
                                    const TradeRunScanParams& params,
                                    TradeRunFeeRateFn feeRate) {
  std::vector<TradeRun> runs;

  const auto p = sanitize(params);
  if (p.manifest.cargoCapacityKg <= 1e-9) return runs;
  if (candidates.empty()) return runs;

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
    if (!found) return runs;
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

  if (stops.size() < 2) return runs;
  const std::size_t nStops = stops.size();

  // Build a unique system node list for routing.
  std::vector<SystemStub> sysNodes;
  sysNodes.reserve(candidates.size() + 1);
  sysNodes.push_back(originStub);
  for (const auto& s : candidates) sysNodes.push_back(s);

  std::sort(sysNodes.begin(), sysNodes.end(), [](const SystemStub& a, const SystemStub& b) {
    return a.id < b.id;
  });
  sysNodes.erase(std::unique(sysNodes.begin(), sysNodes.end(), [](const SystemStub& a, const SystemStub& b) {
    return a.id == b.id;
  }), sysNodes.end());

  const std::size_t nSys = sysNodes.size();
  std::unordered_map<SystemId, std::size_t> sysIndex;
  sysIndex.reserve(nSys * 2);
  for (std::size_t i = 0; i < nSys; ++i) sysIndex[sysNodes[i].id] = i;

  std::vector<RouteCache> routeCache(nSys * nSys);
  auto route = [&](SystemId a, SystemId b) -> const RouteCache& {
    const auto iaIt = sysIndex.find(a);
    const auto ibIt = sysIndex.find(b);
    if (iaIt == sysIndex.end() || ibIt == sysIndex.end()) {
      static RouteCache empty{};
      return empty;
    }

    const std::size_t ia = iaIt->second;
    const std::size_t ib = ibIt->second;
    auto& e = routeCache[ia * nSys + ib];
    if (e.done) return e;
    e.done = true;

    if (a == b) {
      e.reachable = true;
      e.route = {a};
      e.hops = 0;
      e.distLy = 0.0;
      e.cost = 0.0;
      return e;
    }

    if (p.jumpRangeLy <= 1e-9) {
      // No reachability checks: treat as a direct jump.
      e.reachable = true;
      e.route = {a, b};
      e.hops = 1;
      e.distLy = distLy(sysNodes[ia].posLy, sysNodes[ib].posLy);
      e.cost = p.routeCostPerJump + p.routeCostPerLy * e.distLy;
      return e;
    }

    // Fast path: direct jump is reachable.
    {
      const double direct = distLy(sysNodes[ia].posLy, sysNodes[ib].posLy);
      if (direct <= p.jumpRangeLy + 1e-9) {
        e.reachable = true;
        e.route = {a, b};
        e.hops = 1;
        e.distLy = direct;
        e.cost = p.routeCostPerJump + p.routeCostPerLy * e.distLy;
        return e;
      }
    }

    auto r = plotRouteAStarHops(sysNodes, a, b, p.jumpRangeLy);
    if (r.empty()) {
      e.reachable = false;
      return e;
    }

    e.reachable = true;
    e.route = std::move(r);
    e.hops = (int)e.route.size() - 1;
    e.distLy = routeDistanceLy(sysNodes, e.route);
    e.cost = routeCost(sysNodes, e.route, p.routeCostPerJump, p.routeCostPerLy);
    return e;
  };


  std::vector<LegCache> legCache(nStops * nStops);
  auto leg = [&](std::size_t i, std::size_t j) -> const LegCache& {
    auto& e = legCache[i * nStops + j];
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

  std::vector<std::vector<CandidateLeg>> outgoingCache(nStops);
  std::vector<bool> outgoingDone(nStops, false);

  auto bestOutgoing = [&](std::size_t fromIdx) -> const std::vector<CandidateLeg>& {
    if (outgoingDone[fromIdx]) return outgoingCache[fromIdx];
    outgoingDone[fromIdx] = true;

    auto& out = outgoingCache[fromIdx];
    out.reserve(nStops);

    for (std::size_t j = 0; j < nStops; ++j) {
      if (j == fromIdx) continue;
      const auto& e = leg(fromIdx, j);
      const double prof = e.plan.netProfitCr;
      if (!(prof + 1e-9 >= p.minLegProfitCr)) continue;
      if (prof <= 1e-9) continue;
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

  auto pathLess = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i) {
      const auto sa = stops[a[i]].stationId;
      const auto sb = stops[b[i]].stationId;
      if (sa != sb) return sa < sb;
    }
    return a.size() < b.size();
  };

  struct Partial {
    std::vector<std::size_t> path; // sequence of stop indices (start with 0)
    double profitCr{0.0};
    double routeDistLy{0.0};
    double routeCost{0.0};
    int hops{0};
  };

  auto partialScore = [&](const Partial& s) {
    return scoreTotals(s.profitCr, s.routeDistLy, s.hops, s.routeCost, p.scoreMode);
  };


  std::vector<Partial> beam;
  beam.reserve(p.beamWidth);
  beam.push_back(Partial{{0}, 0.0, 0.0, 0.0, 0});

  for (std::size_t depth = 0; depth < p.legs; ++depth) {
    std::vector<Partial> next;

    for (const auto& st : beam) {
      const std::size_t fromIdx = st.path.back();
      const auto& out = bestOutgoing(fromIdx);
      if (out.empty()) continue;

      for (const auto& cand : out) {
        const std::size_t toIdx = cand.to;
        if (toIdx == fromIdx) continue;

        if (p.loopless) {
          bool seen = false;
          for (const auto ix : st.path) {
            if (ix == toIdx) { seen = true; break; }
          }
          if (seen) continue;
        }

        const auto& e = leg(fromIdx, toIdx);
        if (e.plan.netProfitCr <= 1e-9) continue;

        const auto& r = route(stops[fromIdx].systemId, stops[toIdx].systemId);
        if (p.jumpRangeLy > 1e-9 && !r.reachable) continue;

        Partial s = st;
        s.path.push_back(toIdx);
        s.profitCr += e.plan.netProfitCr;
        s.routeDistLy += r.distLy;
        s.routeCost += r.cost;
        s.hops += r.hops;
        next.push_back(std::move(s));
      }
    }

    if (next.empty()) break;

    std::sort(next.begin(), next.end(), [&](const Partial& a, const Partial& b) {
      const double sa = partialScore(a);
      const double sb = partialScore(b);
      if (std::fabs(sa - sb) > 1e-6) return sa > sb;

      if (std::fabs(a.profitCr - b.profitCr) > 1e-6) return a.profitCr > b.profitCr;
      return pathLess(a.path, b.path);
    });

    if (next.size() > p.beamWidth) next.resize(p.beamWidth);
    beam = std::move(next);
  }

  // Convert beam states into full TradeRun results.
  for (const auto& st : beam) {
    if (st.path.size() < 2) continue;
    if (st.path.size() != p.legs + 1) continue;
    if (!(st.profitCr + 1e-9 >= p.minRunProfitCr)) continue;

    TradeRun run{};
    run.legs.reserve(p.legs);

    for (std::size_t i = 0; i + 1 < st.path.size(); ++i) {
      const std::size_t a = st.path[i];
      const std::size_t b = st.path[i + 1];

      const auto& e = leg(a, b);
      const auto& r = route(stops[a].systemId, stops[b].systemId);
      if (p.jumpRangeLy > 1e-9 && !r.reachable) {
        run.legs.clear();
        break;
      }

      run.legs.push_back(buildLeg(u, stops[a], stops[b], e, r));
    }

    if (run.legs.empty()) continue;

    run.totalProfitCr = st.profitCr;
    run.totalRouteDistanceLy = st.routeDistLy;
    run.totalRouteCost = st.routeCost;
    run.totalHops = st.hops;

    run.profitPerLy = run.totalProfitCr / std::max(1e-6, run.totalRouteDistanceLy);
    run.profitPerHop = run.totalProfitCr / (double)std::max(1, run.totalHops);
    run.profitPerCost = run.totalProfitCr / std::max(1e-6, run.totalRouteCost);

    runs.push_back(std::move(run));
  }

  if (runs.empty()) return runs;

  std::sort(runs.begin(), runs.end(), [&](const TradeRun& a, const TradeRun& b) {
    const double sa = scoreTotals(a.totalProfitCr, a.totalRouteDistanceLy, a.totalHops, a.totalRouteCost, p.scoreMode);
    const double sb = scoreTotals(b.totalProfitCr, b.totalRouteDistanceLy, b.totalHops, b.totalRouteCost, p.scoreMode);
    if (std::fabs(sa - sb) > 1e-6) return sa > sb;
    if (std::fabs(a.totalProfitCr - b.totalProfitCr) > 1e-6) return a.totalProfitCr > b.totalProfitCr;

    // Deterministic tiebreak: lexicographic station sequence.
    const std::size_t na = a.legs.size();
    const std::size_t nb = b.legs.size();
    const std::size_t n = std::min(na, nb);
    for (std::size_t i = 0; i < n; ++i) {
      if (a.legs[i].toStation != b.legs[i].toStation) return a.legs[i].toStation < b.legs[i].toStation;
    }
    return na < nb;
  });

  if (runs.size() > p.maxResults) runs.resize(p.maxResults);
  return runs;
}

std::vector<TradeRun> planTradeRuns(Universe& u,
                                    const SystemStub& originStub,
                                    const Station& originStation,
                                    double timeDays,
                                    double radiusLy,
                                    std::size_t maxSystems,
                                    const TradeRunScanParams& params,
                                    TradeRunFeeRateFn feeRate) {
  const auto candidates = u.queryNearby(originStub.posLy, radiusLy, maxSystems);
  return planTradeRuns(u, originStub, originStation, timeDays, candidates, params, std::move(feeRate));
}

// -----------------------------------------------------------------------------
// Parallel planner
// -----------------------------------------------------------------------------

namespace {

struct StopSnapshot {
  SystemId systemId{0};
  core::u64 stationKey{0};

  // Copied station data (so we don't hold pointers into Universe caches).
  Station station{};

  std::string systemName;
  math::Vec3d systemPosLy{0,0,0};
};

static TradeRunLeg buildLeg(const StopSnapshot& from,
                            const StopSnapshot& to,
                            const LegCache& e,
                            const RouteCache& r) {
  TradeRunLeg leg{};
  leg.fromSystem = from.systemId;
  leg.fromStation = from.station.id;
  leg.toSystem = to.systemId;
  leg.toStation = to.station.id;

  leg.fromSystemName = from.systemName;
  leg.fromStationName = from.station.name;
  leg.toSystemName = to.systemName;
  leg.toStationName = to.station.name;

  leg.feeFrom = e.feeFrom;
  leg.feeTo = e.feeTo;
  leg.manifest = e.plan;

  leg.route = r.route;
  leg.routeHops = r.hops;
  leg.routeDistanceLy = r.distLy;
  leg.routeCost = r.cost;
  return leg;
}

} // namespace

std::vector<TradeRun> planTradeRunsParallel(core::JobSystem& jobs,
                                            Universe& u,
                                            const SystemStub& originStub,
                                            const Station& originStation,
                                            double timeDays,
                                            const std::vector<SystemStub>& candidates,
                                            const TradeRunScanParams& params,
                                            TradeRunFeeRateFn feeRate) {
  std::vector<TradeRun> runs;

  const auto p = sanitize(params);
  if (p.manifest.cargoCapacityKg <= 1e-9) return runs;
  if (candidates.empty()) return runs;

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
    if (!found) return runs;
  }

  // Snapshot station nodes on the calling thread (Universe caches are not thread-safe).
  std::vector<StopSnapshot> stops;
  stops.reserve(p.maxStations);

  {
    const auto& sys = u.getSystem(originStub.id, &originStub);
    StopSnapshot origin{};
    origin.systemId = originStub.id;
    origin.stationKey = (core::u64)originStation.id;
    origin.station = sys.stations[originStationIdx];
    origin.systemName = originStub.name;
    origin.systemPosLy = originStub.posLy;
    stops.push_back(std::move(origin));
  }

  for (const auto& stub : candidates) {
    const bool isOriginSystem = (stub.id == originStub.id);
    if (!p.includeSameSystem && isOriginSystem) continue;

    const auto& sys = u.getSystem(stub.id, &stub);
    for (const auto& st : sys.stations) {
      if (isOriginSystem && st.id == originStation.id) continue;
      if (stops.size() >= p.maxStations) break;

      StopSnapshot s{};
      s.systemId = stub.id;
      s.stationKey = (core::u64)st.id;
      s.station = st;
      s.systemName = stub.name;
      s.systemPosLy = stub.posLy;
      stops.push_back(std::move(s));
    }
    if (stops.size() >= p.maxStations) break;
  }

  if (stops.size() < 2) return runs;
  const std::size_t nStops = stops.size();

  // Build a unique system node list for routing.
  std::vector<SystemStub> sysNodes;
  sysNodes.reserve(candidates.size() + 1);
  sysNodes.push_back(originStub);
  for (const auto& s : candidates) sysNodes.push_back(s);

  std::sort(sysNodes.begin(), sysNodes.end(), [](const SystemStub& a, const SystemStub& b) {
    return a.id < b.id;
  });
  sysNodes.erase(std::unique(sysNodes.begin(), sysNodes.end(), [](const SystemStub& a, const SystemStub& b) {
    return a.id == b.id;
  }), sysNodes.end());

  const std::size_t nSys = sysNodes.size();
  std::unordered_map<SystemId, std::size_t> sysIndex;
  sysIndex.reserve(nSys * 2);
  for (std::size_t i = 0; i < nSys; ++i) sysIndex[sysNodes[i].id] = i;

  // Snapshot station economies and effective fees sequentially.
  std::vector<double> fee(nStops, 0.0);
  std::vector<econ::StationEconomyState> econStates(nStops);

  for (std::size_t i = 0; i < nStops; ++i) {
    fee[i] = clampFee(feeRate(stops[i].station));
    econStates[i] = u.stationEconomy(stops[i].station, timeDays);
  }

  std::vector<RouteCache> routeCache(nSys * nSys);
  auto route = [&](SystemId a, SystemId b) -> const RouteCache& {
    const auto iaIt = sysIndex.find(a);
    const auto ibIt = sysIndex.find(b);
    if (iaIt == sysIndex.end() || ibIt == sysIndex.end()) {
      static RouteCache empty{};
      return empty;
    }

    const std::size_t ia = iaIt->second;
    const std::size_t ib = ibIt->second;
    auto& e = routeCache[ia * nSys + ib];
    if (e.done) return e;
    e.done = true;

    if (a == b) {
      e.reachable = true;
      e.route = {a};
      e.hops = 0;
      e.distLy = 0.0;
      e.cost = 0.0;
      return e;
    }

    const double direct = distLy(sysNodes[ia].posLy, sysNodes[ib].posLy);

    if (p.jumpRangeLy <= 1e-9) {
      // No reachability checks: treat as a direct jump.
      e.reachable = true;
      e.route = {a, b};
      e.hops = 1;
      e.distLy = direct;
      e.cost = p.routeCostPerJump + p.routeCostPerLy * e.distLy;
      return e;
    }

    // Fast path: direct jump is reachable.
    if (direct <= p.jumpRangeLy + 1e-9) {
      e.reachable = true;
      e.route = {a, b};
      e.hops = 1;
      e.distLy = direct;
      e.cost = p.routeCostPerJump + p.routeCostPerLy * e.distLy;
      return e;
    }

    auto r = plotRouteAStarHops(sysNodes, a, b, p.jumpRangeLy);
    if (r.empty()) {
      e.reachable = false;
      return e;
    }

    e.reachable = true;
    e.route = std::move(r);
    e.hops = (int)e.route.size() - 1;
    e.distLy = routeDistanceLy(sysNodes, e.route);
    e.cost = routeCost(sysNodes, e.route, p.routeCostPerJump, p.routeCostPerLy);
    return e;
  };


  std::vector<LegCache> legCache(nStops * nStops);
  std::vector<core::u8> rowDone(nStops, 0);

  auto computeLeg = [&](std::size_t i, std::size_t j) {
    auto& e = legCache[i * nStops + j];
    e.done = true;
    if (i == j) return;

    econ::CargoManifestParams mp = p.manifest;
    mp.fromFeeRate = fee[i];
    mp.toFeeRate = fee[j];

    e.plan = econ::bestManifestForCargo(econStates[i],
                                       stops[i].station.economyModel,
                                       econStates[j],
                                       stops[j].station.economyModel,
                                       mp);
    e.feeFrom = fee[i];
    e.feeTo = fee[j];
    e.distLy = distLy(stops[i].systemPosLy, stops[j].systemPosLy);
  };

  auto computeRow = [&](std::size_t i) {
    if (rowDone[i]) return;
    rowDone[i] = 1;
    jobs.parallelFor(nStops, [&](std::size_t j) {
      if (j == i) return;
      computeLeg(i, j);
    });
  };

  struct CandidateLeg {
    std::size_t to{0};
    double profit{0.0};
  };

  const auto betterCandidate = [&](const CandidateLeg& a, const CandidateLeg& b) {
    if (std::fabs(a.profit - b.profit) > 1e-6) return a.profit > b.profit;
    const auto& ha = stops[a.to];
    const auto& hb = stops[b.to];
    if (ha.systemId != hb.systemId) return ha.systemId < hb.systemId;
    return ha.stationKey < hb.stationKey;
  };

  std::vector<std::vector<CandidateLeg>> outgoingCache(nStops);
  std::vector<core::u8> outgoingDone(nStops, 0);

  auto bestOutgoing = [&](std::size_t fromIdx) -> const std::vector<CandidateLeg>& {
    if (outgoingDone[fromIdx]) return outgoingCache[fromIdx];
    outgoingDone[fromIdx] = 1;

    computeRow(fromIdx);

    auto& out = outgoingCache[fromIdx];
    out.clear();
    out.reserve(nStops);

    for (std::size_t j = 0; j < nStops; ++j) {
      if (j == fromIdx) continue;
      const auto& e = legCache[fromIdx * nStops + j];
      const double prof = e.plan.netProfitCr;
      if (!(prof + 1e-9 >= p.minLegProfitCr)) continue;
      if (prof <= 1e-9) continue;
      out.push_back(CandidateLeg{j, prof});
    }

    std::sort(out.begin(), out.end(), [&](const CandidateLeg& a, const CandidateLeg& b) {
      return betterCandidate(a, b);
    });

    if (out.size() > p.maxLegCandidates) out.resize(p.maxLegCandidates);
    return out;
  };

  auto pathLess = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i) {
      const auto sa = stops[a[i]].stationKey;
      const auto sb = stops[b[i]].stationKey;
      if (sa != sb) return sa < sb;
    }
    return a.size() < b.size();
  };

  struct Partial {
    std::vector<std::size_t> path; // sequence of stop indices (start with 0)
    double profitCr{0.0};
    double routeDistLy{0.0};
    double routeCost{0.0};
    int hops{0};
  };

  auto partialScore = [&](const Partial& s) {
    return scoreTotals(s.profitCr, s.routeDistLy, s.hops, s.routeCost, p.scoreMode);
  };


  std::vector<Partial> beam;
  beam.reserve(p.beamWidth);
  beam.push_back(Partial{{0}, 0.0, 0.0, 0.0, 0});

  for (std::size_t depth = 0; depth < p.legs; ++depth) {
    std::vector<Partial> next;

    for (const auto& st : beam) {
      const std::size_t fromIdx = st.path.back();
      const auto& out = bestOutgoing(fromIdx);
      if (out.empty()) continue;

      for (const auto& cand : out) {
        const std::size_t toIdx = cand.to;
        if (toIdx == fromIdx) continue;

        if (p.loopless) {
          bool seen = false;
          for (const auto ix : st.path) {
            if (ix == toIdx) { seen = true; break; }
          }
          if (seen) continue;
        }

        const auto& e = legCache[fromIdx * nStops + toIdx];
        if (e.plan.netProfitCr <= 1e-9) continue;

        const auto& r = route(stops[fromIdx].systemId, stops[toIdx].systemId);
        if (p.jumpRangeLy > 1e-9 && !r.reachable) continue;

        Partial s = st;
        s.path.push_back(toIdx);
        s.profitCr += e.plan.netProfitCr;
        s.routeDistLy += r.distLy;
        s.routeCost += r.cost;
        s.hops += r.hops;
        next.push_back(std::move(s));
      }
    }

    if (next.empty()) break;

    std::sort(next.begin(), next.end(), [&](const Partial& a, const Partial& b) {
      const double sa = partialScore(a);
      const double sb = partialScore(b);
      if (std::fabs(sa - sb) > 1e-6) return sa > sb;

      if (std::fabs(a.profitCr - b.profitCr) > 1e-6) return a.profitCr > b.profitCr;
      return pathLess(a.path, b.path);
    });

    if (next.size() > p.beamWidth) next.resize(p.beamWidth);
    beam = std::move(next);
  }

  // Convert beam states into full TradeRun results.
  for (const auto& st : beam) {
    if (st.path.size() < 2) continue;
    if (st.path.size() != p.legs + 1) continue;
    if (!(st.profitCr + 1e-9 >= p.minRunProfitCr)) continue;

    TradeRun run{};
    run.legs.reserve(p.legs);

    for (std::size_t i = 0; i + 1 < st.path.size(); ++i) {
      const std::size_t a = st.path[i];
      const std::size_t b = st.path[i + 1];

      // Ensure this row is computed (so the cached manifest exists).
      computeRow(a);

      const auto& e = legCache[a * nStops + b];
      const auto& r = route(stops[a].systemId, stops[b].systemId);
      if (p.jumpRangeLy > 1e-9 && !r.reachable) {
        run.legs.clear();
        break;
      }

      run.legs.push_back(buildLeg(stops[a], stops[b], e, r));
    }

    if (run.legs.empty()) continue;

    run.totalProfitCr = st.profitCr;
    run.totalRouteDistanceLy = st.routeDistLy;
    run.totalRouteCost = st.routeCost;
    run.totalHops = st.hops;

    run.profitPerLy = run.totalProfitCr / std::max(1e-6, run.totalRouteDistanceLy);
    run.profitPerHop = run.totalProfitCr / (double)std::max(1, run.totalHops);
    run.profitPerCost = run.totalProfitCr / std::max(1e-6, run.totalRouteCost);

    runs.push_back(std::move(run));
  }

  if (runs.empty()) return runs;

  std::sort(runs.begin(), runs.end(), [&](const TradeRun& a, const TradeRun& b) {
    const double sa = scoreTotals(a.totalProfitCr, a.totalRouteDistanceLy, a.totalHops, a.totalRouteCost, p.scoreMode);
    const double sb = scoreTotals(b.totalProfitCr, b.totalRouteDistanceLy, b.totalHops, b.totalRouteCost, p.scoreMode);
    if (std::fabs(sa - sb) > 1e-6) return sa > sb;
    if (std::fabs(a.totalProfitCr - b.totalProfitCr) > 1e-6) return a.totalProfitCr > b.totalProfitCr;

    // Deterministic tiebreak: lexicographic station sequence.
    const std::size_t na = a.legs.size();
    const std::size_t nb = b.legs.size();
    const std::size_t n = std::min(na, nb);
    for (std::size_t i = 0; i < n; ++i) {
      if (a.legs[i].toStation != b.legs[i].toStation) return a.legs[i].toStation < b.legs[i].toStation;
    }
    return na < nb;
  });

  if (runs.size() > p.maxResults) runs.resize(p.maxResults);
  return runs;
}

std::vector<TradeRun> planTradeRunsParallel(core::JobSystem& jobs,
                                            Universe& u,
                                            const SystemStub& originStub,
                                            const Station& originStation,
                                            double timeDays,
                                            double radiusLy,
                                            std::size_t maxSystems,
                                            const TradeRunScanParams& params,
                                            TradeRunFeeRateFn feeRate) {
  auto stubs = u.queryNearbyParallel(jobs, originStub.posLy, radiusLy, maxSystems);
  return planTradeRunsParallel(jobs, u, originStub, originStation, timeDays, stubs, params, std::move(feeRate));
}

} // namespace stellar::sim
