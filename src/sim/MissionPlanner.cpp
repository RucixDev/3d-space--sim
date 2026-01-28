#include "stellar/sim/MissionPlanner.h"

#include "stellar/sim/MissionBriefing.h"
#include "stellar/sim/NavRouteBatch.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace stellar::sim {

namespace {

inline double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

inline double safeHoursLeft(double deadlineDay, double nowDay) {
  if (deadlineDay <= 0.0) return std::numeric_limits<double>::infinity();
  const double d = (deadlineDay - nowDay) * 24.0;
  return d;
}

inline double safeNonNeg(double v) {
  if (!std::isfinite(v)) return 0.0;
  return std::max(0.0, v);
}

inline double etaStopOverheadSec(const MissionObjective& obj, const MissionItineraryParams& params) {
  // Sites usually imply "scan the signal" or "reach the area" rather than a full
  // docking/menu flow, so allow a separate overhead knob.
  if (obj.isSite) {
    const double s = safeNonNeg(params.etaSecondsPerSite);
    if (s > 1e-6) return s;
  }
  return safeNonNeg(params.etaSecondsPerStop);
}

inline double etaTravelDays(int hops, double distLy, const MissionObjective& obj, const MissionItineraryParams& params) {
  const double perJump = safeNonNeg(params.etaSecondsPerJump);
  const double perLy = safeNonNeg(params.etaSecondsPerLy);
  const double tSec = (double)std::max(0, hops) * perJump + safeNonNeg(distLy) * perLy + etaStopOverheadSec(obj, params);
  return tSec / 86400.0;
}

inline double repForFaction(std::span<const FactionReputation> rep, core::u32 factionId) {
  if (factionId == 0) return 0.0;
  for (const auto& r : rep) {
    if (r.factionId == factionId) return r.rep;
  }
  return 0.0;
}

const Station* findStationById(const StarSystem& sys, StationId id) {
  if (id == 0) return nullptr;
  for (const auto& st : sys.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

struct ObjKey {
  SystemId sys{0};
  StationId st{0};
  bool site{false};
};

struct ObjKeyLess {
  bool operator()(const ObjKey& a, const ObjKey& b) const {
    if (a.sys != b.sys) return a.sys < b.sys;
    if (a.st != b.st) return a.st < b.st;
    return (int)a.site < (int)b.site;
  }
};

struct Agg {
  MissionObjective obj{};
  std::vector<core::u64> missionIds;
  double rewardCr{0.0};
  double riskSum{0.0};
  int riskN{0};
  double earliestDeadlineDay{0.0};
};

} // namespace

MissionObjective missionNextObjective(const Mission& m) {
  MissionObjective o{};

  // Multi-leg missions first.
  const bool multiLeg = (m.type == MissionType::MultiDelivery) || (m.type == MissionType::Passenger);
  if (multiLeg && m.viaSystem != 0 && m.viaStation != 0 && m.leg == 0) {
    o.systemId = m.viaSystem;
    o.stationId = m.viaStation;
    o.isSite = false;
    return o;
  }

  // Salvage: the first objective is a mission site inside the destination system.
  if (m.type == MissionType::Salvage && !m.scanned && m.targetNpcId != 0) {
    o.systemId = m.toSystem;
    o.stationId = 0;
    o.isSite = true;
    return o;
  }

  // Bounties are system-level (no station required).
  if (m.type == MissionType::BountyScan || m.type == MissionType::BountyKill) {
    o.systemId = m.toSystem;
    o.stationId = 0;
    o.isSite = true;
    return o;
  }

  // Default: dock at destination.
  o.systemId = m.toSystem;
  o.stationId = m.toStation;
  o.isSite = (o.stationId == 0);
  return o;
}

MissionItineraryResult planMissionItinerary(Universe& universe,
                                           const StarSystem& currentSystem,
                                           double timeDays,
                                           std::span<const Mission> missions,
                                           std::span<const FactionReputation> playerRepWithFaction,
                                           std::span<const SystemSecurityDeltaState> securityDeltas,
                                           const MissionItineraryParams& params) {
  MissionItineraryResult out{};
  out.startSystemId = currentSystem.stub.id;
  out.ok = true;

  if (!std::isfinite(timeDays)) timeDays = 0.0;

  const double maxJumpLy = std::max(0.0, params.maxJumpLy);
  const std::size_t maxSys = std::max<std::size_t>(1, params.maxSystems);
  const int maxStops = std::clamp(params.maxStops, 1, 64);

  // 1) Aggregate missions into objective stops.
  std::map<ObjKey, Agg, ObjKeyLess> aggs;
  for (const auto& m : missions) {
    if (m.completed || m.failed) continue;
    if (m.toSystem == 0) continue;

    MissionObjective obj = missionNextObjective(m);
    if (obj.systemId == 0) continue;

    ObjKey key{};
    key.sys = obj.systemId;
    key.site = obj.isSite;
    if (params.groupBySystem) {
      key.st = 0;
    } else {
      key.st = obj.stationId;
    }

    Agg& a = aggs[key];
    if (a.obj.systemId == 0) {
      a.obj = obj;
    }

    // If grouping by system, attempt to keep a stable stationId only when all
    // missions require docking at the same station.
    if (params.groupBySystem) {
      // Keep the stop as a docking stop only if *every* mission objective in the
      // system is a dock at the *same* station. Any non-dock objective (site / system)
      // collapses the stop to system-level.
      if (obj.isSite || obj.stationId == 0) {
        a.obj.stationId = 0;
        a.obj.isSite = true;
      } else {
        if (a.obj.isSite) {
          a.obj.stationId = 0;
          a.obj.isSite = true;
        } else if (a.obj.stationId == 0) {
          a.obj.stationId = obj.stationId;
          a.obj.isSite = false;
        } else if (a.obj.stationId != obj.stationId) {
          a.obj.stationId = 0;
          a.obj.isSite = true;
        }
      }
    }

    a.missionIds.push_back(m.id);
    a.rewardCr += std::max(0.0, m.reward);

    if (m.deadlineDay > 0.0) {
      if (a.earliestDeadlineDay <= 0.0 || m.deadlineDay < a.earliestDeadlineDay) {
        a.earliestDeadlineDay = m.deadlineDay;
      }
    }

    // Risk model (best-effort; falls back to 0.5 if origin station isn't available).
    double risk01 = 0.5;
    {
      const StarSystem* originSys = nullptr;
      const Station* originSt = nullptr;

      if (m.fromSystem != 0) {
        originSys = &universe.getSystem(m.fromSystem);
        originSt = findStationById(*originSys, m.fromStation);
        if (!originSt && !originSys->stations.empty()) originSt = &originSys->stations.front();
      }

      if (!originSys) {
        originSys = &currentSystem;
      }
      Station dummy{};
      if (!originSt) {
        if (!originSys->stations.empty()) originSt = &originSys->stations.front();
        else originSt = &dummy;
      }

      const double rep = repForFaction(playerRepWithFaction, m.factionId);
      MissionBriefingParams bp{};
      bp.useMarkup = false;
      bp.includeRiskHints = false;
      const auto r = computeMissionRisk(universe, *originSys, *originSt, timeDays, rep, m, securityDeltas, bp);
      risk01 = clamp01(r.overall01);
    }

    a.riskSum += risk01;
    a.riskN += 1;
  }

  if (aggs.empty()) {
    out.ok = true;
    return out;
  }

  // Create a stable candidate list.
  std::vector<Agg> candidates;
  candidates.reserve(aggs.size());
  for (auto& kv : aggs) {
    Agg a = kv.second;
    if (a.riskN > 0) a.riskSum /= (double)a.riskN;
    // Stable, deterministic ordering of mission ids within each stop.
    std::sort(a.missionIds.begin(), a.missionIds.end());
    candidates.push_back(std::move(a));
  }

  // Greedy selection.
  std::vector<bool> used(candidates.size(), false);
  SystemId curId = currentSystem.stub.id;
  const auto& curSys = universe.getSystem(curId, &currentSystem.stub);
  (void)curSys;

  // Deadline scoring needs to account for the fact that earlier itinerary legs
  // consume time. We model this using a simple ETA accumulator (in simulation
  // days) rather than assuming all stops can be reached "right now".
  double etaNowDay = timeDays;

  for (int step = 0; step < maxStops; ++step) {
    // Find any remaining.
    bool any = false;
    for (bool u : used) {
      if (!u) { any = true; break; }
    }
    if (!any) break;

    // Automatic radius based on remaining straight-line distances.
    double rLy = params.queryRadiusLy;
    if (rLy <= 1e-6) {
      const auto& curStub = universe.getSystem(curId).stub;
      double maxDirect = 0.0;
      for (std::size_t i = 0; i < candidates.size(); ++i) {
        if (used[i]) continue;
        const auto& dstStub = universe.getSystem(candidates[i].obj.systemId).stub;
        const double d = (dstStub.posLy - curStub.posLy).length();
        if (d > maxDirect) maxDirect = d;
      }
      // Keep it comfortably above the farthest straight-line target.
      rLy = std::clamp(maxDirect * 1.25 + 60.0, 180.0, 2400.0);
    }

    // Batch routes from the current system to all nearby nodes.
    auto nodes = universe.queryNearby(universe.getSystem(curId).stub.posLy, rLy, maxSys);
    bool hasStart = false;
    for (const auto& s : nodes) {
      if (s.id == curId) { hasStart = true; break; }
    }
    if (!hasStart) {
      nodes.push_back(universe.getSystem(curId).stub);
    }

    const auto batch = computeNavRouteBatchCost(nodes,
                                               curId,
                                               maxJumpLy,
                                               params.costPerJump,
                                               params.costPerLy);

    // Choose the best candidate.
    int bestIdx = -1;
    double bestScore = -1.0;
    double bestCost = std::numeric_limits<double>::infinity();
    int bestHops = 0;
    double bestDist = 0.0;
    bool bestReach = true;

    for (std::size_t i = 0; i < candidates.size(); ++i) {
      if (used[i]) continue;
      const Agg& c = candidates[i];

      const SystemId dstId = c.obj.systemId;
      bool reach = true;
      double cost = 0.0;
      int hops = 0;
      double dist = 0.0;

      if (dstId != curId) {
        if (batch.reachable(dstId)) {
          cost = batch.costTo(dstId);
          hops = batch.hopsTo(dstId);
          dist = batch.distanceTo(dstId);
        } else {
          reach = false;
          const auto& curStub = universe.getSystem(curId).stub;
          const auto& dstStub = universe.getSystem(dstId).stub;
          dist = (dstStub.posLy - curStub.posLy).length();
          const int estHops = (maxJumpLy > 1e-6) ? (int)std::ceil(dist / maxJumpLy) : 999999;
          hops = std::max(1, estHops);
          cost = (double)hops * params.costPerJump + dist * params.costPerLy;
          // Penalize unreachable estimates so they don't dominate.
          cost *= 2.25;
        }
      }

      // Score = reward adjusted by risk/urgency over travel cost.
      const double reward = std::max(0.0, c.rewardCr) * std::max(0.0, params.rewardWeight);
      const double riskFactor = std::clamp(1.0 - std::max(0.0, params.riskWeight) * clamp01(c.riskSum), 0.15, 1.0);

      double urgency = 1.0;
      const double hrsLeft = safeHoursLeft(c.earliestDeadlineDay, etaNowDay);
      if (std::isfinite(hrsLeft)) {
        double slackHrs = hrsLeft;
        if (params.etaAwareUrgency) {
          const double tDays = etaTravelDays(hops, dist, c.obj, params);
          slackHrs = hrsLeft - tDays * 24.0;
        }

        if (slackHrs <= 0.0) urgency = 0.05;
        else urgency = 1.0 + std::max(0.0, params.urgencyWeight) * (1.0 / (0.5 + slackHrs));
      }

      const double denom = 1.0 + std::max(0.0, cost);
      const double score = (reward * riskFactor * urgency) / denom;

      // Deterministic tie-breaking.
      const bool better =
        (score > bestScore + 1e-12) ||
        (std::abs(score - bestScore) <= 1e-12 && cost < bestCost - 1e-9) ||
        (std::abs(score - bestScore) <= 1e-12 && std::abs(cost - bestCost) <= 1e-9 && c.obj.systemId < (bestIdx >= 0 ? candidates[(std::size_t)bestIdx].obj.systemId : std::numeric_limits<SystemId>::max()));

      if (better) {
        bestIdx = (int)i;
        bestScore = score;
        bestCost = cost;
        bestHops = hops;
        bestDist = dist;
        bestReach = reach;
      }
    }

    if (bestIdx < 0) break;

    used[(std::size_t)bestIdx] = true;

    MissionItineraryStop s{};
    s.objective = candidates[(std::size_t)bestIdx].obj;
    s.missionIds = candidates[(std::size_t)bestIdx].missionIds;
    s.rewardCr = candidates[(std::size_t)bestIdx].rewardCr;
    s.avgRisk01 = clamp01(candidates[(std::size_t)bestIdx].riskSum);
    s.earliestDeadlineDay = candidates[(std::size_t)bestIdx].earliestDeadlineDay;
    s.reachable = bestReach;
    s.hopsFromPrev = bestHops;
    s.distanceLyFromPrev = bestDist;
    s.costFromPrev = bestCost;

    // ETA analytics for this selected leg.
    s.etaTravelDaysFromPrev = etaTravelDays(bestHops, bestDist, s.objective, params);
    etaNowDay += s.etaTravelDaysFromPrev;
    s.etaDay = etaNowDay;
    if (s.earliestDeadlineDay > 0.0) {
      s.etaSlackHours = (s.earliestDeadlineDay - s.etaDay) * 24.0;
    } else {
      s.etaSlackHours = std::numeric_limits<double>::infinity();
    }

    out.totalCost += bestCost;
    out.totalHops += bestHops;
    out.totalDistanceLy += bestDist;
    if (!bestReach) out.unreachableStops += 1;

    out.stops.push_back(std::move(s));
    curId = out.stops.back().objective.systemId;
  }

  // If we couldn't emit any stops, still report ok (the caller can interpret).
  out.ok = true;
  return out;
}

} // namespace stellar::sim
