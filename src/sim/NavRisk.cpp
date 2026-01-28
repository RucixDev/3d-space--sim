#include "stellar/sim/NavRisk.h"

#include "stellar/sim/SystemConditions.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

namespace stellar::sim {

namespace {

inline double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

inline double clamp(double v, double lo, double hi) {
  if (v < lo) return lo;
  if (v > hi) return hi;
  return v;
}

static const SystemSecurityDeltaState* findDelta(SystemId id,
                                                 const std::unordered_map<SystemId, const SystemSecurityDeltaState*>& map) {
  const auto it = map.find(id);
  if (it == map.end()) return nullptr;
  return it->second;
}

static double riskFromProfile(const SystemSecurityProfile& p, const NavRiskParams& params) {
  const double piracy = clamp01(p.piracy01);
  const double insecurity = clamp01(1.0 - p.security01);
  const double lowTraffic = clamp01(1.0 - p.traffic01);
  const double contest = clamp01(p.contest01);

  const double wP = std::max(0.0, params.piracyWeight);
  const double wI = std::max(0.0, params.insecurityWeight);
  const double wT = std::max(0.0, params.lowTrafficWeight);
  const double wC = std::max(0.0, params.contestWeight);

  const double sumW = wP + wI + wT + wC;
  double r = 0.0;
  if (sumW > 1e-12) {
    r = (wP * piracy + wI * insecurity + wT * lowTraffic + wC * contest) / sumW;
  }

  const double lo = std::min(params.minRisk, params.maxRisk);
  const double hi = std::max(params.minRisk, params.maxRisk);
  return clamp(r, lo, hi);
}

} // namespace

NavSystemRiskSignals computeNavSystemRiskSignals(const Universe& universe,
                                                 const StarSystem& sys,
                                                 double timeDays,
                                                 const SystemSecurityDeltaState* delta,
                                                 const NavRiskParams& params) {
  NavSystemRiskSignals out{};

  // Configure whether we want deterministic SystemEvents in the snapshot.
  SystemEventParams evp = params.eventParams;
  if (!params.includeSystemEvents) {
    // Force no events without changing other tuning values.
    evp.eventChance = 0.0;
  }

  const SystemSecurityDeltaState* useDelta = (params.includeSecurityDeltas ? delta : nullptr);

  const auto snap = snapshotSystemConditions(universe.seed(),
                                            sys,
                                            timeDays,
                                            useDelta,
                                            params.dynamicsParams,
                                            evp);

  out.security01 = snap.effective.security01;
  out.piracy01 = snap.effective.piracy01;
  out.traffic01 = snap.effective.traffic01;
  out.contest01 = snap.effective.contest01;

  out.security01 = clamp01(out.security01);
  out.piracy01 = clamp01(out.piracy01);
  out.traffic01 = clamp01(out.traffic01);
  out.contest01 = clamp01(out.contest01);

  out.risk01 = riskFromProfile(snap.effective, params);
  return out;
}

std::vector<double> computeNavRisk01ForNodes(const Universe& universe,
                                             std::span<const SystemStub> nodes,
                                             double timeDays,
                                             std::span<const SystemSecurityDeltaState> securityDeltas,
                                             const NavRiskParams& params) {
  std::vector<double> risk;
  risk.resize(nodes.size(), 0.5);

  // Build a fast delta lookup table.
  std::unordered_map<SystemId, const SystemSecurityDeltaState*> deltaMap;
  if (params.includeSecurityDeltas && !securityDeltas.empty()) {
    deltaMap.reserve(securityDeltas.size());
    for (const auto& d : securityDeltas) {
      if (d.systemId == 0) continue;
      // Keep the first entry for determinism.
      if (deltaMap.find(d.systemId) == deltaMap.end()) {
        deltaMap.emplace(d.systemId, &d);
      }
    }
  }

  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const SystemId id = nodes[i].id;
    if (id == 0) {
      risk[i] = 0.5;
      continue;
    }

    const auto* delta = findDelta(id, deltaMap);
    const StarSystem& sys = universe.getSystem(id, &nodes[i]);
    const auto sig = computeNavSystemRiskSignals(universe, sys, timeDays, delta, params);
    risk[i] = clamp01(sig.risk01);
  }

  return risk;
}

double estimateRouteAvgRisk01(std::span<const SystemStub> nodes,
                              std::span<const SystemId> route,
                              std::span<const double> risk01PerNode) {
  if (route.empty() || nodes.empty()) return 0.0;

  std::unordered_map<SystemId, std::size_t> idx;
  idx.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (nodes[i].id != 0) idx[nodes[i].id] = i;
  }

  auto riskAt = [&](std::size_t i) -> double {
    if (i >= nodes.size()) return 0.0;
    if (risk01PerNode.size() != nodes.size()) return 0.0;
    return clamp01(risk01PerNode[i]);
  };

  if (route.size() == 1) {
    const auto it = idx.find(route[0]);
    if (it == idx.end()) return 0.0;
    return riskAt(it->second);
  }

  double weighted = 0.0;
  double totalD = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) continue;
    const std::size_t ia = itA->second;
    const std::size_t ib = itB->second;
    const double d = (nodes[ia].posLy - nodes[ib].posLy).length();
    if (!(d > 1e-12)) continue;
    const double r = 0.5 * (riskAt(ia) + riskAt(ib));
    weighted += r * d;
    totalD += d;
  }
  if (!(totalD > 1e-12)) return 0.0;
  return clamp01(weighted / totalD);
}

double estimateRouteMaxRisk01(std::span<const SystemStub> nodes,
                              std::span<const SystemId> route,
                              std::span<const double> risk01PerNode) {
  if (route.empty() || nodes.empty()) return 0.0;
  if (risk01PerNode.size() != nodes.size()) return 0.0;

  std::unordered_map<SystemId, std::size_t> idx;
  idx.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (nodes[i].id != 0) idx[nodes[i].id] = i;
  }

  double mx = 0.0;
  for (const auto id : route) {
    const auto it = idx.find(id);
    if (it == idx.end()) continue;
    const std::size_t i = it->second;
    mx = std::max(mx, clamp01(risk01PerNode[i]));
  }
  return clamp01(mx);
}

} // namespace stellar::sim
