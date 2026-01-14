#include "stellar/proc/TradeRoutes.h"

#include "stellar/econ/Commodity.h"

#include "stellar/proc/HyperlaneRouter.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {

static double clamp01(double x) {
  return std::clamp(x, 0.0, 1.0);
}

static double distLy(const math::Vec3d& a, const math::Vec3d& b) {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  const double dz = a.z - b.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static double marketMass(const TradeProfile& p) {
  // The core "gravity" term for macro trade intensity.
  //
  // Population drives market size, hub drives connectivity, and wealth
  // represents credit throughput.
  const double pop = clamp01(p.population);
  const double hub = clamp01(p.hub);
  const double wealth = clamp01(p.wealth);

  // Keep a non-zero floor so tiny frontier systems still have some signal.
  const double m = 0.18 + 0.62 * pop + 0.32 * hub + 0.18 * wealth;
  return std::max(0.05, m);
}

static double distanceWeight(double distanceLy, const TradeRouteParams& params) {
  const double soft = std::max(0.1, params.distanceSofteningLy);
  const double d = std::max(soft, std::max(0.0, distanceLy));
  const double exp = std::max(0.25, params.distanceExponent);
  return 1.0 / std::pow(d, exp);
}

static double valueDensityWeight(econ::CommodityId cid) {
  // A gentle, game-y bias toward higher value-density goods.
  //
  // We intentionally keep the exponent small so bulk commodities still matter.
  const auto& def = econ::commodityDef(cid);
  const double mass = std::max(0.05, def.massKg);
  const double density = std::max(1.0, def.basePrice / mass);
  return std::pow(density, 0.10); // ~[1.0 .. ~2.0]
}

static std::pair<econ::CommodityId, double> bestAffinity(const std::array<double, econ::kCommodityCount>& a,
                                                        const std::array<double, econ::kCommodityCount>& b) {
  econ::CommodityId bestC = econ::CommodityId::Food;
  double best = -1.0;

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const auto cid = static_cast<econ::CommodityId>(i);
    const double aa = clamp01(a[i]);
    const double bb = clamp01(b[i]);

    const double raw = aa * bb;
    const double w = valueDensityWeight(cid);
    const double s = raw * w;

    if (s > best + 1e-15) {
      best = s;
      bestC = cid;
    }
  }

  if (!std::isfinite(best)) {
    best = 0.0;
    bestC = econ::CommodityId::Food;
  }

  return {bestC, std::max(0.0, best)};
}

static double routeScore(double affinity,
                         const TradeProfile& a,
                         const TradeProfile& b,
                         double distanceLy,
                         const TradeRouteParams& params,
                         bool sameFaction) {
  // Macro trade intensity "gravity" + distance falloff.
  const double g = marketMass(a) * marketMass(b);

  // Hubs amplify routes superlinearly when both ends are strong.
  const double hubA = clamp01(a.hub);
  const double hubB = clamp01(b.hub);
  const double hubBoost = 0.45 + 0.55 * std::sqrt(std::max(0.0, hubA * hubB));

  double s = affinity * g * hubBoost * distanceWeight(distanceLy, params);

  // Small bias for same-faction trade corridors (tariffs/bureaucracy lower).
  if (params.sameFactionBonus && sameFaction) {
    s *= 1.06;
  }

  if (!std::isfinite(s)) return 0.0;
  return std::max(0.0, s);
}

static double riskProxy01(const TradeProfile& a, const TradeProfile& b) {
  // A small proxy for "route risk" that can later drive piracy/insurance.
  const double law = 0.5 * (clamp01(a.lawlessness) + clamp01(b.lawlessness));
  // Hub routes are more patrolled; reduce perceived risk a bit when hub is high.
  const double hub = 0.5 * (clamp01(a.hub) + clamp01(b.hub));
  return clamp01(law * (1.0 - 0.35 * hub));
}

struct RouteGeom {
  bool reachable{true};
  double distanceLy{0.0};
  double directDistanceLy{0.0};
  double laneDistanceLy{0.0};
  int laneHops{0};
  double laneBottleneckBandwidth01{0.0};
  double risk01{0.0};
};

static RouteGeom geometricRouteGeom(const sim::SystemStub& origin,
                                   const TradeProfile& originProfile,
                                   const sim::SystemStub& cand,
                                   const TradeProfile& candProfile) {
  RouteGeom g{};
  g.reachable = true;
  g.directDistanceLy = distLy(origin.posLy, cand.posLy);
  g.distanceLy = g.directDistanceLy;
  g.laneDistanceLy = 0.0;
  g.laneHops = 0;
  g.laneBottleneckBandwidth01 = 0.0;
  g.risk01 = riskProxy01(originProfile, candProfile);
  return g;
}

static RouteGeom hyperlaneRouteGeom(HyperlaneRouter& router,
                                   const sim::SystemStub& origin,
                                   const sim::SystemStub& cand) {
  RouteGeom g{};
  g.directDistanceLy = distLy(origin.posLy, cand.posLy);

  const auto m = router.metricsTo(cand.id);
  if (!m.reachable) {
    g.reachable = false;
    return g;
  }

  g.reachable = true;
  g.distanceLy = m.costLy;
  g.laneDistanceLy = m.distanceLy;
  g.laneHops = m.hops;
  g.laneBottleneckBandwidth01 = m.bottleneckBandwidth01;
  g.risk01 = m.risk01;
  return g;
}

TradeRouteSuggestions suggestTradeRoutes(const sim::SystemStub& origin,
                                        const TradeProfile& originProfile,
                                        const std::vector<sim::SystemStub>& candidates,
                                        const std::vector<TradeProfile>& candidateProfiles,
                                        const TradeRouteParams& params) {
  TradeRouteSuggestions out;

  const std::size_t n = std::min(candidates.size(), candidateProfiles.size());
  out.exports.reserve(std::min<std::size_t>(n, params.maxRoutes));
  out.imports.reserve(std::min<std::size_t>(n, params.maxRoutes));

  for (std::size_t i = 0; i < n; ++i) {
    const auto& cand = candidates[i];
    const auto& p = candidateProfiles[i];

    if (cand.id == 0 || cand.id == origin.id) continue;

    const RouteGeom g = geometricRouteGeom(origin, originProfile, cand, p);
    if (params.maxDistanceLy > 0.0 && g.distanceLy > params.maxDistanceLy) continue;

    const bool sameFaction = (origin.factionId != 0 && origin.factionId == cand.factionId);

    // --- Export (origin -> candidate) --------------------------------------
    {
      const auto [cid, affinity] = bestAffinity(originProfile.exportScore, p.importScore);
      if (!(affinity > 0.0)) {
        // no-op
      } else if (params.dropWeakRoutes && affinity < std::max(0.0, params.minAffinity)) {
        // no-op
      } else {
        TradeRoute r;
        r.otherSystem = cand.id;
        r.commodity = cid;
        r.distanceLy = g.distanceLy;
        r.directDistanceLy = g.directDistanceLy;
        r.laneDistanceLy = g.laneDistanceLy;
        r.laneHops = g.laneHops;
        r.laneBottleneckBandwidth01 = g.laneBottleneckBandwidth01;
        r.risk01 = g.risk01;
        r.score = routeScore(affinity, originProfile, p, g.distanceLy, params, sameFaction);
        if (r.score > 0.0) out.exports.push_back(std::move(r));
      }
    }

    // --- Import (candidate -> origin) --------------------------------------
    {
      const auto [cid, affinity] = bestAffinity(p.exportScore, originProfile.importScore);
      if (!(affinity > 0.0)) {
        // no-op
      } else if (params.dropWeakRoutes && affinity < std::max(0.0, params.minAffinity)) {
        // no-op
      } else {
        TradeRoute r;
        r.otherSystem = cand.id;
        r.commodity = cid;
        r.distanceLy = g.distanceLy;
        r.directDistanceLy = g.directDistanceLy;
        r.laneDistanceLy = g.laneDistanceLy;
        r.laneHops = g.laneHops;
        r.laneBottleneckBandwidth01 = g.laneBottleneckBandwidth01;
        r.risk01 = g.risk01;
        r.score = routeScore(affinity, originProfile, p, g.distanceLy, params, sameFaction);
        if (r.score > 0.0) out.imports.push_back(std::move(r));
      }
    }
  }

  auto cmp = [](const TradeRoute& a, const TradeRoute& b) {
    const double ds = a.score - b.score;
    if (std::abs(ds) > 1e-12) return a.score > b.score;

    const double dd = a.distanceLy - b.distanceLy;
    if (std::abs(dd) > 1e-12) return a.distanceLy < b.distanceLy;

    if (a.otherSystem != b.otherSystem) return a.otherSystem < b.otherSystem;
    return (core::u32)a.commodity < (core::u32)b.commodity;
  };

  std::sort(out.exports.begin(), out.exports.end(), cmp);
  std::sort(out.imports.begin(), out.imports.end(), cmp);

  if (params.maxRoutes > 0) {
    if (out.exports.size() > params.maxRoutes) out.exports.resize(params.maxRoutes);
    if (out.imports.size() > params.maxRoutes) out.imports.resize(params.maxRoutes);
  }

  return out;
}

TradeRouteSuggestions suggestTradeRoutes(const sim::SystemStub& origin,
                                        const TradeProfile& originProfile,
                                        const std::vector<sim::SystemStub>& candidates,
                                        const std::vector<TradeProfile>& candidateProfiles,
                                        const HyperlaneNetwork& hyperlanes,
                                        const TradeRouteParams& params,
                                        const HyperlaneTravelParams& travelParams) {
  TradeRouteSuggestions out;

  const std::size_t n = std::min(candidates.size(), candidateProfiles.size());
  out.exports.reserve(std::min<std::size_t>(n, params.maxRoutes));
  out.imports.reserve(std::min<std::size_t>(n, params.maxRoutes));

  HyperlaneRouter router(hyperlanes, candidates);
  router.compute(origin.id, travelParams);

  for (std::size_t i = 0; i < n; ++i) {
    const auto& cand = candidates[i];
    const auto& p = candidateProfiles[i];

    if (cand.id == 0 || cand.id == origin.id) continue;

    const RouteGeom g = hyperlaneRouteGeom(router, origin, cand);
    if (!g.reachable) continue;
    if (params.maxDistanceLy > 0.0 && g.distanceLy > params.maxDistanceLy) continue;

    const bool sameFaction = (origin.factionId != 0 && origin.factionId == cand.factionId);

    // --- Export (origin -> candidate) --------------------------------------
    {
      const auto [cid, affinity] = bestAffinity(originProfile.exportScore, p.importScore);
      if (!(affinity > 0.0)) {
        // no-op
      } else if (params.dropWeakRoutes && affinity < std::max(0.0, params.minAffinity)) {
        // no-op
      } else {
        TradeRoute r;
        r.otherSystem = cand.id;
        r.commodity = cid;
        r.distanceLy = g.distanceLy;
        r.directDistanceLy = g.directDistanceLy;
        r.laneDistanceLy = g.laneDistanceLy;
        r.laneHops = g.laneHops;
        r.laneBottleneckBandwidth01 = g.laneBottleneckBandwidth01;
        r.risk01 = g.risk01;
        r.score = routeScore(affinity, originProfile, p, g.distanceLy, params, sameFaction);
        if (r.score > 0.0) out.exports.push_back(std::move(r));
      }
    }

    // --- Import (candidate -> origin) --------------------------------------
    {
      const auto [cid, affinity] = bestAffinity(p.exportScore, originProfile.importScore);
      if (!(affinity > 0.0)) {
        // no-op
      } else if (params.dropWeakRoutes && affinity < std::max(0.0, params.minAffinity)) {
        // no-op
      } else {
        TradeRoute r;
        r.otherSystem = cand.id;
        r.commodity = cid;
        r.distanceLy = g.distanceLy;
        r.directDistanceLy = g.directDistanceLy;
        r.laneDistanceLy = g.laneDistanceLy;
        r.laneHops = g.laneHops;
        r.laneBottleneckBandwidth01 = g.laneBottleneckBandwidth01;
        r.risk01 = g.risk01;
        r.score = routeScore(affinity, originProfile, p, g.distanceLy, params, sameFaction);
        if (r.score > 0.0) out.imports.push_back(std::move(r));
      }
    }
  }

  auto cmp = [](const TradeRoute& a, const TradeRoute& b) {
    const double ds = a.score - b.score;
    if (std::abs(ds) > 1e-12) return a.score > b.score;

    const double dd = a.distanceLy - b.distanceLy;
    if (std::abs(dd) > 1e-12) return a.distanceLy < b.distanceLy;

    if (a.otherSystem != b.otherSystem) return a.otherSystem < b.otherSystem;
    return (core::u32)a.commodity < (core::u32)b.commodity;
  };

  std::sort(out.exports.begin(), out.exports.end(), cmp);
  std::sort(out.imports.begin(), out.imports.end(), cmp);

  if (params.maxRoutes > 0) {
    if (out.exports.size() > params.maxRoutes) out.exports.resize(params.maxRoutes);
    if (out.imports.size() > params.maxRoutes) out.imports.resize(params.maxRoutes);
  }

  return out;
}

TradeRouteSuggestions suggestTradeRoutes(core::u64 universeSeed,
                                        const sim::SystemStub& origin,
                                        const std::vector<sim::SystemStub>& candidates,
                                        const TradeRouteParams& params) {
  TradeProfile op = generateTradeProfile(universeSeed, origin);

  std::vector<TradeProfile> cps;
  cps.reserve(candidates.size());
  for (const auto& c : candidates) {
    cps.push_back(generateTradeProfile(universeSeed, c));
  }

  return suggestTradeRoutes(origin, op, candidates, cps, params);
}

TradeRouteSuggestions suggestTradeRoutes(core::u64 universeSeed,
                                        const sim::SystemStub& origin,
                                        const std::vector<sim::SystemStub>& candidates,
                                        const HyperlaneNetwork& hyperlanes,
                                        const TradeRouteParams& params,
                                        const HyperlaneTravelParams& travelParams) {
  TradeProfile op = generateTradeProfile(universeSeed, origin);

  std::vector<TradeProfile> cps;
  cps.reserve(candidates.size());
  for (const auto& c : candidates) {
    cps.push_back(generateTradeProfile(universeSeed, c));
  }

  return suggestTradeRoutes(origin, op, candidates, cps, hyperlanes, params, travelParams);
}

} // namespace stellar::proc
