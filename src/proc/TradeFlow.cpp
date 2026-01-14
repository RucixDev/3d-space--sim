#include "stellar/proc/TradeFlow.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace stellar::proc {

static inline double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

static inline double safePow(double a, double b) {
  const double x = std::max(0.0, a);
  const double y = std::max(0.0, b);
  const double r = std::pow(x, y);
  if (!std::isfinite(r)) return std::numeric_limits<double>::infinity();
  return r;
}

struct EdgeKey {
  sim::SystemId a{0};
  sim::SystemId b{0};
};

struct EdgeKeyHash {
  std::size_t operator()(const EdgeKey& k) const noexcept {
    // Mix into 64-bit then truncate.
    core::u64 h = 0;
    h = core::hashCombine(h, static_cast<core::u64>(k.a));
    h = core::hashCombine(h, static_cast<core::u64>(k.b));
    return static_cast<std::size_t>(h);
  }
};

struct EdgeKeyEq {
  bool operator()(const EdgeKey& x, const EdgeKey& y) const noexcept {
    return x.a == y.a && x.b == y.b;
  }
};

static inline EdgeKey makeKey(sim::SystemId u, sim::SystemId v) {
  if (u < v) return EdgeKey{u, v};
  return EdgeKey{v, u};
}

static double economicMass(const TradeProfile& p, const TradeFlowParams& params) {
  const double hub = clamp01(p.hub);
  const double pop = clamp01(p.population);
  const double wea = clamp01(p.wealth);
  const double ind = clamp01(p.industry);

  const double m = params.massFloor
                   + params.hubWeight * hub
                   + params.populationWeight * pop
                   + params.wealthWeight * wea
                   + params.industryWeight * ind;
  return std::clamp(m, 0.01, 3.0);
}

static double commodityCompatibility(const TradeProfile& a, const TradeProfile& b) {
  // Symmetric match: exports of A with imports of B and vice versa.
  double sum = 0.0;
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    sum += clamp01(a.exportScore[i]) * clamp01(b.importScore[i]);
    sum += clamp01(b.exportScore[i]) * clamp01(a.importScore[i]);
  }
  const double denom = std::max(1.0, 2.0 * (double)econ::kCommodityCount);
  return std::clamp(sum / denom, 0.0, 1.0);
}

static double pairPotential(const TradeProfile& a, const TradeProfile& b,
                            double massA, double massB,
                            const TradeFlowParams& params) {
  const double comp = commodityCompatibility(a, b);
  const double w = std::clamp(params.commodityCompatibilityWeight, 0.0, 1.0);
  const double compF = (1.0 - w) + w * comp;
  return std::max(0.0, massA * massB * compF);
}

static std::vector<int> pickSampleSources(const std::vector<double>& masses,
                                          const std::vector<sim::SystemId>& ids,
                                          const TradeFlowParams& params) {
  const int n = (int)masses.size();
  std::vector<int> idx(n);
  for (int i = 0; i < n; ++i) idx[i] = i;

  std::sort(idx.begin(), idx.end(), [&](int a, int b) {
    if (masses[a] != masses[b]) return masses[a] > masses[b];
    return ids[a] < ids[b];
  });

  const int want = std::clamp(params.sampleSourceCount, 1, n);
  idx.resize((std::size_t)want);
  std::sort(idx.begin(), idx.end()); // group Dijkstra runs by stable node index
  return idx;
}

TradeFlowNetwork computeTradeFlow(core::u64 universeSeed,
                                  const std::vector<sim::SystemStub>& stubs,
                                  const std::vector<TradeProfile>& profiles,
                                  const HyperlaneNetwork& hyperlanes,
                                  const TradeFlowParams& params) {
  TradeFlowNetwork out{};

  if (stubs.empty() || profiles.empty()) return out;
  if (stubs.size() != profiles.size()) return out;

  // Build a stable node list (sorted by id) so results are deterministic even
  // if the caller passes differently ordered stubs.
  struct Node {
    sim::SystemId id{0};
    TradeProfile profile{};
    double mass{0.0};
  };

  std::vector<Node> nodes;
  nodes.reserve(stubs.size());
  for (std::size_t i = 0; i < stubs.size(); ++i) {
    const auto id = stubs[i].id;
    if (id == 0) continue;
    Node n{};
    n.id = id;
    n.profile = profiles[i];
    n.mass = economicMass(n.profile, params);
    nodes.push_back(n);
  }

  if (nodes.size() < 2) return out;

  std::sort(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
    return a.id < b.id;
  });

  const int n = (int)nodes.size();
  std::unordered_map<sim::SystemId, int> idToNode;
  idToNode.reserve(nodes.size() * 2);
  for (int i = 0; i < n; ++i) idToNode[nodes[(std::size_t)i].id] = i;

  // Collect stable arrays.
  std::vector<sim::SystemId> ids;
  std::vector<double> masses;
  ids.reserve(nodes.size());
  masses.reserve(nodes.size());
  for (const auto& nd : nodes) {
    ids.push_back(nd.id);
    masses.push_back(nd.mass);
  }

  // Build routing helper.
  HyperlaneRouter router(hyperlanes, stubs);

  // Determine which pairs to consider.
  std::vector<std::pair<int, int>> pairs;
  pairs.reserve((std::size_t)n * 8);

  if (nodes.size() <= params.fullPairLimit) {
    // Full all-pairs (i < j).
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        pairs.emplace_back(i, j);
      }
    }
  } else {
    // Sampled approximation.
    const std::vector<int> sources = pickSampleSources(masses, ids, params);

    // Add direct hyperlane neighbors for each source (captures local corridor flow).
    std::unordered_map<sim::SystemId, std::vector<sim::SystemId>> neigh;
    neigh.reserve(hyperlanes.edges.size());
    for (const auto& e : hyperlanes.edges) {
      neigh[e.a].push_back(e.b);
      neigh[e.b].push_back(e.a);
    }
    for (auto& kv : neigh) {
      auto& v = kv.second;
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    for (int si : sources) {
      const sim::SystemId sid = ids[(std::size_t)si];

      // Always include direct neighbors (if present in the node set).
      if (auto it = neigh.find(sid); it != neigh.end()) {
        for (sim::SystemId nb : it->second) {
          auto itN = idToNode.find(nb);
          if (itN == idToNode.end()) continue;
          const int j = itN->second;
          if (j == si) continue;
          const int a = std::min(si, j);
          const int b = std::max(si, j);
          pairs.emplace_back(a, b);
        }
      }

      // Add random partners across the node set.
      const int k = std::clamp(params.samplePartnersPerSource, 0, 4096);
      core::u64 s = core::hashCombine(universeSeed, core::fnv1a64("trade_flow_pairs_v1"));
      s = core::hashCombine(s, static_cast<core::u64>(sid));
      for (int t = 0; t < k; ++t) {
        const core::u64 h = core::hashCombine(s, static_cast<core::u64>(t));
        const int j = (int)(h % (core::u64)n);
        if (j == si) continue;
        const int a = std::min(si, j);
        const int b = std::max(si, j);
        pairs.emplace_back(a, b);
      }
    }

    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
  }

  // Group targets by source so we can reuse Dijkstra results.
  std::vector<std::vector<int>> targets((std::size_t)n);
  for (const auto& pr : pairs) {
    const int a = pr.first;
    const int b = pr.second;
    if (a < 0 || b < 0 || a >= n || b >= n) continue;
    if (a == b) continue;
    targets[(std::size_t)a].push_back(b);
  }
  for (auto& v : targets) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }

  // Accumulate routed flow per hyperlane edge.
  std::unordered_map<EdgeKey, double, EdgeKeyHash, EdgeKeyEq> edgeFlow;
  edgeFlow.reserve(hyperlanes.edges.size() * 2 + 64);

  double totalFlow = 0.0;
  std::vector<sim::SystemId> path;
  path.reserve(64);

  const double exp = std::max(0.01, params.gravityExponent);
  const double minCost = std::max(1e-6, params.minCostLy);

  for (int i = 0; i < n; ++i) {
    if (targets[(std::size_t)i].empty()) continue;

    const sim::SystemId originId = ids[(std::size_t)i];
    if (!router.compute(originId, params.travelParams)) continue;

    const Node& A = nodes[(std::size_t)i];
    const double massA = masses[(std::size_t)i];

    for (int j : targets[(std::size_t)i]) {
      if (j <= i || j >= n) continue;
      const sim::SystemId targetId = ids[(std::size_t)j];

      const auto met = router.metricsTo(targetId);
      if (!met.reachable) continue;

      const double cost = std::max(minCost, met.costLy);
      const Node& B = nodes[(std::size_t)j];
      const double massB = masses[(std::size_t)j];

      const double pot = pairPotential(A.profile, B.profile, massA, massB, params);
      if (pot <= 0.0) continue;

      const double denom = safePow(cost, exp);
      if (!std::isfinite(denom) || denom <= 1e-12) continue;

      const double f = pot / denom;
      if (!std::isfinite(f) || f <= 0.0) continue;

      if (!router.buildPathTo(targetId, path)) continue;
      if (path.size() < 2) continue;

      totalFlow += f;

      for (std::size_t k = 1; k < path.size(); ++k) {
        const sim::SystemId u = path[k - 1];
        const sim::SystemId v = path[k];
        if (u == 0 || v == 0 || u == v) continue;
        const EdgeKey key = makeKey(u, v);
        edgeFlow[key] += f;
      }
    }
  }

  // Compute per-node incident traffic by summing incident edge flow.
  std::vector<double> nodeTraffic((std::size_t)n, 0.0);
  double maxEdgeFlow = 0.0;

  for (const auto& kv : edgeFlow) {
    const EdgeKey& k = kv.first;
    const double f = kv.second;
    if (!std::isfinite(f) || f <= 0.0) continue;

    maxEdgeFlow = std::max(maxEdgeFlow, f);

    auto ita = idToNode.find(k.a);
    auto itb = idToNode.find(k.b);
    if (ita != idToNode.end()) nodeTraffic[(std::size_t)ita->second] += f;
    if (itb != idToNode.end()) nodeTraffic[(std::size_t)itb->second] += f;
  }

  double maxTraffic = 0.0;
  for (double t : nodeTraffic) maxTraffic = std::max(maxTraffic, t);

  // Build outputs.
  out.totalFlow = totalFlow;
  out.nodes.reserve(nodes.size());
  for (int i = 0; i < n; ++i) {
    TradeFlowNode tn{};
    tn.id = ids[(std::size_t)i];
    tn.traffic = std::max(0.0, nodeTraffic[(std::size_t)i]);
    tn.traffic01 = (maxTraffic > 1e-12) ? std::clamp(tn.traffic / maxTraffic, 0.0, 1.0) : 0.0;
    out.nodes.push_back(tn);
  }

  out.edges.reserve(edgeFlow.size());
  for (const auto& kv : edgeFlow) {
    const EdgeKey& k = kv.first;
    const double f = kv.second;
    if (!std::isfinite(f) || f <= 0.0) continue;
    TradeFlowEdge te{};
    te.a = k.a;
    te.b = k.b;
    if (te.a > te.b) std::swap(te.a, te.b);
    te.flow = f;
    te.flow01 = (maxEdgeFlow > 1e-12) ? std::clamp(f / maxEdgeFlow, 0.0, 1.0) : 0.0;
    out.edges.push_back(te);
  }

  std::sort(out.nodes.begin(), out.nodes.end(), [](const TradeFlowNode& a, const TradeFlowNode& b) {
    if (a.traffic != b.traffic) return a.traffic > b.traffic;
    return a.id < b.id;
  });

  std::sort(out.edges.begin(), out.edges.end(), [](const TradeFlowEdge& a, const TradeFlowEdge& b) {
    if (a.flow != b.flow) return a.flow > b.flow;
    if (a.a != b.a) return a.a < b.a;
    return a.b < b.b;
  });

  return out;
}

TradeFlowNetwork computeTradeFlow(core::u64 universeSeed,
                                  const std::vector<sim::SystemStub>& stubs,
                                  const HyperlaneNetwork& hyperlanes,
                                  const TradeFlowParams& params) {
  std::vector<TradeProfile> profiles;
  profiles.reserve(stubs.size());
  for (const auto& s : stubs) {
    profiles.push_back(generateTradeProfile(universeSeed, s));
  }
  return computeTradeFlow(universeSeed, stubs, profiles, hyperlanes, params);
}

} // namespace stellar::proc
