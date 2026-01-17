#include "stellar/proc/HyperlaneKRoutes.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace stellar::proc {

namespace {

static double clamp01(double x) {
  return std::clamp(x, 0.0, 1.0);
}

static double edgeCostLy(double distanceLy, double bandwidth01, double risk01, const HyperlaneTravelParams& p) {
  const double dist = std::max(0.0, distanceLy);

  const double risk = clamp01(risk01);
  const double riskMul = 1.0 + std::max(0.0, p.riskWeight) * risk;

  const double bw = clamp01(bandwidth01);
  const double minF = std::clamp(p.minBandwidthFactor, 0.05, 1.0);
  const double bwF = minF + (1.0 - minF) * bw;
  const double bias = clamp01(p.bandwidthBias);
  const double denom = (1.0 - bias) + bias * bwF;

  const double c = (dist * riskMul) / std::max(0.05, denom);
  if (!std::isfinite(c)) return dist;
  return std::max(0.0, c);
}

static std::size_t hashCombine(std::size_t h, std::size_t v) {
  // A small, decent hash combine (boost-style).
  return h ^ (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
}

struct Edge {
  sim::SystemId from{0};
  sim::SystemId to{0};
  bool operator==(const Edge& o) const { return from == o.from && to == o.to; }
};

struct EdgeHash {
  std::size_t operator()(const Edge& e) const {
    std::size_t h = 0;
    h = hashCombine(h, std::hash<sim::SystemId>()(e.from));
    h = hashCombine(h, std::hash<sim::SystemId>()(e.to));
    return h;
  }
};

static bool pathHasPrefix(const std::vector<sim::SystemId>& path, const std::vector<sim::SystemId>& prefix) {
  if (prefix.size() > path.size()) return false;
  return std::equal(prefix.begin(), prefix.end(), path.begin());
}

static bool pathLexLess(const std::vector<sim::SystemId>& a, const std::vector<sim::SystemId>& b) {
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

struct Adj {
  std::size_t to{0};
  double distanceLy{0.0};
  double bandwidth01{0.0};
  double risk01{0.0};
};

struct Graph {
  std::vector<sim::SystemId> ids;
  std::vector<math::Vec3d> posLy;
  std::unordered_map<sim::SystemId, std::size_t> idx;
  std::vector<std::vector<Adj>> adj;
};

static Graph buildGraph(const std::vector<sim::SystemStub>& nodes, const HyperlaneNetwork& net) {
  Graph g{};

  g.ids.reserve(nodes.size());
  std::unordered_map<sim::SystemId, math::Vec3d> posMap;
  posMap.reserve(nodes.size() * 2);

  for (const auto& s : nodes) {
    if (s.id == 0) continue;
    g.ids.push_back(s.id);
    posMap[s.id] = s.posLy;
  }

  std::sort(g.ids.begin(), g.ids.end());
  g.ids.erase(std::unique(g.ids.begin(), g.ids.end()), g.ids.end());

  g.idx.reserve(g.ids.size() * 2);
  g.posLy.reserve(g.ids.size());

  for (std::size_t i = 0; i < g.ids.size(); ++i) {
    const sim::SystemId id = g.ids[i];
    g.idx[id] = i;
    auto it = posMap.find(id);
    g.posLy.push_back(it != posMap.end() ? it->second : math::Vec3d{0.0, 0.0, 0.0});
  }

  g.adj.assign(g.ids.size(), {});

  for (const auto& e : net.edges) {
    auto itA = g.idx.find(e.a);
    auto itB = g.idx.find(e.b);
    if (itA == g.idx.end() || itB == g.idx.end()) continue;
    const std::size_t ia = itA->second;
    const std::size_t ib = itB->second;
    if (ia == ib) continue;

    Adj a2b;
    a2b.to = ib;
    a2b.distanceLy = e.distanceLy;
    a2b.bandwidth01 = e.bandwidth01;
    a2b.risk01 = e.risk01;
    g.adj[ia].push_back(a2b);

    Adj b2a = a2b;
    b2a.to = ia;
    g.adj[ib].push_back(b2a);
  }

  // Deterministic neighbor iteration.
  for (auto& v : g.adj) {
    std::sort(v.begin(), v.end(), [&](const Adj& a, const Adj& b) {
      if (a.to != b.to) return g.ids[a.to] < g.ids[b.to];
      if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
      if (a.bandwidth01 != b.bandwidth01) return a.bandwidth01 > b.bandwidth01;
      return a.risk01 < b.risk01;
    });
  }

  return g;
}

struct AStarSolveResult {
  std::vector<sim::SystemId> path;
  double cost{0.0};
  int expansions{0};
  bool reached{false};
};

static const Adj* findAdj(const Graph& g, std::size_t u, std::size_t v) {
  if (u >= g.adj.size()) return nullptr;
  for (const auto& a : g.adj[u]) {
    if (a.to == v) return &a;
  }
  return nullptr;
}

static HyperlanePathMetrics computeMetrics(const Graph& g,
                                          const std::vector<sim::SystemId>& path,
                                          const HyperlaneTravelParams& travel) {
  HyperlanePathMetrics m{};
  if (path.empty()) return m;
  if (path.size() == 1) {
    m.reachable = true;
    m.hops = 0;
    m.costLy = 0.0;
    m.distanceLy = 0.0;
    m.risk01 = 0.0;
    m.bottleneckBandwidth01 = 1.0;
    return m;
  }

  double cost = 0.0;
  double dist = 0.0;
  double riskNot = 1.0;
  double bottleneck = 1.0;

  for (std::size_t i = 0; i + 1 < path.size(); ++i) {
    auto itA = g.idx.find(path[i]);
    auto itB = g.idx.find(path[i + 1]);
    if (itA == g.idx.end() || itB == g.idx.end()) {
      return HyperlanePathMetrics{};
    }
    const std::size_t ia = itA->second;
    const std::size_t ib = itB->second;

    const Adj* a = findAdj(g, ia, ib);
    if (!a) return HyperlanePathMetrics{};

    const double d = std::max(0.0, a->distanceLy);
    const double bw = clamp01(a->bandwidth01);
    const double risk = clamp01(a->risk01);

    dist += d;
    cost += edgeCostLy(d, bw, risk, travel);
    bottleneck = std::min(bottleneck, bw);
    riskNot *= (1.0 - risk);
  }

  m.reachable = true;
  m.hops = (int)path.size() - 1;
  m.distanceLy = std::max(0.0, dist);
  m.costLy = std::max(0.0, cost);
  m.bottleneckBandwidth01 = clamp01(bottleneck);
  m.risk01 = clamp01(1.0 - std::clamp(riskNot, 0.0, 1.0));
  return m;
}

static AStarSolveResult aStarSolve(const Graph& g,
                                  sim::SystemId startId,
                                  sim::SystemId goalId,
                                  const HyperlaneTravelParams& travel,
                                  const std::unordered_set<sim::SystemId>* bannedNodes,
                                  const std::unordered_set<Edge, EdgeHash>* bannedEdges,
                                  std::size_t maxExpansions) {
  AStarSolveResult out{};

  if (startId == 0 || goalId == 0) return out;
  if (g.ids.empty()) return out;

  if (bannedNodes) {
    if (bannedNodes->count(startId) != 0) return out;
    if (bannedNodes->count(goalId) != 0) return out;
  }

  auto itS = g.idx.find(startId);
  auto itG = g.idx.find(goalId);
  if (itS == g.idx.end() || itG == g.idx.end()) return out;

  const std::size_t start = itS->second;
  const std::size_t goal = itG->second;
  const std::size_t N = g.ids.size();

  if (start == goal) {
    out.reached = true;
    out.path = {startId};
    out.cost = 0.0;
    out.expansions = 0;
    return out;
  }

  std::vector<int> cameFrom(N, -1);
  std::vector<double> gScore(N, std::numeric_limits<double>::infinity());
  std::vector<char> closed(N, 0);

  // Admissible heuristic: straight-line distance to the goal.
  // edgeCostLy(distance, ...) >= distance, and sum(edgeDistance) is always >= direct distance.
  auto heuristic = [&](std::size_t i) -> double {
    if (i == goal) return 0.0;
    const double d = (g.posLy[i] - g.posLy[goal]).length();
    return std::max(0.0, d);
  };

  struct QN {
    double f{0.0};
    double g{0.0};
    int hops{0};
    std::size_t i{0};
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      if (a.hops != b.hops) return a.hops > b.hops;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;

  gScore[start] = 0.0;
  open.push(QN{heuristic(start), 0.0, 0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (cur.i >= N) continue;
    if (closed[cur.i]) continue;

    closed[cur.i] = 1;
    ++expansions;

    if (cur.i == goal) {
      std::vector<sim::SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(g.ids[(std::size_t)at]);
      }
      std::reverse(path.begin(), path.end());

      out.reached = true;
      out.path = std::move(path);
      out.cost = gScore[goal];
      out.expansions = (int)expansions;
      return out;
    }

    const sim::SystemId curId = g.ids[cur.i];

    for (const Adj& e : g.adj[cur.i]) {
      const std::size_t j = e.to;
      if (j >= N) continue;
      if (closed[j]) continue;

      const sim::SystemId nid = g.ids[j];
      if (bannedNodes && bannedNodes->count(nid) != 0) continue;
      if (bannedEdges && bannedEdges->count(Edge{curId, nid}) != 0) continue;

      const double ec = edgeCostLy(e.distanceLy, e.bandwidth01, e.risk01, travel);
      const double tentative = gScore[cur.i] + ec;

      if (!std::isfinite(tentative)) continue;

      if (tentative + 1e-12 < gScore[j]) {
        gScore[j] = tentative;
        cameFrom[j] = (int)cur.i;
        const double f = tentative + heuristic(j);
        open.push(QN{f, tentative, cur.hops + 1, j});
      }
    }
  }

  out.expansions = (int)expansions;
  return out;
}

} // namespace

std::vector<HyperlaneKRoute> plotKHyperlaneRoutesAStarCost(const std::vector<sim::SystemStub>& nodes,
                                                          const HyperlaneNetwork& net,
                                                          sim::SystemId startId,
                                                          sim::SystemId goalId,
                                                          const HyperlaneTravelParams& travel,
                                                          std::size_t k,
                                                          std::size_t maxExpansionsPerSolve) {
  std::vector<HyperlaneKRoute> out;
  if (k == 0) return out;
  if (startId == 0 || goalId == 0) return out;
  if (nodes.empty()) return out;

  const Graph g = buildGraph(nodes, net);
  if (g.ids.empty()) return out;

  if (g.idx.find(startId) == g.idx.end() || g.idx.find(goalId) == g.idx.end()) return out;

  // First shortest path.
  {
    auto base = aStarSolve(g, startId, goalId, travel, nullptr, nullptr, maxExpansionsPerSolve);
    if (base.path.empty()) return out;

    HyperlanePathMetrics m = computeMetrics(g, base.path, travel);
    if (!m.reachable) return out;
    out.push_back(HyperlaneKRoute{std::move(base.path), m});
  }

  struct Candidate {
    std::vector<sim::SystemId> path;
    HyperlanePathMetrics metrics;
  };

  const auto containsPath = [](const std::vector<HyperlaneKRoute>& routes, const std::vector<sim::SystemId>& p) {
    for (const auto& r : routes) {
      if (r.path == p) return true;
    }
    return false;
  };

  const auto containsPathCand = [](const std::vector<Candidate>& cands, const std::vector<sim::SystemId>& p) {
    for (const auto& c : cands) {
      if (c.path == p) return true;
    }
    return false;
  };

  const auto betterCand = [](const Candidate& a, const Candidate& b) {
    const double da = a.metrics.costLy;
    const double db = b.metrics.costLy;
    if (std::abs(da - db) > 1e-12) return da < db;
    return pathLexLess(a.path, b.path);
  };

  std::vector<Candidate> candidates;

  // Yen's algorithm:
  //  - out[0] is the shortest.
  //  - candidates holds the next-best deviations.
  for (std::size_t kth = 1; kth < k; ++kth) {
    const auto& prev = out[kth - 1].path;
    if (prev.size() < 2) break;

    for (std::size_t i = 0; i + 1 < prev.size(); ++i) {
      const sim::SystemId spurNode = prev[i];

      // Root path includes spur node.
      std::vector<sim::SystemId> root(prev.begin(), prev.begin() + (long long)i + 1);

      // Ban nodes in the root path *except* the spur node to enforce loopless paths.
      std::unordered_set<sim::SystemId> bannedNodes;
      bannedNodes.reserve(i);
      for (std::size_t r = 0; r < i; ++r) bannedNodes.insert(root[r]);

      // Ban edges that would recreate any previously found shortest path with this root.
      std::unordered_set<Edge, EdgeHash> bannedEdges;
      for (const auto& found : out) {
        const auto& p = found.path;
        if (p.size() > i + 1 && pathHasPrefix(p, root)) {
          bannedEdges.insert(Edge{p[i], p[i + 1]});
        }
      }

      const auto spur = aStarSolve(g, spurNode, goalId, travel,
                                  &bannedNodes,
                                  &bannedEdges,
                                  maxExpansionsPerSolve);

      if (spur.path.empty()) continue;

      // Combine root + spur (skip spurNode duplicate).
      std::vector<sim::SystemId> total = root;
      if (spur.path.size() > 1) {
        total.insert(total.end(), spur.path.begin() + 1, spur.path.end());
      }

      if (total.size() < 2) continue;
      if (containsPath(out, total)) continue;
      if (containsPathCand(candidates, total)) continue;

      HyperlanePathMetrics m = computeMetrics(g, total, travel);
      if (!m.reachable) continue;

      candidates.push_back(Candidate{std::move(total), m});
    }

    if (candidates.empty()) break;

    // Pick the best candidate.
    std::size_t bestIdx = 0;
    for (std::size_t i = 1; i < candidates.size(); ++i) {
      if (betterCand(candidates[i], candidates[bestIdx])) bestIdx = i;
    }

    Candidate best = std::move(candidates[bestIdx]);
    candidates.erase(candidates.begin() + (long long)bestIdx);

    out.push_back(HyperlaneKRoute{std::move(best.path), best.metrics});
  }

  return out;
}

} // namespace stellar::proc
