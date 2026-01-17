#include "stellar/proc/HyperlaneCentrality.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <vector>

namespace stellar::proc {

namespace {

static double clamp01(double x) {
  return std::clamp(x, 0.0, 1.0);
}

// Matches the travel metric used by HyperlaneRouter / HyperlaneKRoutes.
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

static bool distLess(double a, double b, double relEps) {
  const double eps = std::max(0.0, relEps);
  const double s = std::max(1.0, std::max(std::abs(a), std::abs(b)));
  return a < b - eps * s;
}

static bool distEq(double a, double b, double relEps) {
  const double eps = std::max(0.0, relEps);
  const double s = std::max(1.0, std::max(std::abs(a), std::abs(b)));
  return std::abs(a - b) <= eps * s;
}

struct Pred {
  int v{0};
  std::size_t edgeIdx{0};
};

struct Adj {
  int to{0};
  double w{0.0};
  std::size_t edgeIdx{0};
};

struct Graph {
  std::vector<sim::SystemId> ids;
  std::unordered_map<sim::SystemId, int> idx;
  std::vector<std::vector<Adj>> adj;
};

static Graph buildGraph(const std::vector<sim::SystemStub>& nodes,
                        const HyperlaneNetwork& net,
                        const HyperlaneTravelParams& travel) {
  Graph g{};

  g.ids.reserve(nodes.size() + net.edges.size() * 2);
  for (const auto& s : nodes) {
    if (s.id != 0) g.ids.push_back(s.id);
  }
  for (const auto& e : net.edges) {
    if (e.a != 0) g.ids.push_back(e.a);
    if (e.b != 0) g.ids.push_back(e.b);
  }

  std::sort(g.ids.begin(), g.ids.end());
  g.ids.erase(std::unique(g.ids.begin(), g.ids.end()), g.ids.end());

  g.idx.reserve(g.ids.size() * 2);
  for (int i = 0; i < (int)g.ids.size(); ++i) {
    g.idx[g.ids[(std::size_t)i]] = i;
  }

  g.adj.assign(g.ids.size(), {});

  for (std::size_t ei = 0; ei < net.edges.size(); ++ei) {
    const auto& e = net.edges[ei];
    auto ita = g.idx.find(e.a);
    auto itb = g.idx.find(e.b);
    if (ita == g.idx.end() || itb == g.idx.end()) continue;
    const int ia = ita->second;
    const int ib = itb->second;
    if (ia == ib) continue;

    const double w = edgeCostLy(e.distanceLy, e.bandwidth01, e.risk01, travel);

    g.adj[(std::size_t)ia].push_back(Adj{ib, w, ei});
    g.adj[(std::size_t)ib].push_back(Adj{ia, w, ei});
  }

  // Deterministic neighbor iteration (important for repeatable floating tie cases).
  for (auto& v : g.adj) {
    std::sort(v.begin(), v.end(), [&](const Adj& a, const Adj& b) {
      if (a.to != b.to) return g.ids[(std::size_t)a.to] < g.ids[(std::size_t)b.to];
      if (a.w != b.w) return a.w < b.w;
      return a.edgeIdx < b.edgeIdx;
    });
  }

  return g;
}

static std::vector<int> chooseSourceIndices(const Graph& g, std::size_t k, core::u64 seed) {
  std::vector<int> out;
  const std::size_t n = g.ids.size();
  if (n == 0) return out;

  if (k == 0 || k >= n) {
    out.reserve(n);
    for (int i = 0; i < (int)n; ++i) out.push_back(i);
    return out;
  }

  struct Cand {
    core::u64 key{0};
    int idx{0};
  };

  struct CandLess {
    // Max-heap: top() is the largest key (worst of bottom-K).
    bool operator()(const Cand& a, const Cand& b) const { return a.key < b.key; }
  };

  std::priority_queue<Cand, std::vector<Cand>, CandLess> heap;
  heap = {};

  const core::u64 base = core::hashCombine(core::fnv1a64("hyperlane_betweenness_sources/v1"), seed);

  for (int i = 0; i < (int)n; ++i) {
    const core::u64 key = core::hashCombine(base, (core::u64)g.ids[(std::size_t)i]);
    Cand c{key, i};
    if (heap.size() < k) {
      heap.push(c);
    } else if (c.key < heap.top().key) {
      heap.pop();
      heap.push(c);
    }
  }

  out.reserve(heap.size());
  while (!heap.empty()) {
    out.push_back(heap.top().idx);
    heap.pop();
  }

  // Stable ordering improves determinism when users compare visualizations.
  std::sort(out.begin(), out.end(), [&](int a, int b) {
    const core::u64 ka = core::hashCombine(base, (core::u64)g.ids[(std::size_t)a]);
    const core::u64 kb = core::hashCombine(base, (core::u64)g.ids[(std::size_t)b]);
    if (ka != kb) return ka < kb;
    return g.ids[(std::size_t)a] < g.ids[(std::size_t)b];
  });

  return out;
}

} // namespace

HyperlaneBetweennessResult estimateHyperlaneBetweennessCentrality(const std::vector<sim::SystemStub>& nodes,
                                                                  const HyperlaneNetwork& net,
                                                                  const HyperlaneBetweennessParams& params) {
  HyperlaneBetweennessResult out{};

  // Always size edge output to match `net.edges` order.
  out.edgeBetweenness.assign(net.edges.size(), 0.0);

  const Graph g = buildGraph(nodes, net, params.travel);
  out.nodeIds = g.ids;
  out.nodeBetweenness.assign(g.ids.size(), 0.0);

  const int N = (int)g.ids.size();
  if (N <= 0 || net.edges.empty()) {
    // For empty edge set, betweenness is trivially 0.
    out.maxEdge = 0.0;
    out.maxNode = 0.0;
    out.sourcesUsed = 0;
    return out;
  }

  const std::vector<int> sources = chooseSourceIndices(g, params.sampleSources, params.sampleSeed);
  out.sourcesUsed = sources.size();
  if (sources.empty()) return out;

  std::vector<double> dist((std::size_t)N, std::numeric_limits<double>::infinity());
  std::vector<double> sigma((std::size_t)N, 0.0);
  std::vector<double> delta((std::size_t)N, 0.0);
  std::vector<char> visited((std::size_t)N, 0);

  std::vector<std::vector<Pred>> pred((std::size_t)N);
  std::vector<int> stack;
  stack.reserve((std::size_t)N);

  struct Q {
    double d{0.0};
    int v{0};
  };

  struct Cmp {
    bool operator()(const Q& a, const Q& b) const {
      if (a.d != b.d) return a.d > b.d;
      return a.v > b.v;
    }
  };

  const double tieEps = params.tieEps;

  for (int sIdx : sources) {
    if (sIdx < 0 || sIdx >= N) continue;

    std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::infinity());
    std::fill(sigma.begin(), sigma.end(), 0.0);
    std::fill(delta.begin(), delta.end(), 0.0);
    std::fill(visited.begin(), visited.end(), 0);
    stack.clear();

    dist[(std::size_t)sIdx] = 0.0;
    sigma[(std::size_t)sIdx] = 1.0;
    pred[(std::size_t)sIdx].clear();

    std::priority_queue<Q, std::vector<Q>, Cmp> pq;
    pq.push(Q{0.0, sIdx});

    std::size_t expansions = 0;

    while (!pq.empty()) {
      const Q cur = pq.top();
      pq.pop();

      const int v = cur.v;
      if (v < 0 || v >= N) continue;
      if (visited[(std::size_t)v]) continue;

      // Stale entry?
      if (distLess(dist[(std::size_t)v], cur.d, tieEps) || distLess(cur.d, dist[(std::size_t)v], tieEps)) {
        // If cur.d != dist[v] within eps, ignore.
        if (!distEq(cur.d, dist[(std::size_t)v], tieEps)) continue;
      }

      visited[(std::size_t)v] = 1;
      stack.push_back(v);

      ++expansions;
      if (params.maxExpansionsPerSource > 0 && expansions > params.maxExpansionsPerSource) {
        // Guard: stop early to keep UI responsive.
        break;
      }

      const double dv = dist[(std::size_t)v];
      if (!std::isfinite(dv)) continue;

      for (const Adj& e : g.adj[(std::size_t)v]) {
        const int w = e.to;
        if (w < 0 || w >= N) continue;
        const double alt = dv + std::max(0.0, e.w);
        if (!std::isfinite(alt)) continue;

        if (distLess(alt, dist[(std::size_t)w], tieEps)) {
          dist[(std::size_t)w] = alt;
          pq.push(Q{alt, w});

          sigma[(std::size_t)w] = sigma[(std::size_t)v];
          pred[(std::size_t)w].clear();
          pred[(std::size_t)w].push_back(Pred{v, e.edgeIdx});
        } else if (distEq(alt, dist[(std::size_t)w], tieEps)) {
          sigma[(std::size_t)w] += sigma[(std::size_t)v];
          pred[(std::size_t)w].push_back(Pred{v, e.edgeIdx});
        }
      }
    }

    // Accumulation: process vertices in reverse order of distance.
    for (int si = (int)stack.size() - 1; si >= 0; --si) {
      const int w = stack[(std::size_t)si];
      const double sw = sigma[(std::size_t)w];
      if (sw <= 0.0) continue;

      for (const Pred& p : pred[(std::size_t)w]) {
        const int v = p.v;
        if (v < 0 || v >= N) continue;
        const double sv = sigma[(std::size_t)v];
        if (sv <= 0.0) continue;

        // Dependency contribution.
        const double c = (sv / sw) * (1.0 + delta[(std::size_t)w]);
        if (!std::isfinite(c)) continue;

        delta[(std::size_t)v] += c;

        if (p.edgeIdx < out.edgeBetweenness.size()) {
          out.edgeBetweenness[p.edgeIdx] += c;
        }
      }

      if (w != sIdx) {
        out.nodeBetweenness[(std::size_t)w] += delta[(std::size_t)w];
      }
    }
  }

  // Scale sampled sources to approximate the full sum over all sources.
  // For undirected graphs, divide by 2 to avoid double counting.
  const double scale = ((double)N / (double)std::max<std::size_t>(1, sources.size())) * 0.5;

  for (double& v : out.edgeBetweenness) v *= scale;
  for (double& v : out.nodeBetweenness) v *= scale;

  out.maxEdge = 0.0;
  for (double v : out.edgeBetweenness) out.maxEdge = std::max(out.maxEdge, v);

  out.maxNode = 0.0;
  for (double v : out.nodeBetweenness) out.maxNode = std::max(out.maxNode, v);

  return out;
}

} // namespace stellar::proc
