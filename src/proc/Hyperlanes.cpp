#include "stellar/proc/Hyperlanes.h"

#include "stellar/core/Hash.h"
#include "stellar/proc/GalaxyHazards.h"
#include "stellar/proc/GalaxyRegions.h"
#include "stellar/proc/TradeProfile.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>

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

static double unitFromHash(core::u64 h) {
  // Map a 64-bit hash to [0,1) using the high 53 bits (double mantissa).
  // This is stable across platforms as long as IEEE-754 doubles are used.
  const core::u64 v = (h >> 11) & ((1ull << 53) - 1ull);
  return static_cast<double>(v) / static_cast<double>(1ull << 53);
}

static double regionNodeMul(GalaxyRegionKind k) {
  // Mild art-direction knob: different region kinds bias lane density/importance.
  switch (k) {
    case GalaxyRegionKind::Core: return 1.10;
    case GalaxyRegionKind::InnerDisc: return 1.02;
    case GalaxyRegionKind::OuterRim: return 0.95;
    case GalaxyRegionKind::Nebula: return 1.00;
    case GalaxyRegionKind::Cluster: return 1.12;
    case GalaxyRegionKind::Rift: return 0.78;
    default: return 1.00;
  }
}

static double regionRiskAdd(GalaxyRegionKind k) {
  switch (k) {
    case GalaxyRegionKind::Core: return -0.05;
    case GalaxyRegionKind::InnerDisc: return -0.02;
    case GalaxyRegionKind::OuterRim: return 0.06;
    case GalaxyRegionKind::Nebula: return 0.10;
    case GalaxyRegionKind::Cluster: return 0.02;
    case GalaxyRegionKind::Rift: return 0.22;
    default: return 0.0;
  }
}

struct CellCoord {
  int x{0}, y{0}, z{0};
  bool operator==(const CellCoord& o) const { return x == o.x && y == o.y && z == o.z; }
};

struct CellCoordHash {
  std::size_t operator()(const CellCoord& c) const noexcept {
    core::u64 h = 0xcbf29ce484222325ull;
    h = core::hashCombine(h, (core::u64)(std::uint32_t)c.x);
    h = core::hashCombine(h, (core::u64)(std::uint32_t)c.y);
    h = core::hashCombine(h, (core::u64)(std::uint32_t)c.z);
    return (std::size_t)h;
  }
};

struct Node {
  sim::SystemId id{0};
  math::Vec3d posLy{0, 0, 0};
  core::u32 factionId{0};

  // Procedural knobs.
  double hub{0.0};
  double pop{0.0};
  double wealth{0.0};
  double lawlessness{0.0};

  GalaxyRegionKind regionKind{GalaxyRegionKind::InnerDisc};
  double regionEdge01{0.0};

  // Derived.
  double nodeMass{0.2};
  double nodeMul{1.0};
};

static Node makeNode(core::u64 universeSeed, const sim::SystemStub& s, double regionCellSizeLy) {
  Node n;
  n.id = s.id;
  n.posLy = s.posLy;
  n.factionId = s.factionId;

  const TradeProfile p = generateTradeProfile(universeSeed, s);
  n.hub = clamp01(p.hub);
  n.pop = clamp01(p.population);
  n.wealth = clamp01(p.wealth);
  n.lawlessness = clamp01(p.lawlessness);

  if (regionCellSizeLy > 0.0) {
    const auto reg = sampleGalaxyRegion(universeSeed, s.posLy, regionCellSizeLy);
    n.regionKind = reg.kind;
    n.regionEdge01 = clamp01(reg.edge01);
    n.nodeMul = regionNodeMul(reg.kind);
  }

  // A tiny "importance" scalar used to bias lane costs.
  // - hub/pop drive connectivity
  // - wealth matters, but less
  // - keep a non-zero floor so frontier systems still connect
  const double hp = std::sqrt(std::max(0.0, n.hub * (0.25 + 0.75 * n.pop)));
  const double w = n.wealth;
  n.nodeMass = std::clamp((0.18 + 0.62 * hp + 0.20 * w) * n.nodeMul, 0.05, 1.35);

  if (!std::isfinite(n.nodeMass)) n.nodeMass = 0.18;
  return n;
}

struct EdgeKey {
  sim::SystemId a{0};
  sim::SystemId b{0};
  bool operator==(const EdgeKey& o) const { return a == o.a && b == o.b; }
};

struct EdgeKeyHash {
  std::size_t operator()(const EdgeKey& k) const noexcept {
    core::u64 h = 0xcbf29ce484222325ull;
    h = core::hashCombine(h, k.a);
    h = core::hashCombine(h, k.b);
    return (std::size_t)h;
  }
};

struct EdgeCand {
  int ia{0};
  int ib{0};
  sim::SystemId a{0};
  sim::SystemId b{0};
  double d{0.0};

  // Lower is better (used by MST).
  double cost{0.0};

  // Output metrics.
  double bandwidth01{0.0};
  double risk01{0.0};
};

struct Dsu {
  std::vector<int> p;
  std::vector<int> r;

  explicit Dsu(int n = 0) { reset(n); }

  void reset(int n) {
    p.resize(n);
    r.assign(n, 0);
    for (int i = 0; i < n; ++i) p[i] = i;
  }

  int find(int x) {
    while (p[x] != x) {
      p[x] = p[p[x]];
      x = p[x];
    }
    return x;
  }

  bool unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a == b) return false;
    if (r[a] < r[b]) std::swap(a, b);
    p[b] = a;
    if (r[a] == r[b]) r[a]++;
    return true;
  }
};

static double edgeBandwidth01(const Node& a, const Node& b) {
  // Macro capacity: high for hub-to-hub corridors, low for frontier links.
  const double hp = std::sqrt(std::max(0.0, (a.hub * a.pop) * (b.hub * b.pop)));
  double bw = 0.10 + 0.90 * hp;
  // Same-faction lanes are typically more maintained.
  if (a.factionId != 0 && a.factionId == b.factionId) bw *= 1.05;
  // Cluster regions tend to be more connected.
  bw *= 0.98 + 0.04 * (regionNodeMul(a.regionKind) + regionNodeMul(b.regionKind));
  return clamp01(bw);
}

static double edgeRisk01(const Node& a, const Node& b, double bandwidth01) {
  // Base risk: average lawlessness.
  double r = 0.18 + 0.70 * (0.5 * (a.lawlessness + b.lawlessness));

  // Region flavor.
  r += 0.5 * (regionRiskAdd(a.regionKind) + regionRiskAdd(b.regionKind));

  // Boundaries feel "contested".
  r += 0.22 * (0.5 * (a.regionEdge01 + b.regionEdge01));

  // High-bandwidth lanes are more patrolled.
  r *= (1.0 - 0.35 * bandwidth01);

  return clamp01(r);
}

static double edgeCost(const Node& a, const Node& b, double dLy) {
  // MST cost: distance divided by importance.
  const double mass = std::sqrt(std::max(0.0, a.nodeMass * b.nodeMass));
  const double denom = 0.22 + 0.78 * mass;

  // Edges near region boundaries are slightly "harder".
  const double edgePenalty = 1.0 + 0.35 * (0.5 * (a.regionEdge01 + b.regionEdge01));

  // Rift regions are sparse; penalize.
  const double riftPenalty = (a.regionKind == GalaxyRegionKind::Rift || b.regionKind == GalaxyRegionKind::Rift) ? 1.35 : 1.0;

  const double c = (dLy / std::max(1e-6, denom)) * edgePenalty * riftPenalty;
  if (!std::isfinite(c)) return dLy;
  return std::max(0.0, c);
}

static void applyHazardsToEdge(core::u64 universeSeed,
                               const HyperlaneParams& hp,
                               const Node& a,
                               const Node& b,
                               EdgeCand& e) {
  if (!hp.hazards.enabled) return;

  GalaxyHazardsParams haz{};
  haz.timeDays = hp.hazards.timeDays;
  haz.driftLyPerDay = hp.hazards.driftLyPerDay;
  haz.nebulaScaleLy = hp.hazards.nebulaScaleLy;
  haz.stormScaleLy = hp.hazards.stormScaleLy;
  // Reuse the same region scale as hyperlane generation when enabled.
  haz.regionCellSizeLy = (hp.regionCellSizeLy > 0.0) ? hp.regionCellSizeLy : 0.0;

  const double h = sampleGalaxyHazardAvgOnSegment(universeSeed, a.posLy, b.posLy, haz, hp.hazards.edgeSamples);

  // Topology bias.
  e.cost *= (1.0 + hp.hazards.hazardToCost * h);

  // Output metrics.
  e.risk01 = clamp01(e.risk01 + hp.hazards.hazardToRisk * h);
  e.bandwidth01 = clamp01(e.bandwidth01 * (1.0 - hp.hazards.hazardToBandwidth * h));
}

static CellCoord cellOf(const math::Vec3d& p, double cellSize) {
  const double s = std::max(1e-6, cellSize);
  CellCoord c;
  c.x = (int)std::floor(p.x / s);
  c.y = (int)std::floor(p.y / s);
  c.z = (int)std::floor(p.z / s);
  return c;
}

HyperlaneNetwork generateHyperlaneNetwork(core::u64 universeSeed,
                                          const std::vector<sim::SystemStub>& stubs,
                                          const HyperlaneParams& params) {
  HyperlaneNetwork out;

  if (stubs.size() < 2) return out;

  HyperlaneParams p = params;
  p.maxNeighborDistanceLy = std::max(1.0, p.maxNeighborDistanceLy);
  p.neighborK = std::clamp(p.neighborK, 1, 16);
  p.minDegree = std::clamp(p.minDegree, 0, 8);
  p.extraEdgeChance = std::clamp(p.extraEdgeChance, 0.0, 1.0);

  if (p.hazards.enabled) {
    p.hazards.edgeSamples = std::clamp(p.hazards.edgeSamples, 1, 11);
    p.hazards.driftLyPerDay = std::max(0.0, p.hazards.driftLyPerDay);
    p.hazards.nebulaScaleLy = std::max(1.0, p.hazards.nebulaScaleLy);
    p.hazards.stormScaleLy = std::max(1.0, p.hazards.stormScaleLy);
    p.hazards.hazardToCost = std::max(0.0, p.hazards.hazardToCost);
    p.hazards.hazardToRisk = std::max(0.0, p.hazards.hazardToRisk);
    p.hazards.hazardToBandwidth = std::clamp(p.hazards.hazardToBandwidth, 0.0, 0.95);
  }

  // Precompute node metadata.
  std::vector<Node> nodes;
  nodes.reserve(stubs.size());
  for (const auto& s : stubs) {
    if (s.id == 0) continue;
    nodes.push_back(makeNode(universeSeed, s, p.regionCellSizeLy));
  }

  const int n = (int)nodes.size();
  if (n < 2) return out;

  // Build a spatial hash grid (cell size = maxNeighborDistanceLy).
  std::unordered_map<CellCoord, std::vector<int>, CellCoordHash> grid;
  grid.reserve((std::size_t)n * 2);
  for (int i = 0; i < n; ++i) {
    const CellCoord c = cellOf(nodes[i].posLy, p.maxNeighborDistanceLy);
    grid[c].push_back(i);
  }

  std::unordered_map<EdgeKey, int, EdgeKeyHash> edgeIndex;
  edgeIndex.reserve((std::size_t)n * (std::size_t)p.neighborK);

  std::vector<EdgeCand> cands;
  cands.reserve((std::size_t)n * (std::size_t)p.neighborK);

  // Candidate edges: each node proposes up to K nearest neighbors.
  for (int i = 0; i < n; ++i) {
    const CellCoord c = cellOf(nodes[i].posLy, p.maxNeighborDistanceLy);

    std::vector<std::pair<double, int>> neigh;
    neigh.reserve(64);

    for (int dz = -1; dz <= 1; ++dz) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;
          for (int j : it->second) {
            if (j == i) continue;
            const double d = distLy(nodes[i].posLy, nodes[j].posLy);
            if (d <= p.maxNeighborDistanceLy + 1e-9) {
              neigh.emplace_back(d, j);
            }
          }
        }
      }
    }

    std::sort(neigh.begin(), neigh.end(), [&](const auto& a, const auto& b) {
      if (a.first != b.first) return a.first < b.first;
      return nodes[a.second].id < nodes[b.second].id;
    });

    if ((int)neigh.size() > p.neighborK) neigh.resize((std::size_t)p.neighborK);

    for (const auto& [d, j] : neigh) {
      const sim::SystemId ida = nodes[i].id;
      const sim::SystemId idb = nodes[j].id;
      const sim::SystemId aId = std::min(ida, idb);
      const sim::SystemId bId = std::max(ida, idb);
      const int ia = (ida <= idb) ? i : j;
      const int ib = (ida <= idb) ? j : i;

      EdgeKey key{aId, bId};
      if (edgeIndex.find(key) != edgeIndex.end()) continue;

      const Node& na = nodes[ia];
      const Node& nb = nodes[ib];

      EdgeCand e;
      e.ia = ia;
      e.ib = ib;
      e.a = aId;
      e.b = bId;
      e.d = d;
      e.bandwidth01 = edgeBandwidth01(na, nb);
      e.risk01 = edgeRisk01(na, nb, e.bandwidth01);
      e.cost = edgeCost(na, nb, d);

      applyHazardsToEdge(universeSeed, p, na, nb, e);

      const int idx = (int)cands.size();
      cands.push_back(e);
      edgeIndex.emplace(key, idx);
    }
  }

  if (cands.empty()) return out;

  // Sort candidates for MST selection.
  std::vector<int> candOrder(cands.size());
  for (std::size_t i = 0; i < candOrder.size(); ++i) candOrder[i] = (int)i;

  std::sort(candOrder.begin(), candOrder.end(), [&](int ia, int ib) {
    const auto& a = cands[(std::size_t)ia];
    const auto& b = cands[(std::size_t)ib];
    if (a.cost != b.cost) return a.cost < b.cost;
    if (a.d != b.d) return a.d < b.d;
    if (a.a != b.a) return a.a < b.a;
    return a.b < b.b;
  });

  std::vector<char> selected(cands.size(), 0);
  std::vector<int> deg(n, 0);

  // Backbone MST.
  {
    Dsu dsu(n);
    int edgesAdded = 0;
    for (int idx : candOrder) {
      const auto& e = cands[(std::size_t)idx];
      if (dsu.unite(e.ia, e.ib)) {
        selected[(std::size_t)idx] = 1;
        deg[e.ia]++;
        deg[e.ib]++;
        edgesAdded++;
        if (edgesAdded >= n - 1) break;
      }
    }

    // Best-effort: if graph is disconnected and forceConnected is enabled,
    // connect components by adding long "bridge" edges between closest pairs.
    //
    // IMPORTANT: this must be deterministic; avoid unordered_map iteration.
    if (p.forceConnected) {
      std::vector<int> root(n);
      for (int i = 0; i < n; ++i) root[i] = dsu.find(i);

      std::vector<int> roots = root;
      std::sort(roots.begin(), roots.end());
      roots.erase(std::unique(roots.begin(), roots.end()), roots.end());

      if (roots.size() > 1) {
        struct Comp {
          sim::SystemId minId{0};
          std::vector<int> members;
        };

        std::vector<Comp> comps;
        comps.reserve(roots.size());

        for (int r : roots) {
          Comp c;
          for (int i = 0; i < n; ++i) {
            if (root[i] == r) c.members.push_back(i);
          }
          std::sort(c.members.begin(), c.members.end(), [&](int a, int b) { return nodes[a].id < nodes[b].id; });
          c.minId = nodes[c.members.front()].id;
          comps.push_back(std::move(c));
        }

        std::sort(comps.begin(), comps.end(), [](const Comp& a, const Comp& b) { return a.minId < b.minId; });

        const double eps = 1e-9;
        while (comps.size() > 1) {
          double bestD = 1e300;
          sim::SystemId bestAId = 0;
          sim::SystemId bestBId = 0;
          int bestIa = -1;
          int bestIb = -1;
          std::size_t bestCi = 0;
          std::size_t bestCj = 1;

          for (std::size_t ci = 0; ci < comps.size(); ++ci) {
            for (std::size_t cj = ci + 1; cj < comps.size(); ++cj) {
              for (int ia2 : comps[ci].members) {
                for (int ib2 : comps[cj].members) {
                  const double d = distLy(nodes[ia2].posLy, nodes[ib2].posLy);
                  const sim::SystemId ida = nodes[ia2].id;
                  const sim::SystemId idb = nodes[ib2].id;
                  const sim::SystemId aId = std::min(ida, idb);
                  const sim::SystemId bId = std::max(ida, idb);

                  bool better = false;
                  if (d + eps < bestD) {
                    better = true;
                  } else if (std::abs(d - bestD) <= eps) {
                    if (aId < bestAId || (aId == bestAId && bId < bestBId)) better = true;
                  }

                  if (better) {
                    bestD = d;
                    bestAId = aId;
                    bestBId = bId;
                    bestIa = (ida <= idb) ? ia2 : ib2;
                    bestIb = (ida <= idb) ? ib2 : ia2;
                    bestCi = ci;
                    bestCj = cj;
                  }
                }
              }
            }
          }

          if (bestIa < 0 || bestIb < 0 || !std::isfinite(bestD)) break;

          EdgeKey key{bestAId, bestBId};
          int idx = -1;
          auto itIdx = edgeIndex.find(key);
          if (itIdx != edgeIndex.end()) {
            idx = itIdx->second;
          } else {
            EdgeCand e;
            e.ia = bestIa;
            e.ib = bestIb;
            e.a = bestAId;
            e.b = bestBId;
            e.d = bestD;
            e.bandwidth01 = edgeBandwidth01(nodes[bestIa], nodes[bestIb]);
            e.risk01 = edgeRisk01(nodes[bestIa], nodes[bestIb], e.bandwidth01);
            e.cost = edgeCost(nodes[bestIa], nodes[bestIb], bestD) * 1.10; // slight penalty for long bridges

            applyHazardsToEdge(universeSeed, p, nodes[bestIa], nodes[bestIb], e);

            idx = (int)cands.size();
            cands.push_back(e);
            selected.push_back(0);
            edgeIndex.emplace(key, idx);
          }

          if (idx >= 0 && !selected[(std::size_t)idx]) {
            selected[(std::size_t)idx] = 1;
            deg[cands[(std::size_t)idx].ia]++;
            deg[cands[(std::size_t)idx].ib]++;
          }

          // Merge components bestCi and bestCj into bestCi, then re-sort.
          auto& aMem = comps[bestCi].members;
          auto& bMem = comps[bestCj].members;
          aMem.insert(aMem.end(), bMem.begin(), bMem.end());
          std::sort(aMem.begin(), aMem.end(), [&](int a, int b) { return nodes[a].id < nodes[b].id; });
          comps[bestCi].minId = nodes[aMem.front()].id;
          comps.erase(comps.begin() + (std::ptrdiff_t)bestCj);

          std::sort(comps.begin(), comps.end(), [](const Comp& a, const Comp& b) { return a.minId < b.minId; });
        }
      }
    }
  }

  // Build adjacency candidate lists for min-degree augmentation.
  if (p.minDegree > 0) {
    std::vector<std::vector<int>> inc((std::size_t)n);
    inc.reserve((std::size_t)n);

    for (std::size_t ei = 0; ei < cands.size(); ++ei) {
      const auto& e = cands[ei];
      inc[(std::size_t)e.ia].push_back((int)ei);
      inc[(std::size_t)e.ib].push_back((int)ei);
    }

    for (int i = 0; i < n; ++i) {
      auto& v = inc[(std::size_t)i];
      std::sort(v.begin(), v.end(), [&](int ea, int eb) {
        const auto& a = cands[(std::size_t)ea];
        const auto& b = cands[(std::size_t)eb];
        if (a.cost != b.cost) return a.cost < b.cost;
        // Prefer higher bandwidth for the same cost.
        if (a.bandwidth01 != b.bandwidth01) return a.bandwidth01 > b.bandwidth01;
        if (a.a != b.a) return a.a < b.a;
        return a.b < b.b;
      });

      while (deg[i] < p.minDegree) {
        bool added = false;
        for (int ei : v) {
          if (selected[(std::size_t)ei]) continue;
          const auto& e = cands[(std::size_t)ei];
          selected[(std::size_t)ei] = 1;
          deg[e.ia]++;
          deg[e.ib]++;
          added = true;
          break;
        }
        if (!added) break;
      }
    }
  }

  // Probabilistic extras.
  if (p.extraEdgeChance > 0.0) {
    const core::u64 salt = core::fnv1a64("proc_hyperlane_extra");
    const core::u64 base = core::hashCombine(universeSeed, salt);

    for (std::size_t ei = 0; ei < cands.size(); ++ei) {
      if (selected[ei]) continue;
      const auto& e = cands[ei];

      // Distance-suppressed probability; bandwidth boosts, risk suppresses.
      const double d01 = std::clamp(e.d / std::max(1.0, p.maxNeighborDistanceLy), 0.0, 3.0);
      const double distW = std::exp(-1.35 * d01);
      const double bwW = 0.35 + 0.65 * e.bandwidth01;
      const double riskW = 1.0 - 0.55 * e.risk01;

      double prob = p.extraEdgeChance * distW * bwW * riskW;
      prob = std::clamp(prob, 0.0, 1.0);

      const core::u64 h = core::hashCombine(base, core::hashCombine(e.a, e.b));
      const double u = unitFromHash(h);
      if (u < prob) {
        selected[ei] = 1;
        deg[e.ia]++;
        deg[e.ib]++;
      }
    }
  }

  // Emit final edge list.
  out.edges.reserve(cands.size());
  for (std::size_t ei = 0; ei < cands.size(); ++ei) {
    if (!selected[ei]) continue;
    const auto& e = cands[ei];

    HyperlaneEdge he;
    he.a = e.a;
    he.b = e.b;
    he.distanceLy = e.d;
    he.bandwidth01 = clamp01(e.bandwidth01);
    he.risk01 = clamp01(e.risk01);

    if (he.a != 0 && he.b != 0 && he.a < he.b) {
      out.edges.push_back(std::move(he));
    }
  }

  // Deterministic output ordering.
  std::sort(out.edges.begin(), out.edges.end(), [](const HyperlaneEdge& a, const HyperlaneEdge& b) {
    if (a.bandwidth01 != b.bandwidth01) return a.bandwidth01 > b.bandwidth01;
    if (a.distanceLy != b.distanceLy) return a.distanceLy < b.distanceLy;
    if (a.a != b.a) return a.a < b.a;
    return a.b < b.b;
  });

  if (p.maxEdges > 0 && out.edges.size() > p.maxEdges) {
    out.edges.resize(p.maxEdges);
  }

  return out;
}

} // namespace stellar::proc
