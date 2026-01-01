#include "stellar/sim/NavRoute.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace stellar::sim {

namespace {

static double systemDistanceLy(const SystemStub& a, const SystemStub& b) {
  return (a.posLy - b.posLy).length();
}

struct CellCoord {
  long long x{0};
  long long y{0};
  long long z{0};

  bool operator==(const CellCoord& o) const { return x == o.x && y == o.y && z == o.z; }
};

static std::size_t hashCombine(std::size_t h, std::size_t v) {
  // A small, decent hash combine (boost-style).
  return h ^ (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
}

struct CellHash {
  std::size_t operator()(const CellCoord& c) const {
    std::size_t h = 0;
    h = hashCombine(h, std::hash<long long>()(c.x));
    h = hashCombine(h, std::hash<long long>()(c.y));
    h = hashCombine(h, std::hash<long long>()(c.z));
    return h;
  }
};

static CellCoord cellFor(const math::Vec3d& p, double cellSize) {
  // cellSize is assumed > 0.
  return CellCoord{
    (long long)std::floor(p.x / cellSize),
    (long long)std::floor(p.y / cellSize),
    (long long)std::floor(p.z / cellSize),
  };
}

static std::unordered_map<SystemId, std::size_t> buildIndex(const std::vector<SystemStub>& nodes) {
  std::unordered_map<SystemId, std::size_t> idx;
  idx.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const auto id = nodes[i].id;
    if (id != 0) idx[id] = i;
  }
  return idx;
}

static void setStats(RoutePlanStats* out, const RoutePlanStats& s) {
  if (out) *out = s;
}

static std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash>
buildGrid(const std::vector<SystemStub>& nodes, double cellSize) {
  std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash> grid;
  grid.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    grid[cellFor(nodes[i].posLy, cellSize)].push_back(i);
  }
  return grid;
}

struct Edge {
  SystemId from{0};
  SystemId to{0};
  bool operator==(const Edge& o) const { return from == o.from && to == o.to; }
};

struct EdgeHash {
  std::size_t operator()(const Edge& e) const {
    std::size_t h = 0;
    h = hashCombine(h, std::hash<SystemId>()(e.from));
    h = hashCombine(h, std::hash<SystemId>()(e.to));
    return h;
  }
};

static bool pathHasPrefix(const std::vector<SystemId>& path, const std::vector<SystemId>& prefix) {
  if (prefix.size() > path.size()) return false;
  return std::equal(prefix.begin(), prefix.end(), path.begin());
}

static bool pathLexLess(const std::vector<SystemId>& a, const std::vector<SystemId>& b) {
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

static double computeDistanceLy(const std::unordered_map<SystemId, std::size_t>& idx,
                               const std::vector<SystemStub>& nodes,
                               const std::vector<SystemId>& route) {
  if (route.size() < 2) return 0.0;

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    sum += systemDistanceLy(nodes[itA->second], nodes[itB->second]);
  }
  return sum;
}

static double computeCost(const std::unordered_map<SystemId, std::size_t>& idx,
                          const std::vector<SystemStub>& nodes,
                          const std::vector<SystemId>& route,
                          double costPerJump,
                          double costPerLy) {
  if (route.size() < 2) return 0.0;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    const double d = systemDistanceLy(nodes[itA->second], nodes[itB->second]);
    sum += costPerJump + costPerLy * d;
  }
  return sum;
}

struct AStarSolveResult {
  std::vector<SystemId> path;
  RoutePlanStats stats;
};

static AStarSolveResult aStarCostSolve(const std::vector<SystemStub>& nodes,
                                      const std::unordered_map<SystemId, std::size_t>& idx,
                                      const std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash>& grid,
                                      SystemId startId,
                                      SystemId goalId,
                                      double maxJumpLy,
                                      double costPerJump,
                                      double costPerLy,
                                      const std::unordered_set<SystemId>* bannedNodes,
                                      const std::unordered_set<Edge, EdgeHash>* bannedEdges,
                                      std::size_t maxExpansions) {
  AStarSolveResult out{};
  RoutePlanStats stats{};

  if (startId == 0 || goalId == 0) {
    out.stats = stats;
    return out;
  }
  if (maxJumpLy <= 0.0) {
    out.stats = stats;
    return out;
  }
  if (nodes.empty()) {
    out.stats = stats;
    return out;
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  if (bannedNodes) {
    if (bannedNodes->count(startId) != 0) {
      out.stats = stats;
      return out;
    }
    if (bannedNodes->count(goalId) != 0) {
      out.stats = stats;
      return out;
    }
  }

  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) {
    out.stats = stats;
    return out;
  }

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  if (start == goal) {
    stats.reached = true;
    stats.hops = 0;
    stats.distanceLy = 0.0;
    stats.cost = 0.0;
    out.stats = stats;
    out.path = {startId};
    return out;
  }

  std::vector<int> cameFrom(N, -1);
  std::vector<double> gScore(N, std::numeric_limits<double>::infinity());
  std::vector<char> closed(N, 0);

  auto heuristic = [&](std::size_t i) -> double {
    if (i == goal) return 0.0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    const double minHops = std::ceil(d / maxJumpLy);
    return minHops * costPerJump + d * costPerLy;
  };

  struct QN {
    double f{0.0};
    double g{0.0};
    std::size_t i{0};
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0.0;
  open.push(QN{heuristic(start), 0.0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    ++expansions;
    ++stats.visited;
    stats.expansions = (int)expansions;

    if (cur.i == goal) {
      std::vector<SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());

      stats.reached = true;
      stats.hops = path.size() > 0 ? (int)path.size() - 1 : 0;
      stats.distanceLy = computeDistanceLy(idx, nodes, path);
      stats.cost = gScore[goal];

      out.path = std::move(path);
      out.stats = stats;
      return out;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, maxJumpLy);

    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;

          for (const std::size_t j : it->second) {
            if (j == cur.i) continue;
            if (closed[j]) continue;

            const SystemId nid = nodes[j].id;
            if (bannedNodes && bannedNodes->count(nid) != 0) continue;
            if (bannedEdges && bannedEdges->count(Edge{nodes[cur.i].id, nid}) != 0) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const double legCost = costPerJump + costPerLy * d;
            const double tentative = gScore[cur.i] + legCost;

            if (tentative + 1e-12 < gScore[j]) {
              gScore[j] = tentative;
              cameFrom[j] = (int)cur.i;
              const double f = tentative + heuristic(j);
              open.push(QN{f, tentative, j});
            }
          }
        }
      }
    }
  }

  out.stats = stats;
  return out;
}

} // namespace

std::vector<SystemId> plotRouteAStarHops(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        RoutePlanStats* outStats,
                                        std::size_t maxExpansions) {
  RoutePlanStats stats{};

  if (startId == 0 || goalId == 0) {
    setStats(outStats, stats);
    return {};
  }
  if (maxJumpLy <= 0.0) {
    setStats(outStats, stats);
    return {};
  }
  if (nodes.empty()) {
    setStats(outStats, stats);
    return {};
  }

  const auto idx = buildIndex(nodes);
  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) {
    setStats(outStats, stats);
    return {};
  }

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  if (start == goal) {
    stats.reached = true;
    stats.hops = 0;
    stats.distanceLy = 0.0;
    stats.cost = 0.0;
    setStats(outStats, stats);
    return {startId};
  }

  // Spatial hash grid for neighbor queries.
  // Cell size is maxJumpLy; any reachable neighbor must be in the same or adjacent cell.
  const auto grid = buildGrid(nodes, maxJumpLy);

  std::vector<int> cameFrom(N, -1);
  std::vector<int> gScore(N, std::numeric_limits<int>::max());
  std::vector<char> closed(N, 0);

  auto heuristic = [&](std::size_t i) -> int {
    if (i == goal) return 0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    return (int)std::ceil(d / maxJumpLy);
  };

  struct QN {
    int f{0};
    int g{0};
    std::size_t i{0};
  };

  // Prefer lower f, then lower g, then lower index (stable/deterministic tie-break).
  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0;
  open.push(QN{heuristic(start), 0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    ++expansions;
    ++stats.visited;
    stats.expansions = (int)expansions;

    if (cur.i == goal) {
      std::vector<SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());

      stats.reached = true;
      stats.hops = path.size() > 0 ? (int)path.size() - 1 : 0;
      stats.distanceLy = routeDistanceLy(nodes, path);
      stats.cost = (double)stats.hops;

      setStats(outStats, stats);
      return path;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, maxJumpLy);

    // Neighbors: any node within maxJumpLy (search the surrounding 3x3x3 cell block).
    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;

          for (const std::size_t j : it->second) {
            if (j == cur.i) continue;
            if (closed[j]) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const int tentative = gScore[cur.i] + 1;
            if (tentative < gScore[j]) {
              gScore[j] = tentative;
              cameFrom[j] = (int)cur.i;
              const int f = tentative + heuristic(j);
              open.push(QN{f, tentative, j});
            }
          }
        }
      }
    }
  }

  // Not found (or hit expansion cap).
  setStats(outStats, stats);
  return {};
}

std::vector<SystemId> plotRouteAStarCost(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        double costPerJump,
                                        double costPerLy,
                                        RoutePlanStats* outStats,
                                        std::size_t maxExpansions) {
  RoutePlanStats stats{};

  if (startId == 0 || goalId == 0) {
    setStats(outStats, stats);
    return {};
  }
  if (maxJumpLy <= 0.0) {
    setStats(outStats, stats);
    return {};
  }
  if (nodes.empty()) {
    setStats(outStats, stats);
    return {};
  }

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  // Degenerate: no optimization signal. Fall back to the classic hop planner.
  if (costPerJump <= 0.0 && costPerLy <= 0.0) {
    return plotRouteAStarHops(nodes, startId, goalId, maxJumpLy, outStats, maxExpansions);
  }

  const auto idx = buildIndex(nodes);
  auto itS = idx.find(startId);
  auto itG = idx.find(goalId);
  if (itS == idx.end() || itG == idx.end()) {
    setStats(outStats, stats);
    return {};
  }

  const std::size_t start = itS->second;
  const std::size_t goal  = itG->second;
  const std::size_t N = nodes.size();

  if (start == goal) {
    stats.reached = true;
    stats.hops = 0;
    stats.distanceLy = 0.0;
    stats.cost = 0.0;
    setStats(outStats, stats);
    return {startId};
  }

  // Spatial hash grid for neighbor queries.
  const auto grid = buildGrid(nodes, maxJumpLy);

  std::vector<int> cameFrom(N, -1);
  std::vector<double> gScore(N, std::numeric_limits<double>::infinity());
  std::vector<char> closed(N, 0);

  auto heuristic = [&](std::size_t i) -> double {
    if (i == goal) return 0.0;
    const double d = systemDistanceLy(nodes[i], nodes[goal]);
    const double minHops = std::ceil(d / maxJumpLy);
    return minHops * costPerJump + d * costPerLy;
  };

  struct QN {
    double f{0.0};
    double g{0.0};
    std::size_t i{0};
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      if (a.f != b.f) return a.f > b.f;
      if (a.g != b.g) return a.g > b.g;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;
  gScore[start] = 0.0;
  open.push(QN{heuristic(start), 0.0, start});

  std::size_t expansions = 0;

  while (!open.empty() && expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;

    ++expansions;
    ++stats.visited;
    stats.expansions = (int)expansions;

    if (cur.i == goal) {
      std::vector<SystemId> path;
      for (int at = (int)goal; at != -1; at = cameFrom[(std::size_t)at]) {
        path.push_back(nodes[(std::size_t)at].id);
      }
      std::reverse(path.begin(), path.end());

      stats.reached = true;
      stats.hops = path.size() > 0 ? (int)path.size() - 1 : 0;
      stats.distanceLy = routeDistanceLy(nodes, path);
      stats.cost = gScore[goal];

      setStats(outStats, stats);
      return path;
    }

    const CellCoord c = cellFor(nodes[cur.i].posLy, maxJumpLy);

    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const CellCoord cc{c.x + dx, c.y + dy, c.z + dz};
          auto it = grid.find(cc);
          if (it == grid.end()) continue;

          for (const std::size_t j : it->second) {
            if (j == cur.i) continue;
            if (closed[j]) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const double legCost = costPerJump + costPerLy * d;
            const double tentative = gScore[cur.i] + legCost;

            if (tentative + 1e-12 < gScore[j]) {
              gScore[j] = tentative;
              cameFrom[j] = (int)cur.i;
              const double f = tentative + heuristic(j);
              open.push(QN{f, tentative, j});
            }
          }
        }
      }
    }
  }

  setStats(outStats, stats);
  return {};
}

std::vector<KRoute> plotKRoutesAStarCost(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        double costPerJump,
                                        double costPerLy,
                                        std::size_t k,
                                        std::size_t maxExpansionsPerSolve) {
  std::vector<KRoute> out;
  if (k == 0) return out;

  if (startId == 0 || goalId == 0) return out;
  if (maxJumpLy <= 0.0) return out;
  if (nodes.empty()) return out;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  // Degenerate: treat as hop-minimizing.
  if (costPerJump <= 0.0 && costPerLy <= 0.0) {
    costPerJump = 1.0;
    costPerLy = 0.0;
  }

  const auto idx = buildIndex(nodes);
  if (idx.find(startId) == idx.end() || idx.find(goalId) == idx.end()) return out;

  const auto grid = buildGrid(nodes, maxJumpLy);

  // First shortest path.
  {
    const auto base = aStarCostSolve(nodes, idx, grid, startId, goalId, maxJumpLy,
                                    costPerJump, costPerLy,
                                    nullptr, nullptr,
                                    maxExpansionsPerSolve);
    if (base.path.empty()) return out;

    out.push_back(KRoute{base.path, base.stats.hops, base.stats.distanceLy, base.stats.cost});
  }

  struct Candidate {
    std::vector<SystemId> path;
    int hops{0};
    double distanceLy{0.0};
    double cost{0.0};
  };

  const auto containsPath = [](const std::vector<KRoute>& routes, const std::vector<SystemId>& p) {
    for (const auto& r : routes) {
      if (r.path == p) return true;
    }
    return false;
  };

  const auto containsPathCand = [](const std::vector<Candidate>& cands, const std::vector<SystemId>& p) {
    for (const auto& c : cands) {
      if (c.path == p) return true;
    }
    return false;
  };

  const auto betterCand = [](const Candidate& a, const Candidate& b) {
    const double da = a.cost;
    const double db = b.cost;
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
      const SystemId spurNode = prev[i];

      // Root path includes spur node.
      std::vector<SystemId> root(prev.begin(), prev.begin() + (long long)i + 1);

      // Ban nodes in the root path *except* the spur node to enforce loopless paths.
      std::unordered_set<SystemId> bannedNodes;
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

      const auto spur = aStarCostSolve(nodes, idx, grid,
                                      spurNode, goalId,
                                      maxJumpLy,
                                      costPerJump, costPerLy,
                                      &bannedNodes,
                                      &bannedEdges,
                                      maxExpansionsPerSolve);

      if (spur.path.empty()) continue;

      // Combine root + spur (skip spurNode duplicate).
      std::vector<SystemId> total = root;
      if (spur.path.size() > 1) {
        total.insert(total.end(), spur.path.begin() + 1, spur.path.end());
      }

      if (total.size() < 2) continue;
      if (containsPath(out, total)) continue;
      if (containsPathCand(candidates, total)) continue;

      const double dist = computeDistanceLy(idx, nodes, total);
      const double cst  = computeCost(idx, nodes, total, costPerJump, costPerLy);

      const int hops = (int)total.size() - 1;
      candidates.push_back(Candidate{std::move(total), hops, dist, cst});
    }

    if (candidates.empty()) break;

    // Pick the best candidate.
    std::size_t bestIdx = 0;
    for (std::size_t i = 1; i < candidates.size(); ++i) {
      if (betterCand(candidates[i], candidates[bestIdx])) bestIdx = i;
    }

    Candidate best = std::move(candidates[bestIdx]);
    candidates.erase(candidates.begin() + (long long)bestIdx);

    out.push_back(KRoute{std::move(best.path), best.hops, best.distanceLy, best.cost});
  }

  return out;
}

std::vector<KRoute> plotKRoutesAStarHops(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        std::size_t k,
                                        std::size_t maxExpansionsPerSolve) {
  return plotKRoutesAStarCost(nodes, startId, goalId, maxJumpLy, 1.0, 0.0, k, maxExpansionsPerSolve);
}

double routeDistanceLy(const std::vector<SystemStub>& nodes,
                       const std::vector<SystemId>& route) {
  if (route.size() < 2) return 0.0;

  const auto idx = buildIndex(nodes);

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    sum += systemDistanceLy(nodes[itA->second], nodes[itB->second]);
  }
  return sum;
}

double routeCost(const std::vector<SystemStub>& nodes,
                 const std::vector<SystemId>& route,
                 double costPerJump,
                 double costPerLy) {
  if (route.size() < 2) return 0.0;

  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  const auto idx = buildIndex(nodes);

  double sum = 0.0;
  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto itA = idx.find(route[i]);
    const auto itB = idx.find(route[i + 1]);
    if (itA == idx.end() || itB == idx.end()) break;
    const double d = systemDistanceLy(nodes[itA->second], nodes[itB->second]);
    sum += costPerJump + costPerLy * d;
  }
  return sum;
}

bool validateRoute(const std::vector<SystemStub>& nodes,
                   const std::vector<SystemId>& route,
                   double maxJumpLy,
                   std::string* outError) {
  if (route.empty()) {
    if (outError) *outError = "route is empty";
    return false;
  }
  if (maxJumpLy <= 0.0) {
    if (outError) *outError = "maxJumpLy must be > 0";
    return false;
  }

  const auto idx = buildIndex(nodes);

  for (const auto id : route) {
    if (idx.find(id) == idx.end()) {
      if (outError) *outError = "route references unknown system id: " + std::to_string((unsigned long long)id);
      return false;
    }
  }

  for (std::size_t i = 0; i + 1 < route.size(); ++i) {
    const auto a = idx.find(route[i]);
    const auto b = idx.find(route[i + 1]);
    if (a == idx.end() || b == idx.end()) continue;
    const double d = systemDistanceLy(nodes[a->second], nodes[b->second]);
    if (d > maxJumpLy + 1e-9) {
      if (outError) {
        *outError = "jump " + std::to_string(i) + " exceeds maxJumpLy (" +
                    std::to_string(d) + " > " + std::to_string(maxJumpLy) + ")";
      }
      return false;
    }
  }

  return true;
}

} // namespace stellar::sim
