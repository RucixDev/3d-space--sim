#include "stellar/sim/NavRouteBatch.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <unordered_map>

namespace stellar::sim {

namespace {

static double systemDistanceLy(const SystemStub& a, const SystemStub& b) {
  return (a.posLy - b.posLy).length();
}

struct CellCoord {
  long long x{0};
  long long y{0};
  long long z{0};

  bool operator==(const CellCoord& o) const {
    return x == o.x && y == o.y && z == o.z;
  }
};

static std::size_t hashCombine(std::size_t h, std::size_t v) {
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

static std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash>
buildGrid(const std::vector<SystemStub>& nodes, double cellSize) {
  std::unordered_map<CellCoord, std::vector<std::size_t>, CellHash> grid;
  grid.reserve(nodes.size());

  for (std::size_t i = 0; i < nodes.size(); ++i) {
    grid[cellFor(nodes[i].posLy, cellSize)].push_back(i);
  }

  return grid;
}

} // namespace

NavRouteBatch computeNavRouteBatchCost(const std::vector<SystemStub>& nodes,
                                      SystemId startId,
                                      double maxJumpLy,
                                      double costPerJump,
                                      double costPerLy,
                                      std::size_t maxExpansions) {
  NavRouteBatch out{};

  if (startId == 0) return out;
  if (maxJumpLy <= 0.0) return out;
  if (nodes.empty()) return out;

  if (!std::isfinite(costPerJump)) costPerJump = 0.0;
  if (!std::isfinite(costPerLy)) costPerLy = 0.0;
  if (costPerJump < 0.0) costPerJump = 0.0;
  if (costPerLy < 0.0) costPerLy = 0.0;

  out.startId = startId;
  out.maxJumpLy = maxJumpLy;
  out.costPerJump = costPerJump;
  out.costPerLy = costPerLy;

  out.index = buildIndex(nodes);
  const auto itS = out.index.find(startId);
  if (itS == out.index.end()) {
    // Caller did not include startId in nodes.
    return out;
  }

  const std::size_t start = itS->second;
  const std::size_t N = nodes.size();

  out.parent.assign(N, -1);
  out.bestCost.assign(N, std::numeric_limits<double>::infinity());
  out.hops.assign(N, std::numeric_limits<int>::max());
  out.distanceLy.assign(N, std::numeric_limits<double>::infinity());

  std::vector<char> closed(N, 0);
  const auto grid = buildGrid(nodes, maxJumpLy);

  struct QN {
    double g{0.0};
    int hops{0};
    std::size_t i{0};
  };

  struct Cmp {
    bool operator()(const QN& a, const QN& b) const {
      // min-heap behavior via reversed comparator
      if (a.g != b.g) return a.g > b.g;
      if (a.hops != b.hops) return a.hops > b.hops;
      return a.i > b.i;
    }
  };

  std::priority_queue<QN, std::vector<QN>, Cmp> open;

  out.bestCost[start] = 0.0;
  out.hops[start] = 0;
  out.distanceLy[start] = 0.0;
  open.push(QN{0.0, 0, start});

  int expansions = 0;
  int visited = 0;

  while (!open.empty() && (std::size_t)expansions < maxExpansions) {
    const QN cur = open.top();
    open.pop();

    if (closed[cur.i]) continue;
    closed[cur.i] = 1;
    ++expansions;
    ++visited;

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
            if (nid == 0) continue;

            const double d = systemDistanceLy(nodes[cur.i], nodes[j]);
            if (d > maxJumpLy + 1e-9) continue;

            const double legCost = costPerJump + costPerLy * d;
            const double tentative = out.bestCost[cur.i] + legCost;
            const int tentativeHops = out.hops[cur.i] + 1;
            const double tentativeDist = out.distanceLy[cur.i] + d;

            bool better = (tentative + 1e-12 < out.bestCost[j]);
            const bool tie = (!better && std::abs(tentative - out.bestCost[j]) <= 1e-12);

            if (tie) {
              // Deterministic tie-break:
              //  1) prefer fewer hops
              //  2) then prefer smaller parent SystemId
              if (tentativeHops < out.hops[j]) {
                better = true;
              } else if (tentativeHops == out.hops[j]) {
                const SystemId newParent = nodes[cur.i].id;
                SystemId oldParent = std::numeric_limits<SystemId>::max();
                if (out.parent[j] >= 0) oldParent = nodes[(std::size_t)out.parent[j]].id;
                if (newParent < oldParent) better = true;
              }
            }

            if (better) {
              out.bestCost[j] = tentative;
              out.parent[j] = (int)cur.i;
              out.hops[j] = tentativeHops;
              out.distanceLy[j] = tentativeDist;
              open.push(QN{tentative, tentativeHops, j});
            }
          }
        }
      }
    }
  }

  out.stats.expansions = expansions;
  out.stats.visited = visited;
  return out;
}

std::vector<SystemId> routeFromBatch(const std::vector<SystemStub>& nodes,
                                    const NavRouteBatch& batch,
                                    SystemId goalId) {
  if (goalId == 0) return {};
  if (batch.startId == 0) return {};

  const auto itGoal = batch.index.find(goalId);
  const auto itStart = batch.index.find(batch.startId);
  if (itGoal == batch.index.end() || itStart == batch.index.end()) return {};

  const std::size_t goal = itGoal->second;
  if (goal >= batch.bestCost.size()) return {};
  if (!std::isfinite(batch.bestCost[goal])) return {};

  std::vector<SystemId> path;
  path.reserve(64);

  int at = (int)goal;
  for (std::size_t steps = 0; steps < nodes.size(); ++steps) {
    if (at < 0) break;
    const std::size_t i = (std::size_t)at;
    if (i >= nodes.size()) break;
    path.push_back(nodes[i].id);
    at = batch.parent[i];
  }

  std::reverse(path.begin(), path.end());

  if (path.empty()) return {};
  if (path.front() != batch.startId) return {};
  if (path.back() != goalId) return {};

  return path;
}

} // namespace stellar::sim
