#pragma once

#include "stellar/sim/Celestial.h"

#include <cstddef>
#include <cmath>
#include <limits>
#include <span>
#include <unordered_map>
#include <vector>

namespace stellar::sim {

// Single-source shortest paths over the implicit "range graph" induced by maxJumpLy.
//
// Compared to calling A* repeatedly for many different goals, this is often much
// faster when you need routes from one origin to many destinations.
struct NavRouteBatchStats {
  int expansions{0}; // nodes popped from the open set
  int visited{0};    // nodes finalized
};

struct NavRouteBatch {
  SystemId startId{0};
  double maxJumpLy{0.0};
  double costPerJump{0.0};
  double costPerLy{0.0};
  double riskWeightPerLy{0.0};

  // Maps SystemId -> index into the input `nodes` vector used for the solve.
  std::unordered_map<SystemId, std::size_t> index;

  // Per-node result arrays (aligned to the input `nodes` vector).
  std::vector<int> parent;       // parent node index, or -1 if none
  std::vector<double> bestCost;  // cost from start, or +inf if unreachable
  std::vector<int> hops;         // hops along the chosen best-cost path
  std::vector<double> distanceLy;// straight-line distance along the chosen best-cost path

  NavRouteBatchStats stats{};

  bool has(SystemId id) const { return index.find(id) != index.end(); }

  bool reachable(SystemId id) const {
    const auto it = index.find(id);
    if (it == index.end()) return false;
    const std::size_t i = it->second;
    if (i >= bestCost.size()) return false;
    return std::isfinite(bestCost[i]);
  }

  double costTo(SystemId id) const {
    const auto it = index.find(id);
    if (it == index.end()) return std::numeric_limits<double>::infinity();
    const std::size_t i = it->second;
    if (i >= bestCost.size()) return std::numeric_limits<double>::infinity();
    return bestCost[i];
  }

  int hopsTo(SystemId id) const {
    const auto it = index.find(id);
    if (it == index.end()) return 0;
    const std::size_t i = it->second;
    if (i >= hops.size()) return 0;
    return hops[i];
  }

  double distanceTo(SystemId id) const {
    const auto it = index.find(id);
    if (it == index.end()) return std::numeric_limits<double>::infinity();
    const std::size_t i = it->second;
    if (i >= distanceLy.size()) return std::numeric_limits<double>::infinity();
    return distanceLy[i];
  }
};

// Compute shortest paths from `startId` to every node in `nodes`.
//
// cost model: legCost = costPerJump + costPerLy * legDistanceLy
//
// Returns a NavRouteBatch containing parent pointers for route reconstruction.
NavRouteBatch computeNavRouteBatchCost(const std::vector<SystemStub>& nodes,
                                      SystemId startId,
                                      double maxJumpLy,
                                      double costPerJump,
                                      double costPerLy,
                                      std::size_t maxExpansions = 250000);

// Risk-aware variant.
//
// cost model: legCost = costPerJump + costPerLy * d + riskWeightPerLy * avgRisk01 * d
// where avgRisk01 is the average of the endpoint risk values for the leg.
//
// If risk01PerNode.size() != nodes.size(), all risks are treated as 0.
NavRouteBatch computeNavRouteBatchCostRisk(const std::vector<SystemStub>& nodes,
                                          SystemId startId,
                                          double maxJumpLy,
                                          double costPerJump,
                                          double costPerLy,
                                          double riskWeightPerLy,
                                          std::span<const double> risk01PerNode,
                                          std::size_t maxExpansions = 250000);

// Convenience wrappers.
inline NavRouteBatch computeNavRouteBatchHops(const std::vector<SystemStub>& nodes,
                                             SystemId startId,
                                             double maxJumpLy,
                                             std::size_t maxExpansions = 250000) {
  return computeNavRouteBatchCost(nodes, startId, maxJumpLy,
                                 /*costPerJump=*/1.0,
                                 /*costPerLy=*/0.0,
                                 maxExpansions);
}

inline NavRouteBatch computeNavRouteBatchDistance(const std::vector<SystemStub>& nodes,
                                                 SystemId startId,
                                                 double maxJumpLy,
                                                 std::size_t maxExpansions = 250000) {
  return computeNavRouteBatchCost(nodes, startId, maxJumpLy,
                                 /*costPerJump=*/0.0,
                                 /*costPerLy=*/1.0,
                                 maxExpansions);
}

// Reconstruct a route from the batch result.
// Returns an empty vector if no path exists or if goalId/startId are missing.
std::vector<SystemId> routeFromBatch(const std::vector<SystemStub>& nodes,
                                    const NavRouteBatch& batch,
                                    SystemId goalId);

} // namespace stellar::sim
