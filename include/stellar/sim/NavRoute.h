#pragma once

#include "stellar/sim/Celestial.h"

#include <cstddef>
#include <string>
#include <vector>

namespace stellar::sim {

// Optional diagnostics about a planned route.
struct RoutePlanStats {
  int expansions{0};      // nodes popped from the open set
  int visited{0};         // nodes closed/expanded
  int hops{0};            // route.size() - 1
  double distanceLy{0.0}; // sum of straight-line leg distances
  double cost{0.0};       // total A* path cost for the chosen cost model
  bool reached{false};
};

// A* route planner that minimizes hop count (each jump cost = 1).
// The heuristic is ceil(remainingDistance / maxJumpLy), which is admissible.
//
// NOTE: `nodes` must include stubs for both startId and goalId.
// A common pattern is to pass the vector returned by Universe::queryNearby(...).
//
// Returns a list of SystemId including start and goal, or empty if no route found.
std::vector<SystemId> plotRouteAStarHops(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        RoutePlanStats* outStats = nullptr,
                                        std::size_t maxExpansions = 250000);

// A* route planner with a simple weighted per-leg cost model:
//
//   legCost = costPerJump + costPerLy * legDistanceLy
//
// This allows different "navigation feels" without changing the graph:
//
//   - Min-distance:  costPerJump = 0,         costPerLy = 1
//   - Fuel-like:     costPerJump = fuelBase,  costPerLy = fuelPerLy
//   - Hybrid:        increase costPerJump to prefer fewer hops
//
// Heuristic: ceil(remainingDistance/maxJumpLy) * costPerJump + remainingDistance * costPerLy
// This is admissible as long as costs are non-negative.
//
// Returns a list of SystemId including start and goal, or empty if no route found.
std::vector<SystemId> plotRouteAStarCost(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        double costPerJump,
                                        double costPerLy,
                                        RoutePlanStats* outStats = nullptr,
                                        std::size_t maxExpansions = 250000);

// Result for K-shortest route planning.
struct KRoute {
  std::vector<SystemId> path;
  int hops{0};
  double distanceLy{0.0};
  double cost{0.0};
};

// K-shortest *loopless* routes between startId and goalId using Yen's algorithm.
//
// - Each route is a list of SystemId including start and goal.
// - Results are ordered by increasing total cost, then deterministically by the
//   path's SystemId sequence.
// - For k=1, this behaves similarly to plotRouteAStarCost (but without returning
//   detailed per-solve expansion diagnostics).
std::vector<KRoute> plotKRoutesAStarCost(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        double costPerJump,
                                        double costPerLy,
                                        std::size_t k,
                                        std::size_t maxExpansionsPerSolve = 250000);

// Convenience wrapper for hop-minimizing K-shortest paths.
std::vector<KRoute> plotKRoutesAStarHops(const std::vector<SystemStub>& nodes,
                                        SystemId startId,
                                        SystemId goalId,
                                        double maxJumpLy,
                                        std::size_t k,
                                        std::size_t maxExpansionsPerSolve = 250000);

// Helper: total straight-line length of a route in ly.
// Returns 0 for empty/single-node routes.
double routeDistanceLy(const std::vector<SystemStub>& nodes,
                       const std::vector<SystemId>& route);

// Helper: total cost of a route using the same cost model as plotRouteAStarCost.
// Returns 0 for empty/single-node routes.
double routeCost(const std::vector<SystemStub>& nodes,
                 const std::vector<SystemId>& route,
                 double costPerJump,
                 double costPerLy);

// Helper: validate that a route is contiguous (each hop <= maxJumpLy) and that all ids exist.
// Useful for tests/tooling.
bool validateRoute(const std::vector<SystemStub>& nodes,
                   const std::vector<SystemId>& route,
                   double maxJumpLy,
                   std::string* outError = nullptr);

} // namespace stellar::sim
