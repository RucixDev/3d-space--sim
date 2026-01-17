#pragma once

#include "stellar/proc/HyperlaneRouter.h"

#include <cstddef>
#include <vector>

namespace stellar::proc {

// -----------------------------------------------------------------------------
// K-shortest hyperlane routes (loopless)
// -----------------------------------------------------------------------------
//
// Why this exists:
//  - HyperlaneRouter provides single-source shortest paths for procedural systems.
//  - UI tools (Procedural Galaxy Lab) and gameplay systems often want *alternatives*:
//      * "fastest" vs "safest" trade-offs
//      * route variety for missions / NPC logistics
//      * debugging hazard modulation in Hyperlanes
//
// Implementation notes:
//  - Uses Yen's algorithm to enumerate K loopless shortest paths.
//  - Uses an A* point-to-point solver as the shortest-path subroutine.
//  - Deterministic ordering: (cost asc) then (lexicographic SystemId path).

struct HyperlaneKRoute {
  std::vector<sim::SystemId> path;
  HyperlanePathMetrics metrics;
};

// Compute K loopless routes between startId and goalId over the supplied hyperlane
// edge set.
//
// Requirements:
//  - `nodes` must contain stubs for both startId and goalId.
//  - `net` edges should refer to ids present in `nodes` (unknown endpoints are ignored).
//
// Returns:
//  - Up to K routes ordered by increasing metrics.costLy (then by path ids).
//  - Empty when no route exists.
std::vector<HyperlaneKRoute> plotKHyperlaneRoutesAStarCost(const std::vector<sim::SystemStub>& nodes,
                                                          const HyperlaneNetwork& net,
                                                          sim::SystemId startId,
                                                          sim::SystemId goalId,
                                                          const HyperlaneTravelParams& travel,
                                                          std::size_t k,
                                                          std::size_t maxExpansionsPerSolve = 250000);

} // namespace stellar::proc
