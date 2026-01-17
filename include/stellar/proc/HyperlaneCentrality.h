#pragma once

#include "stellar/core/Types.h"
#include "stellar/proc/HyperlaneRouter.h"
#include "stellar/proc/Hyperlanes.h"
#include "stellar/sim/Celestial.h"

#include <cstddef>
#include <vector>

namespace stellar::proc {

// -----------------------------------------------------------------------------
// Hyperlane betweenness centrality (graph analytics)
// -----------------------------------------------------------------------------
//
// Betweenness centrality measures how often a node/edge lies on *shortest paths*
// between pairs of nodes.
//
// We use it as a proxy for "traffic" / "chokepoint" importance in the procedural
// hyperlane network.
//
// Implementation:
//  - Weighted Brandes accumulation (Dijkstra per source).
//  - Optional source sampling for fast approximation.
//  - Deterministic source selection via stable hashing.
//
// Notes:
//  - For undirected graphs, we apply the standard 0.5 factor to avoid double
//    counting (because each unordered pair is counted from both endpoints when
//    summing over all sources).
//
// References:
//  - Ulrik Brandes, "A Faster Algorithm for Betweenness Centrality", 2001.

struct HyperlaneBetweennessResult {
  // Edge betweenness aligned with `net.edges` order.
  std::vector<double> edgeBetweenness;

  // Node betweenness aligned with `nodeIds`.
  std::vector<double> nodeBetweenness;
  std::vector<sim::SystemId> nodeIds;

  // Diagnostics / convenience.
  double maxEdge{0.0};
  double maxNode{0.0};
  std::size_t sourcesUsed{0};
};

struct HyperlaneBetweennessParams {
  // Shortest-path weight model.
  HyperlaneTravelParams travel{};

  // Number of source nodes to sample.
  //  - 0 => exact (all sources)
  //  - 1..N => approximate using this many sources.
  std::size_t sampleSources{64};

  // Extra seed mixed into deterministic source selection.
  core::u64 sampleSeed{0};

  // Optional guard to bound work in pathological cases.
  //  - 0 => unbounded
  std::size_t maxExpansionsPerSource{0};

  // Relative epsilon used when comparing floating distances for path ties.
  // (Smaller => fewer ties; larger => more multi-shortest paths.)
  double tieEps{1e-9};
};

// Estimate betweenness centrality for a hyperlane network.
//
// Inputs:
//  - nodes: system stubs for the graph nodes.
//  - net: undirected lane network (each edge appears once, a<b).
//  - params: weighting + sampling controls.
//
// Returns:
//  - edgeBetweenness sized `net.edges.size()` (same order).
//  - nodeBetweenness sized `nodeIds.size()` with nodeIds sorted ascending.
HyperlaneBetweennessResult estimateHyperlaneBetweennessCentrality(
  const std::vector<sim::SystemStub>& nodes,
  const HyperlaneNetwork& net,
  const HyperlaneBetweennessParams& params = {});

} // namespace stellar::proc
