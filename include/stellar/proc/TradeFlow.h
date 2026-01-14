#pragma once

#include "stellar/core/Types.h"
#include "stellar/proc/HyperlaneRouter.h"
#include "stellar/proc/TradeProfile.h"

#include <cstddef>
#include <vector>

namespace stellar::proc {

// -----------------------------------------------------------------------------
// TradeFlow: macro interstellar traffic estimation (headless)
// -----------------------------------------------------------------------------
//
// TradeRoutes / TradeIntel provide *which* systems are likely trade partners and
// what commodities fit the relationship. However, a galaxy also benefits from a
// higher-level concept of *corridors*:
//
//   - Which hyperlanes are expected to be busy?
//   - Which systems are transit hubs vs quiet backwaters?
//
// This module computes a deterministic, macro-scale traffic field over a local
// hyperlane network by:
//   1) Assigning each system an "economic mass" derived from TradeProfile.
//   2) Creating pairwise "trade potential" using a gravity-like model:
//        flow_ij ~ mass_i * mass_j * compatibility_ij / travelCost_ij^k
//   3) Routing that flow along the lowest-cost hyperlane path (HyperlaneRouter).
//   4) Accumulating per-edge flow and per-node incident traffic.
//
// Notes:
//  - This is intentionally *not* a full economic simulation.
//  - Output is useful for UI probes, encounter spawning heuristics, and future
//    dynamic systems (pirate belts, patrol density, mission seeding, etc.).

struct TradeFlowParams {
  // If stubs.size() <= fullPairLimit, compute all i<j pairs.
  // Otherwise, compute an approximate field by sampling pairs from a smaller
  // set of economically important source systems.
  std::size_t fullPairLimit{512};

  // When sampling, number of source systems to use (picked by economic mass).
  int sampleSourceCount{96};

  // When sampling, number of random partners per source.
  int samplePartnersPerSource{96};

  // Gravity model distance exponent (typical empirical values are ~1).
  // Larger => more localized traffic.
  double gravityExponent{1.15};

  // Floor for travel cost in the gravity term.
  double minCostLy{1.0};

  // Commodity compatibility weight.
  // Pair factor = (1 - w) + w * compatibility_ij.
  double commodityCompatibilityWeight{0.70};

  // Economic mass weights (all fields in [0,1]).
  double hubWeight{0.55};
  double populationWeight{0.30};
  double wealthWeight{0.10};
  double industryWeight{0.05};
  double massFloor{0.15};

  // Routing settings (risk/bandwidth affect travel cost).
  HyperlaneTravelParams travelParams{};
};

struct TradeFlowNode {
  sim::SystemId id{0};

  // Sum of incident lane flow (absolute units).
  double traffic{0.0};

  // traffic normalized by max traffic in the output node set.
  double traffic01{0.0};
};

struct TradeFlowEdge {
  // Undirected endpoints (a < b).
  sim::SystemId a{0};
  sim::SystemId b{0};

  // Aggregate lane flow (absolute units).
  double flow{0.0};

  // flow normalized by max flow in the output edge set.
  double flow01{0.0};
};

struct TradeFlowNetwork {
  std::vector<TradeFlowNode> nodes;
  std::vector<TradeFlowEdge> edges;

  // Sum of pairwise routed flow added to the network (absolute units).
  double totalFlow{0.0};
};

// Compute a deterministic macro traffic field over a local hyperlane network.
//
// Requirements:
//  - stubs.size() == profiles.size()
//  - hyperlanes edges should reference systems present in stubs (any missing
//    endpoints are ignored by the router).
TradeFlowNetwork computeTradeFlow(core::u64 universeSeed,
                                  const std::vector<sim::SystemStub>& stubs,
                                  const std::vector<TradeProfile>& profiles,
                                  const HyperlaneNetwork& hyperlanes,
                                  const TradeFlowParams& params = {});

// Convenience overload: generate TradeProfile for each stub.
TradeFlowNetwork computeTradeFlow(core::u64 universeSeed,
                                  const std::vector<sim::SystemStub>& stubs,
                                  const HyperlaneNetwork& hyperlanes,
                                  const TradeFlowParams& params = {});

} // namespace stellar::proc
