#pragma once

#include "stellar/core/Types.h"
#include "stellar/proc/Hyperlanes.h"

#include <unordered_map>
#include <vector>

namespace stellar::proc {

// -----------------------------------------------------------------------------
// HyperlaneRouter
// -----------------------------------------------------------------------------
//
// A small, headless routing helper for proc::HyperlaneNetwork.
//
// Purpose:
//  - Convert a local hyperlane edge set into a useful *travel metric*.
//  - Enable procedural systems (trade routes, mission seeding, NPC logistics)
//    to prefer high-bandwidth / low-risk corridors.
//
// Design notes:
//  - Single-source shortest paths (Dijkstra) over an undirected weighted graph.
//  - Deterministic tie-breaking (cost, distance, hops, bottleneck, predecessor).
//  - "cost" is a distance-like metric (units: ly) produced from edge distance
//    modulated by risk and bandwidth.

struct HyperlaneTravelParams {
  // Additional cost applied for risk.
  //  - 0.0 => ignore risk
  //  - 1.0 => cost *= (1 + risk01)
  double riskWeight{0.60};

  // How strongly bandwidth reduces cost.
  //  - 0.0 => ignore bandwidth
  //  - 1.0 => divide by bandwidth factor
  double bandwidthBias{0.65};

  // Minimum bandwidth factor used when bandwidthBias == 1.
  // Prevents extreme costs for low-bandwidth lanes.
  double minBandwidthFactor{0.35};
};

struct HyperlanePathMetrics {
  bool reachable{false};

  // "Effective" travel cost in ly-like units (this is what routing minimizes).
  double costLy{0.0};

  // Physical length of the chosen lane path (sum of edge distances).
  double distanceLy{0.0};

  // Compound risk of the path in [0,1].
  // Computed as: 1 - Π(1 - edgeRisk).
  double risk01{0.0};

  // Minimum bandwidth along the path in [0,1].
  double bottleneckBandwidth01{0.0};

  // Number of hyperlane edges in the path.
  int hops{0};
};

class HyperlaneRouter {
public:
  HyperlaneRouter() = default;
  HyperlaneRouter(const HyperlaneNetwork& net, const std::vector<sim::SystemStub>& stubs);

  // Rebuild internal adjacency from a new node/edge set.
  void reset(const HyperlaneNetwork& net, const std::vector<sim::SystemStub>& stubs);

  // Run a single-source route computation.
  // Returns false if origin is not present.
  bool compute(sim::SystemId originId, const HyperlaneTravelParams& params = {});

  // Query metrics for a destination.
  // If unreachable, reachable=false and other fields are zero.
  HyperlanePathMetrics metricsTo(sim::SystemId targetId) const;

  // Reconstruct the chosen path from the last compute() call.
  //
  // Output is a list of system ids: [origin, ..., target]. Returns false when
  //  - compute() has not been called,
  //  - target is unknown,
  //  - or the target is unreachable.
  bool buildPathTo(sim::SystemId targetId, std::vector<sim::SystemId>& outPath) const;

private:
  struct Adj {
    int to{0};
    double distanceLy{0.0};
    double bandwidth01{0.0};
    double risk01{0.0};
  };

  sim::SystemId origin_{0};
  bool computed_{false};

  std::vector<sim::SystemId> ids_;
  std::unordered_map<sim::SystemId, int> idToIndex_;
  std::vector<std::vector<Adj>> adj_;

  // Results for the last compute() call.
  std::vector<double> bestCost_;
  std::vector<double> bestDistance_;
  std::vector<double> bestRiskNot_; // Π(1 - risk)
  std::vector<double> bestBottleneckBw_;
  std::vector<int> bestHops_;
  std::vector<int> bestPrev_;
};

} // namespace stellar::proc
