#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Celestial.h"

#include <cstddef>
#include <vector>

namespace stellar::proc {

// -----------------------------------------------------------------------------
// Procedural Hyperlanes (galaxy-scale lane network)
// -----------------------------------------------------------------------------
//
// Motivation:
//  - The universe already supports on-demand procedural generation of systems.
//  - Many gameplay/UI systems benefit from a *sparse* galaxy navigation graph:
//      * trade corridor visualization
//      * mission seeding (origin/destination plausibility)
//      * NPC logistics (future)
//      * "civilized" vs "frontier" travel feel
//
// This module generates an undirected, deterministic lane network for an
// arbitrary set of system stubs.
//
// Design goals:
//  - Deterministic (depends only on universeSeed + stub fields + params)
//  - Headless (proc-only)
//  - Cheap for local queries (O(n log n) backbone; candidates are kNN-limited)
//
// NOTE: The produced graph is local to the provided stub set. If you query a
// different neighborhood, boundary edges may differ (by design).

struct HyperlaneEdge {
  sim::SystemId a{0};
  sim::SystemId b{0};

  // Straight-line system-to-system distance.
  double distanceLy{0.0};

  // Macro "capacity" signal in [0,1]. Higher means the lane is likely to carry
  // more traffic and be better patrolled.
  double bandwidth01{0.0};

  // Macro risk proxy in [0,1]. Higher means more piracy / instability.
  double risk01{0.0};
};

struct HyperlaneHazardParams {
  // Enable the galaxy-scale hazard field modulation.
  //
  // When enabled, hyperlane cost/bandwidth/risk are influenced by a smooth
  // hazard field (nebulae / storms / rifts). This can create emergent
  // “safe corridors” and “stormy dead-zones” without hardcoding any layout.
  bool enabled{false};

  // Hazard field time (days). The field drifts deterministically with time.
  double timeDays{0.0};

  // Drift speed (ly/day). Higher values make hazards evolve faster.
  double driftLyPerDay{0.65};

  // Field scales (ly). Larger => smoother, broader structures.
  double nebulaScaleLy{650.0};
  double stormScaleLy{180.0};

  // Samples taken along an edge to approximate “integrated” hazard.
  int edgeSamples{3};

  // How strongly hazards bias MST selection (network topology).
  //  cost *= (1 + hazardToCost * hazardAvg)
  double hazardToCost{0.45};

  // How strongly hazards increase risk.
  //  risk = clamp01(risk + hazardToRisk * hazardAvg)
  double hazardToRisk{0.28};

  // How strongly hazards reduce bandwidth.
  //  bw = clamp01(bw * (1 - hazardToBandwidth * hazardAvg))
  double hazardToBandwidth{0.20};
};

struct HyperlaneParams {
  // Candidate neighbor search radius (ly). Edges longer than this are not
  // considered unless forceConnected is enabled.
  double maxNeighborDistanceLy{16.0};

  // Each system proposes edges to its K nearest neighbors.
  int neighborK{4};

  // If true, the generator attempts to connect all components by adding
  // longer "bridge" edges (best-effort).
  bool forceConnected{true};

  // Ensure each node has at least this degree (best-effort).
  int minDegree{2};

  // After building the backbone + min-degree edges, add additional edges with
  // a deterministic hash-based acceptance test.
  double extraEdgeChance{0.22};

  // Cell size (ly) used for sampling the GalaxyRegions layer. Set <= 0 to
  // disable region modulation.
  double regionCellSizeLy{900.0};

  // Optional galaxy-scale hazard modulation (nebulae / storms).
  HyperlaneHazardParams hazards{};

  // Optional cap on returned edges (0 means "no cap").
  std::size_t maxEdges{0};
};

struct HyperlaneNetwork {
  std::vector<HyperlaneEdge> edges;
};

// Build a deterministic hyperlane network for the provided system stubs.
//
// Output ordering is deterministic:
//   bandwidth desc, distance asc, a id asc, b id asc
//
// Each edge appears once with a < b.
HyperlaneNetwork generateHyperlaneNetwork(core::u64 universeSeed,
                                          const std::vector<sim::SystemStub>& stubs,
                                          const HyperlaneParams& params = {});

} // namespace stellar::proc
