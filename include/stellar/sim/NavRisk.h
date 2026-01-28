#pragma once

#include "stellar/core/Types.h"

#include "stellar/sim/Celestial.h"
#include "stellar/sim/SystemEvents.h"
#include "stellar/sim/SystemSecurityDynamics.h"

#include <span>
#include <vector>

namespace stellar::sim {

class Universe;
struct StarSystem;

// -----------------------------------------------------------------------------
// Navigation Risk Model (headless)
// -----------------------------------------------------------------------------
//
// This module provides a simple, deterministic "travel risk" signal derived
// from system security conditions (security / piracy / traffic / contest).
//
// Unlike MissionRisk (which is contract-profile aware), this is intended for
// route planning and cockpit tooling: "Should I avoid that system right now?"
//
// The model is deliberately lightweight:
//  - deterministically derived from Universe seed + (optional) persistent deltas
//  - no persistent simulation step required
//  - callable from tests/tools/headless builds
//
// Typical usage:
//  1) computeNavRisk01ForNodes() for a nav graph neighborhood
//  2) pass risk01PerNode into plotRouteAStarCostRisk() / computeNavRouteBatchCostRisk()

struct NavSystemRiskSignals {
  // Raw effective condition channels.
  double security01{0.5};
  double piracy01{0.5};
  double traffic01{0.5};
  double contest01{0.0};

  // Combined, normalized risk signal in [0,1].
  double risk01{0.5};
};

struct NavRiskParams {
  // Relative weights for the combined risk score.
  // The weights do not need to sum to 1; the implementation normalizes.
  double piracyWeight{0.50};      // more piracy => more risk
  double insecurityWeight{0.25};  // lower security => more risk
  double lowTrafficWeight{0.15};  // lower traffic => more risk
  double contestWeight{0.10};     // contested systems are more volatile

  // Clamp final risk.
  double minRisk{0.0};
  double maxRisk{1.0};

  // Whether to include deterministic SystemEvents in the effective profile.
  bool includeSystemEvents{true};

  // Whether to apply persistent security deltas from the provided delta list.
  bool includeSecurityDeltas{true};

  // Tunables for the underlying system conditions snapshot.
  SystemSecurityDynamicsParams dynamicsParams{};
  SystemEventParams eventParams{};
};

// Compute travel-risk signals for a specific system.
//
// If `delta` is null, only the deterministic baseline (and optional event)
// layers are considered.
NavSystemRiskSignals computeNavSystemRiskSignals(const Universe& universe,
                                                 const StarSystem& sys,
                                                 double timeDays,
                                                 const SystemSecurityDeltaState* delta = nullptr,
                                                 const NavRiskParams& params = {});

// Compute risk01 values aligned to `nodes`.
//
// - risk01PerNode[i] corresponds to nodes[i]
// - Unknown system ids (0) yield risk 0.5.
std::vector<double> computeNavRisk01ForNodes(const Universe& universe,
                                             std::span<const SystemStub> nodes,
                                             double timeDays,
                                             std::span<const SystemSecurityDeltaState> securityDeltas = {},
                                             const NavRiskParams& params = {});

// Distance-weighted average risk along a route.
//
// If `risk01PerNode.size() != nodes.size()`, unknown entries are treated as 0.
double estimateRouteAvgRisk01(std::span<const SystemStub> nodes,
                              std::span<const SystemId> route,
                              std::span<const double> risk01PerNode);

// Maximum node risk along a route.
double estimateRouteMaxRisk01(std::span<const SystemStub> nodes,
                              std::span<const SystemId> route,
                              std::span<const double> risk01PerNode);

} // namespace stellar::sim
