#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/FactionProfile.h"
#include "stellar/sim/System.h"

#include <utility>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// System Security Model (deterministic, headless)
// -----------------------------------------------------------------------------
//
// Many parts of the prototype (encounter scheduling, mission rewards, traffic)
// want a small set of "world state" knobs that vary meaningfully per system,
// but remain deterministic across runs.
//
// This module derives a few 0..1 metrics from:
//   - which factions own stations in the system ("control" and "contestedness")
//   - the controlling faction's procedural traits (FactionProfile)
//
// The outputs are intentionally *not* a full simulation. They're compact
// signals that other systems can consume without pulling in the whole Universe.

// Weighted ownership breakdown for a star system.
struct SystemControl {
  core::u32 controllingFactionId{0};

  // Fraction of station "weight" owned by the controlling faction.
  // 1.0 => totally controlled, 0.5 => perfectly contested.
  double controlFrac{0.0};

  // 0..1 measure of how contested the system is.
  // 0 => one faction dominates, 1 => highly contested.
  double contest01{0.0};

  // Sorted by weight (desc), then factionId (asc).
  std::vector<std::pair<core::u32, double>> factionWeights{};
};

// Compute a simple deterministic "control" score from station ownership.
// If the system has no stations, controllingFactionId is 0.
SystemControl computeSystemControl(const StarSystem& sys);

// Compact, derived system-level security/risk metrics.
struct SystemSecurityProfile {
  core::u32 controllingFactionId{0};
  FactionProfile traits{};

  double controlFrac{0.0};
  double contest01{0.0};

  // 0..1 values where 0.5 is the "baseline".
  // These defaults intentionally match the old encounter tuning when used
  // as multipliers.
  double security01{0.5}; // higher => more effective security response
  double piracy01{0.5};   // higher => more piracy pressure
  double traffic01{0.5};  // higher => more civilian/trader traffic
};

// Deterministically derive a SystemSecurityProfile from a generated StarSystem.
SystemSecurityProfile systemSecurityProfile(core::u64 universeSeed, const StarSystem& sys);

} // namespace stellar::sim
