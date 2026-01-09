#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/ResourceField.h"

namespace stellar::sim {

// Deterministic traits derived from (universeSeed, asteroidId, fieldKind).
struct MiningAsteroidTraits {
  // True if this asteroid can contain volatile seams that may fracture during extraction.
  bool volatilePocket{false};

  // If volatilePocket is true, this is the remaining-fraction threshold where a fracture triggers.
  // (0 means "no fracture event".)
  double fractureFrac{0.0};
};

// Compute stable asteroid traits (no RNG state; purely hash-based).
MiningAsteroidTraits miningAsteroidTraits(core::u64 universeSeed,
                                         core::u64 asteroidId,
                                         ResourceFieldKind fieldKind);

// Efficiency curve for mining lasers.
//
// Returns a multiplier in [minEfficiency, 1]. At close range (<= fullEfficiencyFrac * rangeKm)
// the multiplier is 1.0 and it linearly falls off toward minEfficiency at max range.
double miningEfficiency(double distKm,
                        double rangeKm,
                        double fullEfficiencyFrac = 0.30,
                        double minEfficiency = 0.25);

struct MiningHitInput {
  core::u64 universeSeed{0};
  core::u64 asteroidId{0};
  ResourceFieldKind fieldKind{ResourceFieldKind::OreBelt};

  // Shooter -> hit distance.
  double distKm{0.0};
  double rangeKm{180000.0};

  // Baseline units extracted per hit at full efficiency (before prospect bonus).
  double baseUnitsPerHit{10.0};

  // True if the asteroid has been prospected.
  bool prospected{false};

  // Asteroid state (for fracture events).
  double baseUnits{120.0};
  double remainingUnits{120.0};
  bool fractureAlreadyTriggered{false};
};

struct MiningHitResult {
  // How many units were extracted this hit (already clamped by remainingUnits).
  double extractedUnits{0.0};

  // Efficiency multiplier used for this hit.
  double efficiency{0.0};

  // Traits (mirrors miningAsteroidTraits result).
  bool volatilePocket{false};
  double fractureFrac{0.0};

  // True if this hit crosses the deterministic fracture threshold.
  bool fractureTriggered{false};
};

MiningHitResult computeMiningHit(const MiningHitInput& in);

} // namespace stellar::sim
