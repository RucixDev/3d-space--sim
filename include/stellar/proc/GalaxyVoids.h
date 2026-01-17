#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

namespace stellar::proc {

// -----------------------------------------------------------------------------
// GalaxyVoids
// -----------------------------------------------------------------------------
//
// A streaming-safe, deterministic "void bubble" influence field.
//
// This is the conceptual inverse of GalaxyClusters:
//  - clusters boost local density
//  - voids suppress local density
//
// The intent is stylized plausibility: large-scale cavities / shells that can
// create negative-space structure ("holes" and "filaments") in the generated
// galaxy without requiring any precomputed global state.
//
// Implementation sketch:
//  - Space is partitioned into a coarse 3D grid of cells (cellSizeLy).
//  - Each cell may contain one jittered void center, enabled with probability
//    chancePerCell (deterministic per cell).
//  - Influence at a point is computed by scanning a bounded neighborhood of
//    cells and taking the strongest radial falloff contribution.
//
// Typical use (in GalaxyGenerator):
//   mean *= clamp01(1 - voidStrength * void01)
//
// All IDs are stable and derived from (universeSeed, cellCoord).

struct GalaxyVoidsParams {
  // Coarse grid cell size (ly). <= 0 disables sampling.
  double cellSizeLy{4500.0};

  // Probability [0..1] that a given cell contains a void center.
  double chancePerCell{0.05};

  // Typical influence radius (ly).
  double radiusLy{2800.0};

  // Per-void radius jitter (multiplicative). 0.40 -> radius in ~[0.6..1.4] * radiusLy.
  double radiusJitter01{0.40};

  // Per-void intrinsic-strength jitter (multiplicative). 0.35 -> strength in ~[0.65..1.35].
  double strengthJitter01{0.35};

  // Falloff exponent (>0). Higher values tighten the bubble core.
  double falloffPower{2.0};
};

struct GalaxyVoidSample {
  bool hasVoid{false};

  // Stable id for the strongest influencing void at this position.
  core::u64 voidId{0};
  core::u64 voidSeed{0};

  // Void center parameters (valid when hasVoid is true).
  math::Vec3d centerLy{0,0,0};
  double radiusLy{0.0};
  double intrinsicStrength{1.0};

  // Influence at the query point in [0,1]. Includes intrinsicStrength.
  double void01{0.0};
};

// Sample the void field at `posLy`.
//
// NOTE: The returned void corresponds to the *strongest influence* at the query
// point (not necessarily the geometrically nearest center when radii vary).
GalaxyVoidSample sampleGalaxyVoids(core::u64 universeSeed,
                                   const math::Vec3d& posLy,
                                   const GalaxyVoidsParams& params = {});

} // namespace stellar::proc
