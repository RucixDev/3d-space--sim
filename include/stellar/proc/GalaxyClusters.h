#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

namespace stellar::proc {

// -----------------------------------------------------------------------------
// GalaxyClusters
// -----------------------------------------------------------------------------
//
// A streaming-safe, deterministic "star cluster" influence field.
//
// The galaxy generator can use this as an additional density multiplier to
// create discrete high-density pockets (open clusters / associations) without
// requiring any global, precomputed list of cluster centers.
//
// Implementation sketch:
//  - Space is partitioned into a coarse 3D grid of cells (cellSizeLy).
//  - Each cell may contain one jittered cluster center, enabled with probability
//    chancePerCell (deterministic per cell).
//  - Influence at a point is computed by scanning a bounded neighborhood of
//    cells and taking the strongest radial falloff contribution.
//
// All IDs are stable and derived from (universeSeed, cellCoord).

struct GalaxyClustersParams {
  // Coarse grid cell size (ly). <= 0 disables sampling.
  double cellSizeLy{2500.0};

  // Probability [0..1] that a given cell contains a cluster center.
  double chancePerCell{0.08};

  // Typical influence radius (ly).
  double radiusLy{1200.0};

  // Per-cluster radius jitter (multiplicative). 0.35 -> radius in ~[0.65..1.35] * radiusLy.
  double radiusJitter01{0.35};

  // Per-cluster intrinsic-strength jitter (multiplicative). 0.35 -> strength in ~[0.65..1.35].
  double strengthJitter01{0.35};

  // Falloff exponent (>0). Higher values tighten the core (stronger central concentration).
  double falloffPower{2.0};
};

struct GalaxyClusterSample {
  bool hasCluster{false};

  // Stable id for the strongest influencing cluster at this position.
  core::u64 clusterId{0};
  core::u64 clusterSeed{0};

  // Cluster center parameters (valid when hasCluster is true).
  math::Vec3d centerLy{0,0,0};
  double radiusLy{0.0};
  double intrinsicStrength{1.0};

  // Influence at the query point in [0,1]. Includes intrinsicStrength.
  double cluster01{0.0};
};

// Sample the cluster field at `posLy`.
//
// NOTE: The returned cluster corresponds to the *strongest influence* at the
// query point (not necessarily the geometrically nearest center when radii vary).
GalaxyClusterSample sampleGalaxyClusters(core::u64 universeSeed,
                                         const math::Vec3d& posLy,
                                         const GalaxyClustersParams& params = {});

} // namespace stellar::proc
