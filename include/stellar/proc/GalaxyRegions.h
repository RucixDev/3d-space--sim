#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

#include <string>

namespace stellar::proc {

// A coarse, galaxy-scale "region" (biome) classification built from a 3D
// cellular/Worley-style Voronoi partition.
//
// Design goals:
//  - Deterministic: depends only on (universeSeed, position, cellSizeLy)
//  - Spatially coherent: nearby points usually share the same region
//  - Cheap: constant-time sampling (checks a fixed neighborhood of cells)
//
// NOTE: This is currently a gameplay/UX-facing label layer (lab window + tools).
// Future rounds can plug region traits into trade/security/mission generation.

enum class GalaxyRegionKind : core::u8 {
  Core = 0,
  InnerDisc = 1,
  OuterRim = 2,
  Nebula = 3,
  Cluster = 4,
  Rift = 5,
  Count
};

const char* galaxyRegionKindName(GalaxyRegionKind k);

struct GalaxyRegionSample {
  // Stable id for the region (derived from the winning Voronoi feature point).
  core::u64 regionId{0};
  core::u64 regionSeed{0};

  GalaxyRegionKind kind{GalaxyRegionKind::OuterRim};

  // Voronoi feature-point position.
  math::Vec3d centerLy{0,0,0};

  // Distance from query position to the region's feature-point.
  double distanceToCenterLy{0.0};

  // A heuristic cell "radius" estimate (not exact; useful for UI/debug).
  double approxRadiusLy{0.0};

  // 0 at region center, ~1 near a boundary between regions.
  double edge01{0.0};

  // Deterministic, stylized region name (e.g. "Alneris Expanse").
  std::string name;
};

// Sample the region at `posLy`.
//
// cellSizeLy controls the typical size of regions (larger -> bigger regions).
// Practical ranges:
//   250..2500 ly  (depending on your galaxy scale)
GalaxyRegionSample sampleGalaxyRegion(core::u64 universeSeed,
                                      const math::Vec3d& posLy,
                                      double cellSizeLy = 900.0);

} // namespace stellar::proc
