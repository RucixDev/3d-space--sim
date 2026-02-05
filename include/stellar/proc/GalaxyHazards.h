#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

namespace stellar::proc {

// -----------------------------------------------------------------------------
// GalaxyHazards
// -----------------------------------------------------------------------------
//
// A deterministic, galaxy-coherent “space weather” / hazard field.
//
// Outputs smooth values in [0,1] that can be used to bias lane generation,
// travel risk, sensor visibility, mission seeding, etc.
//
// The field also supports a simple time parameter (timeDays) that linearly drifts
// the noise domain in a deterministic direction derived from the universe seed.

enum class GalaxyHazardKind : core::u8 {
  Clear = 0,
  Nebula = 1,
  Storm = 2,
  Rift = 3,
};

const char* galaxyHazardKindName(GalaxyHazardKind k);

struct GalaxyHazardsParams {
  // Low-frequency “nebula density” scale (ly).
  double nebulaScaleLy{650.0};

  // Higher-frequency “storm” scale (ly).
  double stormScaleLy{180.0};

  // Octaves for each channel.
  int nebulaOctaves{4};
  int stormOctaves{4};

  // Time (days). The field is evaluated at:
  //   posLy + driftDir * driftLyPerDay * timeDays
  double timeDays{0.0};

  // Drift speed (ly/day).
  double driftLyPerDay{0.65};

  // If > 0, sample the GalaxyRegions layer to bias hazards to match region kinds.
  // Set <= 0 to disable region modulation.
  double regionCellSizeLy{900.0};

  // Strength of region kind bias.
  //  0 => ignore region kinds
  //  1 => strongly match region flavor
  double regionInfluence{0.35};

  // Global post-strength (applied to the combined hazard). Values >1 are allowed,
  // but the final hazard is clamped to [0,1].
  double strength{1.0};
};

struct GalaxyHazardSample {
  GalaxyHazardKind kind{GalaxyHazardKind::Clear};

  // Component channels.
  double nebula01{0.0}; // sensor occlusion / cover
  double storm01{0.0};  // navigation disruption

  // Combined scalar in [0,1].
  double hazard01{0.0};

  // Convenience aliases.
  double sensorOcclusion01{0.0};
  double navDisruption01{0.0};
};

// Deterministically sample the hazard field at a position.
GalaxyHazardSample sampleGalaxyHazards(core::u64 universeSeed,
                                       const math::Vec3d& posLy,
                                       const GalaxyHazardsParams& params = {});

// Convenience: average hazard sampled along a straight segment.
//
// Sampling uses midpoints (avoids hitting exact endpoints). A small sample
// count (3-7) is usually enough for smooth galaxy-scale noise.
double sampleGalaxyHazardAvgOnSegment(core::u64 universeSeed,
                                      const math::Vec3d& aLy,
                                      const math::Vec3d& bLy,
                                      const GalaxyHazardsParams& params = {},
                                      int samples = 3);

// Convenience: average nav disruption (storm01) sampled along a straight segment.
//
// This is useful for hazard-aware route planning and for gameplay systems that
// care specifically about navigation disruption rather than the combined hazard.
//
// Sampling uses midpoints (avoids hitting exact endpoints). A small sample
// count (3-7) is usually enough for smooth galaxy-scale noise.
double sampleGalaxyNavDisruptionAvgOnSegment(core::u64 universeSeed,
                                             const math::Vec3d& aLy,
                                             const math::Vec3d& bLy,
                                             const GalaxyHazardsParams& params = {},
                                             int samples = 3);

// Convenience: average sensor occlusion (nebula01) sampled along a straight segment.
//
// This can be used for sensor range attenuation, stealth/cover models, and
// mission seeding.
double sampleGalaxySensorOcclusionAvgOnSegment(core::u64 universeSeed,
                                               const math::Vec3d& aLy,
                                               const math::Vec3d& bLy,
                                               const GalaxyHazardsParams& params = {},
                                               int samples = 3);

} // namespace stellar::proc
