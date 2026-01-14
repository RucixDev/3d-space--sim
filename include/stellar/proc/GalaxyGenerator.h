#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h"
#include "stellar/sim/Faction.h"

#include <vector>

namespace stellar::proc {

struct GalaxyParams {
  double sectorSizeLy{10.0};            // edge length of a sector cube
  double radiusLy{50000.0};             // approximate disc radius
  double thicknessLy{1000.0};           // disc thickness (full)
  double radialScaleLengthLy{15000.0};  // exponential falloff scale
  double baseMeanSystemsPerSector{5.0}; // at r=0

  // --- Optional minimum-separation placement (blue-noise-ish) ---
  //
  // When > 0, the generator enforces an approximate *minimum distance* between
  // neighboring systems using a deterministic, streaming-safe cell + priority
  // scheme (see GalaxyGenerator.cpp).
  //
  // This helps avoid overly tight clumps without requiring multi-sector state
  // and remains coherent across sector boundaries.
  //
  // 0 disables and uses the legacy in-sector Poisson sampling.
  double minSystemSeparationLy{0.0};

  // --- Optional inhomogeneous structure (disabled by default) ---
  // These parameters let you art-direct the galaxy away from a perfectly smooth
  // exponential disc, while keeping generation deterministic.
  //
  // By default, spiralArmCount=0, densityNoiseStrength=0, clusterStrength=0,
  // barStrength=0, ringStrength=0, warpAmplitudeLy=0, and flareStrength=0, which
  // preserves the legacy distribution (and therefore deterministic regression
  // signatures).

  // Log-spiral arm model (0 disables arms).
  int spiralArmCount{0};
  double spiralArmStrength{0.0};        // 0 disables (typical useful range: 0.25 .. 2.5)
  double spiralPitchDeg{12.0};          // pitch angle of the log spiral (degrees)
  double spiralArmWidthDeg{18.0};       // angular width of each arm (degrees)
  double spiralArmPhaseDeg{0.0};        // rotation offset for the whole pattern (degrees)
  double spiralArmNoiseStrength{0.25};  // 0..1 modulation along arms
  double spiralArmNoiseFreq{0.0015};    // noise frequency (1 / ly)

  // Global "clumpiness" density noise (0 disables).
  double densityNoiseStrength{0.0};     // 0..1 recommended
  double densityNoiseFreq{0.0010};      // noise frequency (1 / ly)

  // Sparse star clusters (streaming-safe RBF blobs).
  //
  // When clusterStrength > 0, the generator samples a deterministic coarse-cell
  // field (see GalaxyClusters) and boosts the local mean density:
  //   mean *= (1 + clusterStrength * cluster01)
  //
  // 0 disables clusters and preserves the legacy distribution.
  double clusterStrength{0.0};
  double clusterCellSizeLy{2500.0};
  double clusterChancePerCell{0.08};
  double clusterRadiusLy{1200.0};
  double clusterRadiusJitter{0.35};
  double clusterStrengthJitter{0.35};
  double clusterFalloffPower{2.0};

  // Galactic bar (elongated inner-disc density enhancement).
  // Strength=0 disables.
  double barStrength{0.0};
  double barAngleDeg{25.0};
  double barLengthLy{14000.0}; // semi-major length
  double barWidthLy{4200.0};   // semi-minor width
  double barPower{2.0};        // falloff sharpness (higher = harder edge)

  // Inner ring (density band around a radius).
  // Strength=0 disables.
  double ringStrength{0.0};
  double ringRadiusLy{9000.0};
  double ringWidthLy{1800.0}; // Gaussian sigma-ish
  double ringPower{1.0};      // exponent on the mask (>=1)

  // Disc warp: a sinusoidal vertical displacement that grows with radius.
  // warpAmplitudeLy=0 disables.
  double warpAmplitudeLy{0.0};
  double warpStartRadiusLy{15000.0};
  double warpPower{1.5};
  int warpLobes{2};           // 1=lopsided, 2=classic S-warp
  double warpPhaseDeg{0.0};
  double warpNoiseStrength{0.25};
  double warpNoiseFreq{0.0005}; // 1/ly

  // Disc flare: increases thickness with radius.
  // flareStrength=0 disables.
  double flareStrength{0.0};
  double flarePower{2.0};
};

struct SectorCoord {
  core::i32 x{0}, y{0}, z{0};

  bool operator==(const SectorCoord& o) const { return x==o.x && y==o.y && z==o.z; }
};

struct Sector {
  SectorCoord coord{};
  std::vector<sim::SystemStub> systems{};
};

struct SectorCoordHash {
  std::size_t operator()(const SectorCoord& c) const noexcept;
};

class GalaxyGenerator {
public:
  GalaxyGenerator(core::u64 seed, GalaxyParams params);

  core::u64 seed() const { return seed_; }
  const GalaxyParams& params() const { return params_; }

  SectorCoord sectorOf(const math::Vec3d& posLy) const;

  // Generate the systems contained in a sector (deterministic).
  Sector generateSector(const SectorCoord& coord, const std::vector<sim::Faction>& factions) const;

  // Create a globally unique system id from (sector, localIndex).
  sim::SystemId makeSystemId(const SectorCoord& coord, core::u32 localIndex) const;

private:
  core::u64 seed_{0};
  GalaxyParams params_{};
};

} // namespace stellar::proc
