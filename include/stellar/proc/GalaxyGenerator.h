#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h"
#include "stellar/sim/Faction.h"

#include <vector>

namespace stellar::proc {

struct GalaxyParams {
  double sectorSizeLy{10.0};          // edge length of a sector cube
  double radiusLy{50000.0};           // approximate disc radius
  double thicknessLy{1000.0};         // disc thickness (full)
  double radialScaleLengthLy{15000.0}; // exponential falloff scale
  double baseMeanSystemsPerSector{5.0}; // at r=0
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
