#pragma once

#include "stellar/econ/Commodity.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h"
#include "stellar/sim/WorldIds.h"

#include <cstddef>
#include <vector>

namespace stellar::sim {

// High-level classification for a generated resource field.
//
// This is primarily used for UI/scan readouts, but it also drives yield distribution.
// Keep numeric values stable (persisted implicitly via deterministic generation).
enum class ResourceFieldKind : core::u8 {
  OreBelt = 0,
  MetalPocket = 1,
  IceField = 2,
};

const char* resourceFieldKindName(ResourceFieldKind k);

// Deterministically generated "resource field" signal/site.
//
// The id is stable across runs and can be used to derive asteroid ids.
struct ResourceFieldSite {
  core::u64 id{0};
  ResourceFieldKind kind{ResourceFieldKind::OreBelt};

  // Position in system-space (km). Caller supplies the anchor position that the
  // site is placed relative to (usually a station).
  math::Vec3d posKm{0, 0, 0};

  // Richness multiplier applied to asteroid yield (roughly ~[0.75,1.40]).
  double richness{1.0};

  // Composition summary used for scan/HUD readouts.
  econ::CommodityId primary{econ::CommodityId::Ore};
  econ::CommodityId secondary{econ::CommodityId::Metals};
  double secondaryChance{0.0}; // probability that an asteroid yields `secondary`
};

// Deterministically generated asteroid node inside a resource field.
struct ResourceAsteroid {
  core::u64 id{0};
  core::u64 fieldId{0};

  math::Vec3d posKm{0, 0, 0};
  double radiusKm{5000.0};

  econ::CommodityId yield{econ::CommodityId::Ore};

  // Baseline yield capacity before depletion/persistence overrides.
  double baseUnits{120.0};
};

struct ResourceFieldPlan {
  std::vector<ResourceFieldSite> fields;
  std::vector<ResourceAsteroid> asteroids;
};

// Generate persistent resource fields for a system.
//
// Design goals:
//  - Stable deterministic IDs (tagged with kDeterministicWorldIdBit)
//  - Positions are stable relative to the caller-supplied anchor position
//    (so they "move" with orbiting stations if the anchor moves)
//  - Yield mixes are stable and suitable for scan/HUD readouts
//
// NOTE: This function does not apply depletion persistence; callers should apply
// any saved remaining-units overrides by asteroid id.
ResourceFieldPlan generateResourceFields(core::u64 universeSeed,
                                        SystemId systemId,
                                        const math::Vec3d& anchorPosKm,
                                        double anchorCommsRangeKm,
                                        int fieldCount = 3);

// Helper: return all asteroids that belong to a given field id.
std::vector<ResourceAsteroid> filterAsteroidsForField(const std::vector<ResourceAsteroid>& asteroids,
                                                      core::u64 fieldId);

} // namespace stellar::sim
