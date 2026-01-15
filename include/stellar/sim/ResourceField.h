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

// Geometric layout for a generated resource field.
//
// This is primarily a *render/gameplay hint*: it describes how asteroids are
// distributed around ResourceFieldSite::posKm.
//
// Keep numeric values stable (they may be implicitly persisted via deterministic
// world id schemes and save-game depletion tables).
enum class ResourceFieldLayout : core::u8 {
  // Roughly spherical / ellipsoidal cluster.
  Cluster = 0,
  // Torus/ring (optionally an arc segment).
  Torus = 1,
  // Thin sheet/disc (useful for icy debris planes).
  Sheet = 2,
};

const char* resourceFieldKindName(ResourceFieldKind k);
const char* resourceFieldLayoutName(ResourceFieldLayout l);

// Sub-structure features used to make fields feel less uniform.
//
// These are *not* persisted; they are generated deterministically from the
// field id and returned in ResourceFieldPlan for UI/debug and for any gameplay
// systems that want to bias behavior toward "hot" regions.
enum class ResourceFieldFeatureKind : core::u8 {
  Hotspot = 0,   // localized density increase (belt clumps, cluster pockets)
  Gap = 1,       // localized density decrease (belt gaps)
  Streak = 2,    // filament/stream inside a sheet
  Spokes = 3,    // subtle periodic modulation along a belt
};

const char* resourceFieldFeatureKindName(ResourceFieldFeatureKind k);

// Generic feature descriptor.
//
// Interpretation notes:
//  - For Torus features: angleRad is the belt angle (field-local X/Z plane).
//    width is an angular half-width in radians.
//  - For Sheet streaks: angleRad is the streak direction in the X/Z plane.
//    width is a half-width in km.
//  - For Cluster hotspots: localPos is a normalized local-space center (roughly
//    in [-1,1]^3) and width is a gaussian sigma in normalized units.
struct ResourceFieldFeature {
  core::u64 fieldId{0};
  ResourceFieldFeatureKind kind{ResourceFieldFeatureKind::Hotspot};

  double angleRad{0.0};
  double width{0.0};
  double strength01{0.0};

  // Optional extra parameter (e.g., spokes frequency).
  double param{0.0};

  // Optional local-space position (used by Cluster hotspots).
  math::Vec3d localPos{0.0, 0.0, 0.0};
};

// Deterministically generated "resource field" signal/site.
//
// The id is stable across runs and can be used to derive asteroid ids.
struct ResourceFieldSite {
  core::u64 id{0};
  ResourceFieldKind kind{ResourceFieldKind::OreBelt};

  // Spatial layout of asteroids in this field.
  ResourceFieldLayout layout{ResourceFieldLayout::Cluster};

  // Position in system-space (km). Caller supplies the anchor position that the
  // site is placed relative to (usually a station).
  math::Vec3d posKm{0, 0, 0};

  // Local orthonormal basis for the layout.
  //
  // For Torus/Sheet: basisY is the plane normal, basisX/basisZ span the plane.
  // For Cluster: the basis is arbitrary but deterministic (can be used for
  // ellipsoid-like shaping).
  math::Vec3d basisX{1, 0, 0};
  math::Vec3d basisY{0, 1, 0};
  math::Vec3d basisZ{0, 0, 1};

  // Layout parameters (km). Interpretation depends on `layout`:
  //  - Cluster: majorRadiusKm = cluster radius. minorRadiusKm may be used as an
  //            optional "tightness"/secondary axis scale.
  //  - Torus:   majorRadiusKm = ring radius (centerline), minorRadiusKm = tube radius.
  //            arcRad < 2π produces a belt arc.
  //  - Sheet:   majorRadiusKm = in-plane half-width (disc radius), minorRadiusKm =
  //            half-thickness along basisY.
  double majorRadiusKm{60000.0};
  double minorRadiusKm{14000.0};
  double arcRad{6.283185307179586};      // 2π by default
  double arcCenterRad{0.0};              // only used when arcRad < 2π

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

  // Local density/structure weight in [0,1].
  //
  // This is useful for UI/debug visualization and for gameplay systems that
  // want to bias scanning/mining behavior toward pockets/clumps.
  double density01{1.0};

  // Baseline yield capacity before depletion/persistence overrides.
  double baseUnits{120.0};
};

struct ResourceFieldPlan {
  std::vector<ResourceFieldSite> fields;
  std::vector<ResourceFieldFeature> features;
  std::vector<ResourceAsteroid> asteroids;
};

// Generate persistent resource fields for a system.
//
// Design goals:
//  - Stable deterministic IDs (tagged with kDeterministicWorldIdBit)
//  - Positions are stable relative to the caller-supplied anchor position
//    (so they "move" with orbiting stations if the anchor moves)
//  - Optional preferred plane normal lets callers align belts/sheets to a
//    meaningful orbital plane (e.g., the anchor station's orbit normal)
//  - Yield mixes are stable and suitable for scan/HUD readouts
//
// NOTE: This function does not apply depletion persistence; callers should apply
// any saved remaining-units overrides by asteroid id.
ResourceFieldPlan generateResourceFields(core::u64 universeSeed,
                                        SystemId systemId,
                                        const math::Vec3d& anchorPosKm,
                                        double anchorCommsRangeKm,
                                        int fieldCount = 3,
                                        math::Vec3d preferredPlaneNormalKm = {0.0, 0.0, 0.0});

// Helper: return all features that belong to a given field id.
std::vector<ResourceFieldFeature> filterFeaturesForField(const std::vector<ResourceFieldFeature>& features,
                                                         core::u64 fieldId);

// Evaluate the deterministic "structure density" for a position inside a
// resource field.
//
// This is the same density function used by the generator to bias Poisson-ish
// placement and to scale yield.
//
// Returns a value in [0,1] where 0 is a "void" and 1 is a hotspot.
double resourceFieldDensity01(const ResourceFieldSite& site,
                              const std::vector<ResourceFieldFeature>& features,
                              const math::Vec3d& worldPosKm);

// Helper: return all asteroids that belong to a given field id.
std::vector<ResourceAsteroid> filterAsteroidsForField(const std::vector<ResourceAsteroid>& asteroids,
                                                      core::u64 fieldId);

} // namespace stellar::sim
