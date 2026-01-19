#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/System.h"

#include <vector>

namespace stellar::proc {

// Procedural “minor body” architecture.
//
// This module synthesizes asteroid belts / debris discs / trojan swarms in a
// deterministic way from:
//   - universeSeed (global)
//   - system.stub.seed (local)
//   - the current planet orbits/types
//
// Design goals:
//   - Deterministic and stable across runs.
//   - Cheap to evaluate for UI/debug visualization.
//   - Feature-rich: resonance gaps (Kirkwood-style) + trojan swarms.
//
// NOTE: These belts are currently a *procedural overlay* and are not persisted.
// They are intended for UI, debug visualization, and future render/gameplay
// hooks.

enum class AsteroidBeltKind : core::u8 {
  MainBelt = 0,    // between inner terrestrial zone and a giant planet
  OuterBelt = 1,   // Kuiper-like belt / debris disc beyond the outermost planet
  TrojanSwarm = 2, // L4/L5 clusters around a dominant planet
  DebrisDisk = 3,  // generic fallback disc when the system layout is sparse
};

const char* asteroidBeltKindName(AsteroidBeltKind k);

// Mean-motion resonance feature.
//
// m:n is the *period ratio* between planet and minor body:
//   - Interior resonance (belt inside planet): planetPeriod / bodyPeriod = m/n, m>n.
//   - Exterior resonance (belt outside planet): bodyPeriod / planetPeriod = m/n, m>n.
//
// strength01:
//   - > 0.0: a density *gap* (dip) with depth strength01
//   - < 0.0: a density *ridge* (boost) with magnitude |strength01|
struct BeltResonanceFeature {
  int m{0};
  int n{0};
  double aAU{0.0};
  double halfWidthAU{0.0};
  double strength01{0.0};
};

struct AsteroidBelt {
  core::u64 id{0};
  AsteroidBeltKind kind{AsteroidBeltKind::MainBelt};

  // Radial span in AU.
  double aInnerAU{2.0};
  double aOuterAU{3.5};

  // Vertical half-thickness in AU (used for sampling points).
  double thicknessAU{0.02};

  // Local orthonormal basis describing the belt plane.
  // basisY is the plane normal.
  math::Vec3d basisX{1.0, 0.0, 0.0};
  math::Vec3d basisY{0.0, 1.0, 0.0};
  math::Vec3d basisZ{0.0, 0.0, 1.0};

  // Index into system.planets. -1 if none.
  int controllingPlanetIndex{-1};

  // Optional azimuthal modulation (density spokes / clumps).
  int mMode{0};
  double mStrength01{0.0};
  double mPhaseRad{0.0};

  // Trojan swarm parameters (only meaningful when kind==TrojanSwarm).
  // The center angle is typically the planet mean anomaly at epoch.
  double trojanCenterThetaRad{0.0};
  double trojanWidthRad{0.45};
  double trojanRadialSigmaAU{0.06};
  double trojanVerticalSigmaAU{0.02};

  std::vector<BeltResonanceFeature> resonances;
};

struct AsteroidBeltPoint {
  math::Vec3d posAU{0.0, 0.0, 0.0};
  double density01{1.0};
};

struct AsteroidBeltPlan {
  std::vector<AsteroidBelt> belts;
};

// Generate deterministic belts for a system.
AsteroidBeltPlan generateAsteroidBelts(core::u64 universeSeed, const sim::StarSystem& system);

// Density in [0,1] for a point at polar coordinates (aAU, thetaRad) in the belt plane.
// (This is 2D: sampling functions apply the vertical falloff separately.)
double asteroidBeltDensity01(const AsteroidBelt& belt, double aAU, double thetaRad);

// Deterministically sample a belt into a quasi-blue-noise point cloud.
//
// The sampler uses a Mitchell-style best-candidate selection with an importance
// weight derived from asteroidBeltDensity01().
std::vector<AsteroidBeltPoint> sampleAsteroidBeltPoints(core::u64 universeSeed,
                                                        const sim::StarSystem& system,
                                                        const AsteroidBelt& belt,
                                                        int pointCount,
                                                        int candidatesPerPoint = 18);

} // namespace stellar::proc
