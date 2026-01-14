#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/proc/GalaxyGenerator.h"

namespace stellar::proc {

// Deterministic galaxy-scale "macro morphology" helpers.
//
// These functions are intentionally lightweight and streaming-safe: they can be
// sampled at any position without storing global state.
//
// The intent is stylized plausibility rather than physical accuracy; the goal
// is to provide art-direction knobs that remain coherent across sector streaming.
//
// Notes:
//  - bar01 and ring01 are masks in [0,1].
//  - densityMul is the combined multiplier from bar+ring only.
//  - warpZLy is the midplane offset in ly.
//  - thicknessHalfLy is the local half-thickness in ly (includes flare).

struct GalaxyMorphologySample {
  double bar01{0.0};
  double ring01{0.0};

  double densityMul{1.0};

  double warpZLy{0.0};
  double thicknessHalfLy{0.0};
};

// Elliptical bar mask in [0,1] (0 when barStrength==0 or outside the bar).
// The bar is centered on the galactic origin.
//
// The mask is compact-support: it is exactly 0 outside the ellipse.
// Inside the ellipse, it falls off as (1-q)^barPower where
// q = sqrt((x'/barLength)^2 + (y'/barWidth)^2).
//
// x',y' are (x,y) rotated by barAngleDeg.
//
// This is a density *shape* helper; the barStrength multiplier is applied by
// the caller.

double galaxyBarMask01(const GalaxyParams& gp, double xLy, double yLy);

// Radial ring mask in [0,1] (0 when ringStrength==0).
// ringWidthLy is treated like a Gaussian sigma.
double galaxyRingMask01(const GalaxyParams& gp, double rxyLy);

// Local disc half-thickness (ly), including optional flare.
double galaxyThicknessHalfLy(const GalaxyParams& gp, double rxyLy);

// Local disc midplane offset (ly), including optional warp.
double galaxyWarpZLy(core::u64 seed, const GalaxyParams& gp, double xLy, double yLy);

GalaxyMorphologySample sampleGalaxyMorphology(core::u64 seed, const GalaxyParams& gp, const math::Vec3d& posLy);

// Convenience: density multiplier from bar+ring only.
double galaxyMorphologyDensityMul(core::u64 seed, const GalaxyParams& gp, const math::Vec3d& posLy);

} // namespace stellar::proc
