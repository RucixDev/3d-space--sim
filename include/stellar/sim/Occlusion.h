#pragma once

#include "stellar/math/Bvh.h"
#include "stellar/math/Vec3.h"

#include <cstddef>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// OcclusionFieldKm — broad-phase accelerated segment occlusion queries
// -----------------------------------------------------------------------------
//
// Many gameplay systems need a cheap "is there something between A and B" test:
//  - sensors/radar: occlusion attenuation
//  - electronic warfare: jammer line-of-sight
//  - AI: visibility heuristics
//
// This helper builds a BVH over spherical occluders (represented as AABBs for
// broad-phase), then refines hits with an exact segment-sphere test.

struct SphereOccluderKm {
  math::Vec3d centerKm{0.0, 0.0, 0.0};
  double radiusKm{0.0};

  // Contribution strength in [0,1].
  // When the segment intersects this sphere, it is composed into the returned
  // occlusion value with: occ = 1 - (1-occ)*(1-strength01).
  double strength01{1.0};
};

class OcclusionFieldKm {
public:
  void clear();
  bool empty() const;

  std::size_t occluderCount() const { return occluders_.size(); }

  // Build from a list of occluders.
  //
  // Notes:
  // - Invalid/non-finite occluders are filtered out.
  // - leafSize is clamped to [1,64].
  void build(std::vector<SphereOccluderKm> occluders, std::size_t leafSize = 6);

  // Return an occlusion value in [0,1] for the segment [aKm, bKm].
  //
  // `maxSphereTests` caps refinement work in extremely dense scenes.
  // The result is conservative: a lower cap can underestimate occlusion.
  double occlusion01Segment(const math::Vec3d& aKm,
                            const math::Vec3d& bKm,
                            int maxSphereTests = 24) const;

private:
  static bool segmentHitsSphereKm_(const math::Vec3d& aKm,
                                   const math::Vec3d& bKm,
                                   const math::Vec3d& centerKm,
                                   double radiusKm);

  static double clamp01_(double v);

  std::vector<SphereOccluderKm> occluders_;
  math::Bvh3d bvh_;
};

} // namespace stellar::sim
