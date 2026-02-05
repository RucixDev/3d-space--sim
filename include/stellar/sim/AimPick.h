#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Geometry.h"
#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// AimPick (ray -> sphere candidates)
// -----------------------------------------------------------------------------
//
// A small, deterministic helper for selecting "what the player is aiming at".
// The inputs are intentionally generic:
//  - originKm, dirNormalized define the picking ray in simulation space (kilometers)
//  - candidates are spheres (center/radius), optionally filtered by an aim cone
//  - results are returned sorted from nearest to farthest along the ray
//
// This module lives in `sim` (not `apps/`) so gameplay code can reuse the same
// picking behavior in headless contexts (AI, tests, future server sim).

struct AimPickSphere {
  // User-defined identifier for mapping hits back to higher level objects.
  // For convenience, gameplay can set this to an index into a parallel metadata
  // array.
  std::size_t tag{0};

  math::Vec3d centerKm{0, 0, 0};
  double radiusKm{1.0};

  // Optional aim cone filter.
  // If > -0.5, the candidate is only considered when:
  //   dot(dirNormalized, (center-origin).normalized()) >= minAimCos
  // Use -1.0 to disable.
  double minAimCos{-1.0};
};

struct AimPickHit {
  std::size_t tag{0};
  double tEnterKm{0.0};
  double tExitKm{0.0};
};

inline void aimPickSphereRayHitsKm(std::vector<AimPickHit>& outHits,
                                   const math::Vec3d& originKm,
                                   const math::Vec3d& dirNormalized,
                                   double maxRangeKm,
                                   const AimPickSphere* spheres,
                                   std::size_t sphereCount,
                                   double padKm = 0.0) {
  outHits.clear();
  if (!spheres || sphereCount == 0) return;

  const double maxR = std::max(0.0, maxRangeKm);
  const double pad = std::max(0.0, padKm);

  // Ensure direction is normalized (robust against accidental misuse).
  math::Vec3d d = dirNormalized;
  const double dLen2 = d.lengthSq();
  if (dLen2 < 1.0e-12) return;
  if (std::abs(dLen2 - 1.0) > 1.0e-6) {
    d = d / std::sqrt(dLen2);
  }

  outHits.reserve(sphereCount);

  for (std::size_t i = 0; i < sphereCount; ++i) {
    const AimPickSphere& s = spheres[i];
    const double r = std::max(0.0, s.radiusKm) + pad;

    const math::Vec3d toC = s.centerKm - originKm;
    const double dist2 = toC.lengthSq();

    // Aim cone filter (optional).
    if (s.minAimCos > -0.5) {
      const double dist = (dist2 > 1.0e-12) ? std::sqrt(dist2) : 0.0;
      if (dist > 1.0e-9) {
        const double aim = math::dot(d, toC / dist);
        if (aim < s.minAimCos) continue;
      }
    }

    // Cheap early reject: center projection outside [0, maxRange + radius].
    const double tProj = math::dot(toC, d);
    if (tProj < -r) continue;
    if (tProj > maxR + r) continue;

    double tEnter = 0.0, tExit = 0.0;
    if (!math::raySphereIntersect(originKm, d, s.centerKm, r, tEnter, tExit)) continue;
    if (tEnter > maxR) continue;

    AimPickHit h;
    h.tag = s.tag;
    h.tEnterKm = tEnter;
    h.tExitKm = tExit;
    outHits.push_back(h);
  }

  std::sort(outHits.begin(), outHits.end(), [](const AimPickHit& a, const AimPickHit& b) {
    if (a.tEnterKm != b.tEnterKm) return a.tEnterKm < b.tEnterKm;
    return a.tag < b.tag;
  });
}

inline bool aimPickSphereRayNearestKm(AimPickHit& out,
                                      const math::Vec3d& originKm,
                                      const math::Vec3d& dirNormalized,
                                      double maxRangeKm,
                                      const AimPickSphere* spheres,
                                      std::size_t sphereCount,
                                      double padKm = 0.0) {
  if (!spheres || sphereCount == 0) return false;
  const double maxR = std::max(0.0, maxRangeKm);

  // Ensure direction is normalized.
  math::Vec3d d = dirNormalized;
  const double dLen2 = d.lengthSq();
  if (dLen2 < 1.0e-12) return false;
  if (std::abs(dLen2 - 1.0) > 1.0e-6) d = d / std::sqrt(dLen2);

  const double pad = std::max(0.0, padKm);

  bool any = false;
  double bestT = maxR;
  AimPickHit best{};

  for (std::size_t i = 0; i < sphereCount; ++i) {
    const AimPickSphere& s = spheres[i];
    const double r = std::max(0.0, s.radiusKm) + pad;

    const math::Vec3d toC = s.centerKm - originKm;
    const double dist2 = toC.lengthSq();

    if (s.minAimCos > -0.5) {
      const double dist = (dist2 > 1.0e-12) ? std::sqrt(dist2) : 0.0;
      if (dist > 1.0e-9) {
        const double aim = math::dot(d, toC / dist);
        if (aim < s.minAimCos) continue;
      }
    }

    const double tProj = math::dot(toC, d);
    if (tProj < -r) continue;
    if (tProj > maxR + r) continue;

    double tEnter = 0.0, tExit = 0.0;
    if (!math::raySphereIntersect(originKm, d, s.centerKm, r, tEnter, tExit)) continue;
    if (tEnter > maxR) continue;

    if (tEnter < bestT) {
      any = true;
      bestT = tEnter;
      best.tag = s.tag;
      best.tEnterKm = tEnter;
      best.tExitKm = tExit;
    }
  }

  if (!any) return false;
  out = best;
  return true;
}

} // namespace stellar::sim
