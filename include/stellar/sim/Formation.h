#pragma once

#include "stellar/core/Random.h"
#include "stellar/core/Types.h"
#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Formation helpers (deterministic, headless)
// -----------------------------------------------------------------------------
//
// These helpers provide simple "good looking" spatial arrangements around a
// leader ship (convoy escorting, squad travel, etc.).
//
// Design goals:
//  - deterministic (seeded), so behaviour is stable across runs/saves
//  - headless (math-only) so it can be used in tests and non-render builds
//  - flexible enough to be used as a building block by gameplay code
// -----------------------------------------------------------------------------

enum class FormationPattern : core::u8 {
  Trail = 0,   // line/stack behind the leader
  Wedge = 1,   // V shape behind the leader
  Diamond = 2, // compact 3D "diamond" behind the leader
  Ring = 3     // ring around the leader (behind plane)
};

// Tunable parameters controlling formation shape.
struct FormationParams {
  // Overall scale (roughly the desired separation from the leader).
  double radiusKm{36000.0};

  // Lateral/vertical scales relative to radius.
  double lateralScale{0.45};
  double verticalScale{0.22};

  // How far behind the leader to place the formation (in multiples of radius).
  double behindScale{1.0};

  // Row spacing for Trail/Wedge (in multiples of radius).
  double rowSpacingScale{0.65};

  // Small deterministic jitter to avoid perfect symmetry (km).
  double jitterKm{800.0};
};

// A tiny orthonormal basis (right/up/forward).
struct Basis {
  math::Vec3d right{1, 0, 0};
  math::Vec3d up{0, 1, 0};
  math::Vec3d fwd{0, 0, 1};
};

// Build a stable orthonormal basis from a desired forward vector and an up hint.
//
// If forward or up are degenerate/parallel, a safe fallback is chosen.
inline Basis makeBasisFromForward(const math::Vec3d& forwardWorld,
                                 const math::Vec3d& upHintWorld = {0, 1, 0}) {
  math::Vec3d f = forwardWorld;
  if (f.lengthSq() <= 1e-12) f = {0, 0, 1};
  f = f.normalized();

  math::Vec3d up = upHintWorld;
  if (up.lengthSq() <= 1e-12) up = {0, 1, 0};
  up = up.normalized();

  math::Vec3d r = math::cross(up, f);
  if (r.lengthSq() <= 1e-12) {
    // Up is parallel to forward; choose a fallback up.
    up = (std::abs(f.y) < 0.9) ? math::Vec3d{0, 1, 0} : math::Vec3d{1, 0, 0};
    r = math::cross(up, f);
  }
  r = r.normalized();

  // Re-orthogonalize up.
  up = math::cross(f, r).normalized();
  return {r, up, f};
}

inline FormationPattern chooseEscortFormationPattern(int slotCount) {
  if (slotCount <= 1) return FormationPattern::Trail;
  if (slotCount <= 3) return FormationPattern::Wedge;
  if (slotCount <= 6) return FormationPattern::Diamond;
  return FormationPattern::Ring;
}

// Compute an offset in leader-local space (x=right, y=up, z=forward).
//
// The returned vector is in *km* and is designed so that, for most patterns,
// offsets are roughly "behind" the leader (negative Z) so escorts do not
// occlude the leader's forward path.
inline math::Vec3d formationOffsetLocalKm(FormationPattern pat,
                                         int slotIndex,
                                         int slotCount,
                                         const FormationParams& params,
                                         core::u64 seed) {
  slotCount = std::max(1, slotCount);
  slotIndex = std::clamp(slotIndex, 0, slotCount - 1);

  const double R = std::max(0.0, params.radiusKm);
  const double lat = R * params.lateralScale;
  const double vert = R * params.verticalScale;
  const double back = -R * params.behindScale;
  const double rowStep = R * params.rowSpacingScale;

  math::Vec3d off{0, 0, 0};

  // Small deterministic jitter (breaks perfect symmetry without causing drift).
  auto jitter = [&](core::u64 s) -> math::Vec3d {
    if (params.jitterKm <= 1e-9) return {0, 0, 0};
    core::SplitMix64 rng(s);
    auto u = [&]() { return rng.nextDouble() * 2.0 - 1.0; };
    return {
      u() * params.jitterKm,
      u() * params.jitterKm * 0.65,
      u() * params.jitterKm,
    };
  };

  switch (pat) {
    case FormationPattern::Trail: {
      // A line behind the leader with slight alternating offsets.
      const int i = slotIndex;
      const double side = ((i & 1) != 0) ? 1.0 : -1.0;
      const double upSign = ((i & 2) != 0) ? 1.0 : -1.0;
      off.x = side * lat * 0.12;
      off.y = upSign * vert * 0.10;
      off.z = back - rowStep * static_cast<double>(i);
    } break;

    case FormationPattern::Wedge: {
      // A simple V behind the leader.
      const int row = slotIndex / 2 + 1;
      const double side = ((slotIndex & 1) != 0) ? 1.0 : -1.0;
      const double widen = 0.75 + 0.25 * std::clamp((double)slotCount / 6.0, 0.0, 1.0);
      off.x = side * lat * static_cast<double>(row) * widen;
      off.y = (((slotIndex / 4) & 1) != 0) ? -vert * 0.18 : vert * 0.18;
      off.z = back - rowStep * static_cast<double>(row - 1);
    } break;

    case FormationPattern::Diamond: {
      // Compact 3D pattern for small groups; spill extras into a small ring.
      const int i = slotIndex;
      if (slotCount <= 1) {
        off = {0, 0, back};
        break;
      }

      if (i == 0) {
        off = {-lat * 0.65, 0.0, back};
      } else if (i == 1) {
        off = {lat * 0.65, 0.0, back};
      } else if (i == 2) {
        off = {0.0, vert * 0.80, back * 0.85};
      } else if (i == 3) {
        off = {0.0, -vert * 0.80, back * 0.85};
      } else {
        const int extra = i - 4;
        const int extraN = std::max(1, slotCount - 4);
        const double t = (extraN > 0) ? (static_cast<double>(extra) / static_cast<double>(extraN)) : 0.0;
        const double ang = 2.0 * math::kPi * t;
        off.x = std::cos(ang) * lat * 0.85;
        off.y = std::sin(ang) * vert * 0.85;
        off.z = back - rowStep * 0.35;
      }
    } break;

    case FormationPattern::Ring: {
      // Ring on the plane perpendicular to forward, placed slightly behind.
      const double t = (slotCount > 0) ? (static_cast<double>(slotIndex) / static_cast<double>(slotCount)) : 0.0;
      const double ang = 2.0 * math::kPi * t;
      off.x = std::cos(ang) * lat;
      off.y = std::sin(ang) * vert;
      off.z = back;
    } break;
  }

  off += jitter(core::hashCombine(seed, static_cast<core::u64>(slotIndex)));
  return off;
}

// Transform a local formation offset into world space.
inline math::Vec3d formationTargetWorldKm(const math::Vec3d& leaderPosKm,
                                         const math::Vec3d& rightWorld,
                                         const math::Vec3d& upWorld,
                                         const math::Vec3d& forwardWorld,
                                         const math::Vec3d& offsetLocalKm) {
  return leaderPosKm + rightWorld * offsetLocalKm.x + upWorld * offsetLocalKm.y + forwardWorld * offsetLocalKm.z;
}

inline math::Vec3d formationTargetWorldKm(const math::Vec3d& leaderPosKm,
                                         const Basis& basis,
                                         const math::Vec3d& offsetLocalKm) {
  return formationTargetWorldKm(leaderPosKm, basis.right, basis.up, basis.fwd, offsetLocalKm);
}

} // namespace stellar::sim
