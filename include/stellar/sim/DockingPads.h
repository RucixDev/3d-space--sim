#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/System.h"

namespace stellar::sim {

struct DockingPad {
  // 1-based index (pad number shown to the player).
  core::u16 index{0};

  // Docking pose in station-local coordinates.
  // Station local axes follow the Docking model:
  //  +X right, +Y up, +Z outward from the slot.
  math::Vec3d localPosKm{0,0,0};
  math::Quatd localOrient{1,0,0,0};

  // Approximate radius used for UI / occupancy heuristics.
  double radiusKm{0.0};
};

// Deterministic number of pads available inside the hangar.
int stationDockingPadCount(const Station& st);

// Deterministically generate a docking pad pose.
// If index1Based is out of range, a fallback pad is returned.
DockingPad stationDockingPad(core::u64 universeSeed, const Station& st, core::u16 index1Based);

// Find the nearest pad to a point in station-local space.
core::u16 nearestDockingPadIndex(core::u64 universeSeed, const Station& st, const math::Vec3d& relLocalKm);

} // namespace stellar::sim
