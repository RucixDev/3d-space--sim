#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/System.h"

namespace stellar::sim {

// Returns true if the point (ship position in station-local space) is inside the
// station hull approximation, excluding the docking slot tunnel cutout.
bool insideStationHullExceptSlot(const Station& st, const math::Vec3d& relLocalKm);

// Returns true if the ship meets the docking slot conditions used by the
// prototype: clearance granted, inside tunnel, under speed limit, and aligned.
bool dockingSlotConditions(const Station& st,
                           const math::Vec3d& relLocalKm,
                           const math::Vec3d& shipVelRelLocalKmS,
                           const math::Vec3d& shipForwardLocal,
                           const math::Vec3d& shipUpLocal,
                           bool clearanceGranted);

} // namespace stellar::sim
