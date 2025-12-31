#include "stellar/sim/Docking.h"

#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

bool insideStationHullExceptSlot(const Station& st, const math::Vec3d& relLocalKm) {
  // Local station hull approximation: box (wx,wy,wz), with a rectangular slot tunnel cut out.
  const double wx = st.radiusKm * 0.70;
  const double wy = st.radiusKm * 0.70;
  const double wz = st.radiusKm * 1.10;

  const bool insideBox = (std::abs(relLocalKm.x) < wx) &&
                         (std::abs(relLocalKm.y) < wy) &&
                         (std::abs(relLocalKm.z) < wz);
  if (!insideBox) return false;

  // Slot tunnel cutout near +Z face (entrance at +wz)
  const double slotHalfW = st.slotWidthKm * 0.5;
  const double slotHalfH = st.slotHeightKm * 0.5;

  const double zEntrance = wz;
  const double zMin = zEntrance - st.slotDepthKm;

  const bool insideTunnel =
    (std::abs(relLocalKm.x) < slotHalfW) &&
    (std::abs(relLocalKm.y) < slotHalfH) &&
    (relLocalKm.z <= zEntrance) &&
    (relLocalKm.z >= zMin);

  return !insideTunnel;
}

bool dockingSlotConditions(const Station& st,
                          const math::Vec3d& relLocalKm,
                          const math::Vec3d& shipVelRelLocalKmS,
                          const math::Vec3d& shipForwardLocal,
                          const math::Vec3d& shipUpLocal,
                          bool clearanceGranted) {
  if (!clearanceGranted) return false;

  const double wx = st.radiusKm * 0.70;
  const double wy = st.radiusKm * 0.70;
  const double wz = st.radiusKm * 1.10;
  (void)wx;
  (void)wy;

  const double zEntrance = wz;
  const double zMin = zEntrance - st.slotDepthKm;

  // Must be inside tunnel volume (i.e. have flown into the mail-slot)
  const double slotHalfW = st.slotWidthKm * 0.5;
  const double slotHalfH = st.slotHeightKm * 0.5;

  const bool insideTunnel =
    (std::abs(relLocalKm.x) < slotHalfW) &&
    (std::abs(relLocalKm.y) < slotHalfH) &&
    (relLocalKm.z <= zEntrance - 0.05 * st.radiusKm) &&
    (relLocalKm.z >= zMin + 0.10 * st.radiusKm);

  if (!insideTunnel) return false;

  const double relSpeed = shipVelRelLocalKmS.length();
  if (relSpeed > st.maxApproachSpeedKmS) return false;

  // Orientation: ship forward should point into the station (-Z local),
  // and ship up should be roughly +Y local (roll alignment).
  const double fwdAlign = math::dot(shipForwardLocal.normalized(), math::Vec3d{0,0,-1});
  const double upAlign  = math::dot(shipUpLocal.normalized(),      math::Vec3d{0,1,0});
  return (fwdAlign > 0.92 && upAlign > 0.70);
}

} // namespace stellar::sim
