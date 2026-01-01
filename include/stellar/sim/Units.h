#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/Orbit.h"
#include "stellar/sim/System.h"

namespace stellar::sim {

// The sim uses kilometers for local/world space while orbital mechanics are expressed in
// astronomical units (AU) and days. These helpers centralize the conversion factors so
// game, tools, and tests can share consistent numbers.

constexpr double kAU_KM = 149597870.7;      // IAU 2012 AU (km) (double precision)
constexpr double kSecondsPerDay = 86400.0;

inline math::Vec3d auToKm(const math::Vec3d& au) {
  return au * kAU_KM;
}

inline math::Vec3d auPerDayToKmS(const math::Vec3d& auPerDay) {
  return auPerDay * (kAU_KM / kSecondsPerDay);
}

inline math::Vec3d stationPosKm(const Station& st, double timeDays) {
  return auToKm(orbitPosition3DAU(st.orbit, timeDays));
}

inline math::Vec3d stationVelKmS(const Station& st, double timeDays) {
  return auPerDayToKmS(orbitVelocity3DAU(st.orbit, timeDays));
}

inline math::Vec3d planetPosKm(const Planet& p, double timeDays) {
  return auToKm(orbitPosition3DAU(p.orbit, timeDays));
}

inline math::Vec3d planetVelKmS(const Planet& p, double timeDays) {
  return auPerDayToKmS(orbitVelocity3DAU(p.orbit, timeDays));
}

} // namespace stellar::sim
