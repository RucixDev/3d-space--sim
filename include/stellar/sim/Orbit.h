#pragma once

#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h"

namespace stellar::sim {

// Solve Kepler's equation for eccentric anomaly E given mean anomaly M and eccentricity e.
double solveKepler(double meanAnomalyRad, double eccentricity, int iterations = 8);

// Position in orbital plane (AU) at a given time in days.
math::Vec3d orbitPositionAU(const OrbitElements& el, double timeDays);

// Full 3D position (AU) with inclination, node, argument of periapsis applied.
math::Vec3d orbitPosition3DAU(const OrbitElements& el, double timeDays);

// Velocity in orbital plane (AU per day) at a given time in days.
math::Vec3d orbitVelocityAU(const OrbitElements& el, double timeDays);

// Full 3D velocity (AU per day) with inclination, node, argument of periapsis applied.
math::Vec3d orbitVelocity3DAU(const OrbitElements& el, double timeDays);

} // namespace stellar::sim
