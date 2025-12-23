#include "stellar/sim/Orbit.h"

#include "stellar/math/Math.h"

#include <cmath>

namespace stellar::sim {

static double normalizeAngle(double a) {
  const double twoPi = 2.0 * stellar::math::kPi;
  a = std::fmod(a, twoPi);
  if (a < 0) a += twoPi;
  return a;
}

double solveKepler(double meanAnomalyRad, double e, int iterations) {
  const double M = normalizeAngle(meanAnomalyRad);
  double E = (e < 0.8) ? M : stellar::math::kPi;

  for (int i = 0; i < iterations; ++i) {
    const double f = E - e * std::sin(E) - M;
    const double fp = 1.0 - e * std::cos(E);
    E = E - f / fp;
  }
  return E;
}

math::Vec3d orbitPositionAU(const OrbitElements& el, double timeDays) {
  const double n = (2.0 * stellar::math::kPi) / el.periodDays; // mean motion rad/day
  const double M = el.meanAnomalyAtEpochRad + n * (timeDays - el.epochDays);
  const double E = solveKepler(M, el.eccentricity);

  const double cosE = std::cos(E);
  const double sinE = std::sin(E);

  const double a = el.semiMajorAxisAU;
  const double r = a * (1.0 - el.eccentricity * cosE);

  // True anomaly
  const double v = std::atan2(std::sqrt(1.0 - el.eccentricity*el.eccentricity) * sinE,
                              cosE - el.eccentricity);

  return { r * std::cos(v), r * std::sin(v), 0.0 };
}

math::Vec3d orbitPosition3DAU(const OrbitElements& el, double timeDays) {
  const math::Vec3d p = orbitPositionAU(el, timeDays);

  const double cosO = std::cos(el.ascendingNodeRad);
  const double sinO = std::sin(el.ascendingNodeRad);
  const double cosi = std::cos(el.inclinationRad);
  const double sini = std::sin(el.inclinationRad);
  const double cosw = std::cos(el.argPeriapsisRad);
  const double sinw = std::sin(el.argPeriapsisRad);

  // Rotation Rz(O) * Rx(i) * Rz(w)
  const double x = p.x;
  const double y = p.y;

  const double x1 = cosw*x - sinw*y;
  const double y1 = sinw*x + cosw*y;

  const double x2 = x1;
  const double y2 = cosi*y1;
  const double z2 = sini*y1;

  const double x3 = cosO*x2 - sinO*y2;
  const double y3 = sinO*x2 + cosO*y2;
  const double z3 = z2;

  return {x3, y3, z3};
}

} // namespace stellar::sim
