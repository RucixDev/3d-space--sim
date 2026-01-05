#include "stellar/sim/OrbitalMechanics.h"

#include "stellar/math/Math.h"

#include <cmath>

namespace stellar::sim {

TwoBodyOrbit solveTwoBodyOrbit(const math::Vec3d& relPosKm,
                              const math::Vec3d& relVelKmS,
                              double muKm3S2) {
  TwoBodyOrbit out{};
  out.muKm3S2 = muKm3S2;

  const double r2 = relPosKm.lengthSq();
  const double v2 = relVelKmS.lengthSq();
  if (muKm3S2 <= 0.0 || r2 <= 1e-18) {
    return out;
  }

  const double r = std::sqrt(r2);
  const double v = std::sqrt(v2);
  out.rKm = r;
  out.vKmS = v;

  // Specific angular momentum h = r x v (km^2/s)
  const math::Vec3d h = math::cross(relPosKm, relVelKmS);
  const double hMag = std::sqrt(h.lengthSq());
  out.angularMomentumKm2S = h;
  out.angularMomentumMagKm2S = hMag;

  // Specific orbital energy (km^2/s^2)
  const double eps = 0.5 * v2 - muKm3S2 / r;
  out.specificEnergyKm2S2 = eps;

  // Eccentricity vector: e = (v x h)/mu - r/|r|
  const math::Vec3d eVec = (math::cross(relVelKmS, h) * (1.0 / muKm3S2)) - (relPosKm * (1.0 / r));
  const double e = std::sqrt(eVec.lengthSq());
  out.eccentricityVec = eVec;
  out.eccentricity = e;

  // Semi-major axis a = -mu / (2*eps). For eps~0 (parabolic), treat as unbound.
  if (std::abs(eps) > 1e-12) {
    out.semiMajorAxisKm = -muKm3S2 / (2.0 * eps);
  } else {
    out.semiMajorAxisKm = 0.0;
  }

  // Periapsis / apoapsis. For parabolic or hyperbolic orbits, apoapsis is undefined.
  if (out.semiMajorAxisKm != 0.0) {
    out.periapsisKm = out.semiMajorAxisKm * (1.0 - e);
    if (e < 1.0) {
      out.apoapsisKm = out.semiMajorAxisKm * (1.0 + e);
    } else {
      out.apoapsisKm = 0.0;
    }
  }

  // Orbital period for bound ellipses.
  if (out.semiMajorAxisKm > 0.0 && e < 1.0) {
    out.periodSec = 2.0 * math::kPi * std::sqrt((out.semiMajorAxisKm * out.semiMajorAxisKm * out.semiMajorAxisKm) / muKm3S2);
  }

  return out;
}

} // namespace stellar::sim
