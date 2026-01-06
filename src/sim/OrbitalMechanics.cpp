#include "stellar/sim/OrbitalMechanics.h"

#include "stellar/math/Math.h"

#include <cmath>

namespace stellar::sim {

namespace {

inline double wrapTwoPi(double a) {
  const double twoPi = 2.0 * math::kPi;
  a = std::fmod(a, twoPi);
  if (a < 0.0) a += twoPi;
  return a;
}

} // namespace

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

  const math::Vec3d rHat = relPosKm * (1.0 / r);
  out.radialSpeedKmS = math::dot(relVelKmS, rHat);
  {
    const double vt2 = std::max(0.0, v2 - out.radialSpeedKmS * out.radialSpeedKmS);
    out.tangentialSpeedKmS = std::sqrt(vt2);
  }

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

  const double eTol = 1e-8;
  const bool aValid = (std::abs(out.semiMajorAxisKm) > 1e-9);
  if (aValid && out.semiMajorAxisKm > 0.0 && e < 1.0 - 1e-6) {
    out.type = TwoBodyOrbit::Type::Elliptic;
  } else if (aValid && out.semiMajorAxisKm < 0.0 && e > 1.0 + 1e-6) {
    out.type = TwoBodyOrbit::Type::Hyperbolic;
  } else if (muKm3S2 > 0.0) {
    out.type = TwoBodyOrbit::Type::Parabolic;
  }

  // True anomaly
  double nu = 0.0;
  if (e > eTol) {
    const double cosNu = math::clamp(math::dot(eVec, relPosKm) / (e * r), -1.0, 1.0);
    // sin(nu) = (r·v) * |h| / (mu * e * r)
    const double sinNu = math::clamp((math::dot(relPosKm, relVelKmS) * hMag) / (muKm3S2 * e * r), -1.0, 1.0);
    nu = std::atan2(sinNu, cosNu);
    if (out.type == TwoBodyOrbit::Type::Elliptic) {
      nu = wrapTwoPi(nu);
    }
  } else {
    // Circular-ish: periapsis is undefined; define nu=0 at the current state.
    nu = 0.0;
  }
  out.trueAnomalyRad = nu;

  // Anomalies and timing
  if (out.type == TwoBodyOrbit::Type::Elliptic && out.semiMajorAxisKm > 0.0) {
    const double a = out.semiMajorAxisKm;
    const double n = std::sqrt(muKm3S2 / (a * a * a));
    out.meanMotionRadPerSec = n;
    out.periodSec = (n > 0.0) ? (2.0 * math::kPi / n) : 0.0;

    const double cosNu = std::cos(nu);
    const double sinNu = std::sin(nu);
    const double denom = 1.0 + e * cosNu;

    // cosE = (e + cosNu)/(1 + e cosNu)
    const double cosE = (denom != 0.0) ? ((e + cosNu) / denom) : 1.0;
    // sinE = sqrt(1-e^2) sinNu /(1 + e cosNu)
    const double sinE = (denom != 0.0) ? (std::sqrt(std::max(0.0, 1.0 - e * e)) * sinNu / denom) : 0.0;

    double E = std::atan2(sinE, cosE);
    E = wrapTwoPi(E);
    out.eccentricAnomalyRad = E;

    double M = E - e * sinE;
    M = wrapTwoPi(M);
    out.meanAnomalyRad = M;

    const double tSince = (n > 0.0) ? (M / n) : 0.0;
    out.timeSincePeriapsisSec = tSince;
    out.timeToPeriapsisSec = (out.periodSec > 0.0) ? (out.periodSec - tSince) : 0.0;
    // If we're exactly at periapsis (M≈0) treat time-to-peri as 0.
    if (out.timeToPeriapsisSec >= out.periodSec - 1e-6) {
      out.timeToPeriapsisSec = 0.0;
    }

    const double M_apo = math::kPi;
    double dtApo = 0.0;
    if (n > 0.0) {
      if (M <= M_apo) {
        dtApo = (M_apo - M) / n;
      } else {
        dtApo = ((2.0 * math::kPi - M) + M_apo) / n;
      }
    }
    out.timeToApoapsisSec = dtApo;
  } else if (out.type == TwoBodyOrbit::Type::Hyperbolic && out.semiMajorAxisKm < 0.0) {
    const double a = out.semiMajorAxisKm; // negative
    const double n = std::sqrt(muKm3S2 / ((-a) * (-a) * (-a)));
    out.meanMotionRadPerSec = n;

    // Hyperbolic anomaly H via: tanh(H/2) = sqrt((e-1)/(e+1)) * tan(nu/2)
    const double tanHalfNu = std::tan(nu * 0.5);
    const double fac = std::sqrt((e - 1.0) / (e + 1.0));
    double x = fac * tanHalfNu;
    x = math::clamp(x, -0.999999999999, 0.999999999999);
    const double H = 2.0 * std::atanh(x);
    out.eccentricAnomalyRad = H;
    const double M = e * std::sinh(H) - H;
    out.meanAnomalyRad = M;

    const double tSince = (n > 0.0) ? (M / n) : 0.0;
    out.timeSincePeriapsisSec = tSince;
    out.timeToPeriapsisSec = (tSince < 0.0) ? (-tSince) : 0.0;
    out.timeToApoapsisSec = 0.0;
  }

  return out;
}

} // namespace stellar::sim
