#include "stellar/sim/OrbitalMechanics.h"

#include "stellar/math/Math.h"
#include "stellar/math/Quat.h"

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


// -----------------------------------------------------------------------------
// Classical orbital elements
// -----------------------------------------------------------------------------

namespace {

inline double safeAcos(double x) {
  return std::acos(math::clamp(x, -1.0, 1.0));
}

inline bool isFinite(const math::Vec3d& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

inline bool isFinite(double v) { return std::isfinite(v); }

// Eccentric anomaly from true anomaly (elliptic).
inline double eccentricAnomalyFromTrueAnomaly(double e, double nu) {
  // tan(E/2) = sqrt((1-e)/(1+e)) * tan(nu/2)
  const double t = std::tan(nu * 0.5);
  const double fac = std::sqrt(std::max(0.0, (1.0 - e) / (1.0 + e)));
  const double x = fac * t;
  double E = 2.0 * std::atan(x);

  // Keep in [0,2pi).
  return wrapTwoPi(E);
}

inline double meanAnomalyFromTrueAnomaly(double e, double nu) {
  const double E = eccentricAnomalyFromTrueAnomaly(e, nu);
  return wrapTwoPi(E - e * std::sin(E));
}

} // namespace

ClassicalOrbitalElements solveClassicalOrbitElements(const math::Vec3d& relPosKm,
                                                     const math::Vec3d& relVelKmS,
                                                     double muKm3S2) {
  ClassicalOrbitalElements out{};
  out.muKm3S2 = muKm3S2;

  if (!isFinite(relPosKm) || !isFinite(relVelKmS) || !isFinite(muKm3S2) || muKm3S2 <= 0.0) {
    return out;
  }

  const double r2 = relPosKm.lengthSq();
  const double v2 = relVelKmS.lengthSq();
  if (!(r2 > 1e-18) || !(v2 >= 0.0)) {
    return out;
  }

  const double r = std::sqrt(r2);
  const double v = std::sqrt(v2);

  // Angular momentum
  const math::Vec3d h = math::cross(relPosKm, relVelKmS);
  const double h2 = h.lengthSq();
  const double hMag = std::sqrt(h2);
  if (!(hMag > 1e-18) || !isFinite(hMag)) {
    // Radial trajectory; classical elements are ill-defined.
    return out;
  }
  out.angularMomentumHat = h * (1.0 / hMag);

  // Eccentricity vector.
  const math::Vec3d eVec = (math::cross(relVelKmS, h) * (1.0 / muKm3S2)) - (relPosKm * (1.0 / r));
  const double e = std::sqrt(std::max(0.0, eVec.lengthSq()));
  out.eccentricityVec = eVec;
  out.eccentricity = e;

  // Energy and semi-major axis.
  const double eps = 0.5 * v2 - muKm3S2 / r;
  if (std::abs(eps) > 1e-12) {
    out.semiMajorAxisKm = -muKm3S2 / (2.0 * eps);
  } else {
    out.semiMajorAxisKm = 0.0;
  }

  // p = h^2 / mu (valid for all conics)
  out.parameterKm = h2 / muKm3S2;

  // Classify.
  if (out.semiMajorAxisKm > 0.0 && e < 1.0 - 1e-6) {
    out.type = TwoBodyOrbit::Type::Elliptic;
  } else if (out.semiMajorAxisKm < 0.0 && e > 1.0 + 1e-6) {
    out.type = TwoBodyOrbit::Type::Hyperbolic;
  } else {
    out.type = TwoBodyOrbit::Type::Parabolic;
  }

  // Inclination
  out.inclinationRad = safeAcos(out.angularMomentumHat.z);

  // Node vector (pointing to ascending node)
  const math::Vec3d kHat{0, 0, 1};
  const math::Vec3d nVec = math::cross(kHat, h);
  const double nMag = std::sqrt(std::max(0.0, nVec.lengthSq()));

  constexpr double kEquatorialTol = 1e-10;
  constexpr double kCircularTol = 1e-10;

  out.equatorial = (nMag < kEquatorialTol);
  out.circular = (e < kCircularTol);

  // RAAN Ω
  if (!out.equatorial) {
    double Omega = safeAcos(nVec.x / nMag);
    if (nVec.y < 0.0) Omega = 2.0 * math::kPi - Omega;
    out.raanRad = wrapTwoPi(Omega);
  } else {
    out.raanRad = 0.0;
  }

  // Argument of periapsis ω and true anomaly ν (or argument of latitude u).
  if (!out.circular) {
    // Non-circular: periapsis direction defined by eVec.
    if (!out.equatorial) {
      double omega = safeAcos(math::dot(nVec, eVec) / (nMag * e));
      // Use orientation to disambiguate [0, 2pi)
      if (math::dot(math::cross(nVec, eVec), h) < 0.0) omega = 2.0 * math::kPi - omega;
      out.argPeriapsisRad = wrapTwoPi(omega);
    } else {
      // Equatorial: Ω undefined; use longitude of periapsis instead.
      out.argPeriapsisRad = wrapTwoPi(std::atan2(eVec.y, eVec.x));
    }

    // True anomaly ν
    double nu = safeAcos(math::dot(eVec, relPosKm) / (e * r));
    if (math::dot(relPosKm, relVelKmS) < 0.0) nu = 2.0 * math::kPi - nu;
    out.trueAnomalyRad = wrapTwoPi(nu);
  } else {
    // Circular: ω undefined. Treat ν as argument of latitude u.
    out.argPeriapsisRad = 0.0;
    if (!out.equatorial) {
      double u = safeAcos(math::dot(nVec, relPosKm) / (nMag * r));
      if (math::dot(math::cross(nVec, relPosKm), h) < 0.0) u = 2.0 * math::kPi - u;
      out.trueAnomalyRad = wrapTwoPi(u);
    } else {
      // Circular equatorial: all angles degenerate; use true longitude.
      out.trueAnomalyRad = wrapTwoPi(std::atan2(relPosKm.y, relPosKm.x));
    }
  }

  // Mean motion / anomalies.
  if (out.type == TwoBodyOrbit::Type::Elliptic && out.semiMajorAxisKm > 0.0) {
    const double a = out.semiMajorAxisKm;
    const double n = std::sqrt(muKm3S2 / (a * a * a));
    out.meanMotionRadPerSec = n;
    out.periodSec = (n > 0.0) ? (2.0 * math::kPi / n) : 0.0;

    const double E = eccentricAnomalyFromTrueAnomaly(e, out.trueAnomalyRad);
    out.eccentricAnomalyRad = E;
    out.meanAnomalyRad = wrapTwoPi(E - e * std::sin(E));
  } else if (out.type == TwoBodyOrbit::Type::Hyperbolic && out.semiMajorAxisKm < 0.0) {
    const double a = out.semiMajorAxisKm; // negative
    const double n = std::sqrt(muKm3S2 / ((-a) * (-a) * (-a)));
    out.meanMotionRadPerSec = n;

    // Hyperbolic anomaly H via: tanh(H/2) = sqrt((e-1)/(e+1)) * tan(nu/2)
    const double tanHalfNu = std::tan(out.trueAnomalyRad * 0.5);
    const double fac = std::sqrt((e - 1.0) / (e + 1.0));
    double x = fac * tanHalfNu;
    x = math::clamp(x, -0.999999999999, 0.999999999999);
    const double H = 2.0 * std::atanh(x);
    out.eccentricAnomalyRad = H;
    out.meanAnomalyRad = e * std::sinh(H) - H;
  }

  out.valid = true;
  return out;
}

bool stateFromClassicalOrbitElementsAtTrueAnomaly(const ClassicalOrbitalElements& el,
                                                  double trueAnomalyRad,
                                                  math::Vec3d& outPosKm,
                                                  math::Vec3d& outVelKmS) {
  outPosKm = {0, 0, 0};
  outVelKmS = {0, 0, 0};

  if (!el.valid || !isFinite(el.muKm3S2) || el.muKm3S2 <= 0.0) return false;
  if (!isFinite(trueAnomalyRad)) return false;
  if (el.type == TwoBodyOrbit::Type::Parabolic) return false;

  const double e = el.eccentricity;
  const double p = el.parameterKm;
  if (!(p > 0.0) || !isFinite(p)) return false;

  const double nu = trueAnomalyRad;
  const double cosNu = std::cos(nu);
  const double sinNu = std::sin(nu);

  const double denom = 1.0 + e * cosNu;
  if (std::abs(denom) < 1e-18) return false;

  const double r = p / denom;

  // Perifocal frame.
  const math::Vec3d r_pf{r * cosNu, r * sinNu, 0.0};

  const double fac = std::sqrt(el.muKm3S2 / p);
  const math::Vec3d v_pf{-fac * sinNu, fac * (e + cosNu), 0.0};

  // Rotate from perifocal to inertial: Rz(Ω) * Rx(i) * Rz(ω).
  const math::Quatd qO = math::Quatd::fromAxisAngle({0, 0, 1}, el.raanRad);
  const math::Quatd qi = math::Quatd::fromAxisAngle({1, 0, 0}, el.inclinationRad);
  const math::Quatd qw = math::Quatd::fromAxisAngle({0, 0, 1}, el.argPeriapsisRad);
  const math::Quatd q = (qO * qi * qw).normalized();

  outPosKm = q * r_pf;
  outVelKmS = q * v_pf;

  return isFinite(outPosKm) && isFinite(outVelKmS);
}

bool stateFromClassicalOrbitElements(const ClassicalOrbitalElements& el,
                                     math::Vec3d& outPosKm,
                                     math::Vec3d& outVelKmS) {
  return stateFromClassicalOrbitElementsAtTrueAnomaly(el, el.trueAnomalyRad, outPosKm, outVelKmS);
}

double timeToTrueAnomalySecElliptic(const ClassicalOrbitalElements& el,
                                    double targetTrueAnomalyRad) {
  if (!el.valid) return -1.0;
  if (el.type != TwoBodyOrbit::Type::Elliptic) return -1.0;
  if (!(el.meanMotionRadPerSec > 0.0) || !isFinite(el.meanMotionRadPerSec)) return -1.0;
  const double e = el.eccentricity;
  if (!(e >= 0.0 && e < 1.0)) return -1.0;

  const double M0 = wrapTwoPi(el.meanAnomalyRad);
  const double M1 = meanAnomalyFromTrueAnomaly(e, wrapTwoPi(targetTrueAnomalyRad));

  double dM = M1 - M0;
  if (dM < 0.0) dM += 2.0 * math::kPi;

  return dM / el.meanMotionRadPerSec;
}

} // namespace stellar::sim
