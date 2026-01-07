#include "stellar/sim/OrbitalMechanics.h"
#include "stellar/sim/Gravity.h"

#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps = 1e-6) { return std::abs(a - b) <= eps; }

static bool near3(const stellar::math::Vec3d& a,
                  const stellar::math::Vec3d& b,
                  double eps = 1e-6) {
  return near(a.x, b.x, eps) && near(a.y, b.y, eps) && near(a.z, b.z, eps);
}

static double wrapTwoPi(double a) {
  const double twoPi = 2.0 * stellar::math::kPi;
  a = std::fmod(a, twoPi);
  if (a < 0.0) a += twoPi;
  return a;
}

static double angleDiff(double a, double b) {
  // Smallest signed difference in [-pi, pi]
  const double twoPi = 2.0 * stellar::math::kPi;
  double d = wrapTwoPi(a) - wrapTwoPi(b);
  if (d > stellar::math::kPi) d -= twoPi;
  if (d < -stellar::math::kPi) d += twoPi;
  return d;
}

static stellar::math::Vec3d rotZ(const stellar::math::Vec3d& v, double a) {
  const double c = std::cos(a);
  const double s = std::sin(a);
  return {c * v.x - s * v.y, s * v.x + c * v.y, v.z};
}

static stellar::math::Vec3d rotX(const stellar::math::Vec3d& v, double a) {
  const double c = std::cos(a);
  const double s = std::sin(a);
  return {v.x, c * v.y - s * v.z, s * v.y + c * v.z};
}

static stellar::math::Vec3d perifocalToInertial(const stellar::math::Vec3d& v,
                                                double raan,
                                                double inc,
                                                double argPeri) {
  // Rz(Ω) Rx(i) Rz(ω) v
  return rotZ(rotX(rotZ(v, argPeri), inc), raan);
}

static void stateFromElementsManual(const stellar::sim::ClassicalOrbitalElements& el,
                                    double nu,
                                    stellar::math::Vec3d& outR,
                                    stellar::math::Vec3d& outV) {
  const double e = el.eccentricity;
  const double p = el.parameterKm;

  const double cosNu = std::cos(nu);
  const double sinNu = std::sin(nu);

  const double r = p / (1.0 + e * cosNu);

  // Perifocal
  const stellar::math::Vec3d r_pf{r * cosNu, r * sinNu, 0.0};

  const double fac = std::sqrt(el.muKm3S2 / p);
  const stellar::math::Vec3d v_pf{-fac * sinNu, fac * (e + cosNu), 0.0};

  outR = perifocalToInertial(r_pf, el.raanRad, el.inclinationRad, el.argPeriapsisRad);
  outV = perifocalToInertial(v_pf, el.raanRad, el.inclinationRad, el.argPeriapsisRad);
}

static stellar::sim::ClassicalOrbitalElements makeElliptic(double mu,
                                                           double a,
                                                           double e,
                                                           double inc,
                                                           double raan,
                                                           double argPeri,
                                                           double nu) {
  stellar::sim::ClassicalOrbitalElements el{};
  el.valid = true;
  el.muKm3S2 = mu;
  el.type = stellar::sim::TwoBodyOrbit::Type::Elliptic;
  el.semiMajorAxisKm = a;
  el.eccentricity = e;
  el.parameterKm = a * (1.0 - e * e);
  el.inclinationRad = inc;
  el.raanRad = raan;
  el.argPeriapsisRad = argPeri;
  el.trueAnomalyRad = nu;

  el.meanMotionRadPerSec = std::sqrt(mu / (a * a * a));
  el.periodSec = (el.meanMotionRadPerSec > 0.0) ? (2.0 * stellar::math::kPi / el.meanMotionRadPerSec) : 0.0;

  // E and M (elliptic)
  const double t = std::tan(nu * 0.5);
  const double fac = std::sqrt((1.0 - e) / (1.0 + e));
  const double x = fac * t;
  const double E = wrapTwoPi(2.0 * std::atan(x));
  el.eccentricAnomalyRad = E;
  el.meanAnomalyRad = wrapTwoPi(E - e * std::sin(E));

  return el;
}

int test_orbital_elements() {
  int fails = 0;
  const double mu = stellar::sim::muFromEarthMass(1.0);

  // --- General inclined ellipse round-trip ---
  {
    const double a = 12000.0;
    const double e = 0.2;
    const double inc = stellar::math::degToRad(30.0);
    const double raan = stellar::math::degToRad(40.0);
    const double argPeri = stellar::math::degToRad(60.0);
    const double nu = stellar::math::degToRad(10.0);

    auto el = makeElliptic(mu, a, e, inc, raan, argPeri, nu);

    stellar::math::Vec3d rM{}, vM{};
    stateFromElementsManual(el, nu, rM, vM);

    // stateFromClassicalOrbitElements should match the manual matrix rotation.
    stellar::math::Vec3d rS{}, vS{};
    if (!stellar::sim::stateFromClassicalOrbitElementsAtTrueAnomaly(el, nu, rS, vS)) {
      std::cerr << "[test_orbital_elements] stateFromClassicalOrbitElements failed (general ellipse)\n";
      ++fails;
    } else {
      const double errR = (rS - rM).length();
      const double errV = (vS - vM).length();
      if (errR > 1e-6 || errV > 1e-9) {
        std::cerr << "[test_orbital_elements] stateFrom mismatch (general ellipse) errR=" << errR
                  << " km errV=" << errV << " km/s\n";
        ++fails;
      }
    }

    // Solve elements from the generated state.
    const auto solved = stellar::sim::solveClassicalOrbitElements(rM, vM, mu);
    if (!solved.valid || solved.type != stellar::sim::TwoBodyOrbit::Type::Elliptic) {
      std::cerr << "[test_orbital_elements] solveClassicalOrbitElements failed to classify ellipse\n";
      ++fails;
    } else {
      if (!near(solved.semiMajorAxisKm, a, 1e-3)) {
        std::cerr << "[test_orbital_elements] a mismatch got " << solved.semiMajorAxisKm << " expected " << a << "\n";
        ++fails;
      }
      if (!near(solved.eccentricity, e, 3e-6)) {
        std::cerr << "[test_orbital_elements] e mismatch got " << solved.eccentricity << " expected " << e << "\n";
        ++fails;
      }
      if (std::abs(angleDiff(solved.inclinationRad, inc)) > 2e-6) {
        std::cerr << "[test_orbital_elements] inc mismatch got " << solved.inclinationRad << " expected " << inc << "\n";
        ++fails;
      }
      if (std::abs(angleDiff(solved.raanRad, raan)) > 2e-6) {
        std::cerr << "[test_orbital_elements] raan mismatch got " << solved.raanRad << " expected " << raan << "\n";
        ++fails;
      }
      if (std::abs(angleDiff(solved.argPeriapsisRad, argPeri)) > 2e-5) {
        std::cerr << "[test_orbital_elements] argPeri mismatch got " << solved.argPeriapsisRad << " expected " << argPeri << "\n";
        ++fails;
      }
      if (std::abs(angleDiff(solved.trueAnomalyRad, nu)) > 2e-5) {
        std::cerr << "[test_orbital_elements] nu mismatch got " << solved.trueAnomalyRad << " expected " << nu << "\n";
        ++fails;
      }
    }

    // Node sanity: ascending node should be at ω+ν = 0
    const double nuAsc = wrapTwoPi(-argPeri);
    const double dtAsc = stellar::sim::timeToTrueAnomalySecElliptic(el, nuAsc);
    if (!(dtAsc >= 0.0)) {
      std::cerr << "[test_orbital_elements] timeToTrueAnomalySecElliptic failed for ascending node\n";
      ++fails;
    } else {
      stellar::math::Vec3d rAsc{}, vAsc{};
      stateFromElementsManual(el, nuAsc, rAsc, vAsc);
      if (std::abs(rAsc.z) > 1e-5) {
        std::cerr << "[test_orbital_elements] ascending node z not ~0, z=" << rAsc.z << "\n";
        ++fails;
      }
      if (!(vAsc.z > 0.0)) {
        std::cerr << "[test_orbital_elements] ascending node vz should be >0, vz=" << vAsc.z << "\n";
        ++fails;
      }
    }
  }

  // --- Circular inclined orbit should be detected as circular ---
  {
    const double r = 8000.0;
    const double v = std::sqrt(mu / r);
    const stellar::math::Vec3d relPos{r, 0, 0};
    const stellar::math::Vec3d relVel{0, v * std::cos(stellar::math::degToRad(45.0)),
                                      v * std::sin(stellar::math::degToRad(45.0))};

    const auto solved = stellar::sim::solveClassicalOrbitElements(relPos, relVel, mu);
    if (!solved.valid) {
      std::cerr << "[test_orbital_elements] circular inclined solve invalid\n";
      ++fails;
    } else if (!solved.circular) {
      std::cerr << "[test_orbital_elements] expected circular flag for e~0\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_orbital_elements] pass\n";
  return fails;
}
