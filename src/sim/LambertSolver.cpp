#include "stellar/sim/LambertSolver.h"

#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::sim {

// Stumpff functions C(z) and S(z) for universal variables.
// Uses series expansions for |z| ~ 0 to avoid catastrophic cancellation.
static void stumpff(double z, double& C, double& S) {
  const double az = std::abs(z);
  if (az < 1e-6) {
    // Series:
    //   C(z) = 1/2 - z/24 + z^2/720 - ...
    //   S(z) = 1/6 - z/120 + z^2/5040 - ...
    const double z2 = z * z;
    C = 0.5 - z / 24.0 + z2 / 720.0;
    S = (1.0 / 6.0) - z / 120.0 + z2 / 5040.0;
    return;
  }

  if (z > 0.0) {
    const double s = std::sqrt(z);
    C = (1.0 - std::cos(s)) / z;
    S = (s - std::sin(s)) / (s * s * s);
  } else {
    const double s = std::sqrt(-z);
    C = (std::cosh(s) - 1.0) / (-z);
    S = (std::sinh(s) - s) / (s * s * s);
  }
}

// Evaluate time-of-flight for a given z.
// Returns false when y(z) is negative (no real solution for this z).
static bool tofForZ(double z,
                    double r1,
                    double r2,
                    double A,
                    double sqrtMu,
                    double& outT,
                    double& outY) {
  double C = 0.0, S = 0.0;
  stumpff(z, C, S);
  if (!(C > 0.0)) return false;

  const double sqrtC = std::sqrt(C);

  // Vallado-style formulation.
  const double y = r1 + r2 + (A * (z * S - 1.0)) / std::max(1e-12, sqrtC);
  outY = y;
  if (!(y > 0.0)) return false;

  const double x = std::sqrt(y / C);
  const double t = (x * x * x * S + A * std::sqrt(y)) / sqrtMu;

  if (!std::isfinite(t)) return false;
  outT = t;
  return true;
}

LambertResult solveLambertUniversal(const math::Vec3d& r1Km,
                                   const math::Vec3d& r2Km,
                                   double dtSec,
                                   double muKm3S2,
                                   const LambertOptions& opt) {
  LambertResult res{};

  if (!(dtSec > 0.0) || !(muKm3S2 > 0.0)) return res;

  const double r1 = r1Km.length();
  const double r2 = r2Km.length();
  if (!(r1 > 1e-9) || !(r2 > 1e-9)) return res;

  // Transfer angle selection.
  double cosD = math::dot(r1Km, r2Km) / (r1 * r2);
  cosD = std::clamp(cosD, -1.0, 1.0);

  const double sinDMag = std::sqrt(std::max(0.0, 1.0 - cosD * cosD));
  if (sinDMag < 1e-10) {
    // Colinear (0 or 180 degrees) -> degenerate for this simple solver.
    return res;
  }

  // Direction (sign of sinD) from reference normal.
  math::Vec3d refN = opt.refNormal;
  if (refN.lengthSq() < 1e-18) {
    refN = math::cross(r1Km, r2Km);
  }

  double sinSign = 1.0;
  {
    const math::Vec3d c = math::cross(r1Km, r2Km);
    const double s = math::dot(c, refN);
    sinSign = (s >= 0.0) ? 1.0 : -1.0;
    if (!opt.prograde) sinSign = -sinSign;
    if (opt.longWay) sinSign = -sinSign;
  }

  const double sinD = sinDMag * sinSign;

  const double oneMinusCos = std::max(1e-12, 1.0 - cosD);
  const double A = sinD * std::sqrt((r1 * r2) / oneMinusCos);
  if (std::abs(A) < 1e-12) return res;

  const double sqrtMu = std::sqrt(muKm3S2);

  // --- Bracket + bisection on z ---
  // Start at z=0 (parabolic). Depending on dt, the root may be on either side.
  double t0 = 0.0, y0 = 0.0;
  const bool ok0 = tofForZ(0.0, r1, r2, A, sqrtMu, t0, y0);

  // If z=0 is invalid (y<0), we still try to find a valid bracket by scanning.
  auto evalF = [&](double z, double& F) -> bool {
    double t = 0.0, y = 0.0;
    if (!tofForZ(z, r1, r2, A, sqrtMu, t, y)) return false;
    F = t - dtSec;
    return true;
  };

  // Pick an initial bracket [zLo, zHi] such that F(zLo) <= 0 <= F(zHi).
  double zLo = 0.0, zHi = 0.0;
  double fLo = 0.0, fHi = 0.0;

  if (ok0) {
    const double f0 = t0 - dtSec;
    if (std::abs(f0) <= opt.tolSec) {
      // Direct parabolic solution.
      zLo = zHi = 0.0;
      fLo = fHi = f0;
    } else if (f0 < 0.0) {
      // Need a larger time-of-flight -> z > 0.
      zLo = 0.0;
      fLo = f0;
      zHi = 4.0;
      for (int i = 0; i < 64; ++i) {
        if (evalF(zHi, fHi) && fHi >= 0.0) break;
        zHi *= 2.0;
        if (zHi > 1e7) return res;
      }
      if (!(fHi >= 0.0)) return res;
    } else {
      // Need a smaller time-of-flight -> z < 0.
      zHi = 0.0;
      fHi = f0;
      zLo = -4.0;
      for (int i = 0; i < 64; ++i) {
        if (evalF(zLo, fLo) && fLo <= 0.0) break;
        zLo *= 2.0;
        if (zLo < -1e7) return res;
      }
      if (!(fLo <= 0.0)) return res;
    }
  } else {
    // Scan outward to find a bracket.
    // We start with small |z| and expand geometrically.
    bool haveLo = false, haveHi = false;
    for (int i = 0; i < 64 && !(haveLo && haveHi); ++i) {
      const double step = std::pow(2.0, i) * 0.5; // 0.5, 1, 2, 4, ...

      double f = 0.0;
      if (!haveHi && evalF(step, f)) {
        if (f >= 0.0) { zHi = step; fHi = f; haveHi = true; }
      }
      if (!haveLo && evalF(-step, f)) {
        if (f <= 0.0) { zLo = -step; fLo = f; haveLo = true; }
      }
    }
    if (!(haveLo && haveHi)) return res;
  }

  // If we landed exactly at z=0, we still need velocities; run the normal solve below.

  // Bisection.
  double zMid = 0.0;
  double fMid = 0.0;
  bool midOk = false;
  int iters = 0;

  for (iters = 0; iters < std::max(8, opt.maxIterations); ++iters) {
    zMid = 0.5 * (zLo + zHi);
    midOk = evalF(zMid, fMid);
    if (!midOk) {
      // If mid is invalid (y<0), move upward (toward larger z) since y(z)
      // tends to increase with z for this formulation.
      zLo = zMid;
      continue;
    }

    if (std::abs(fMid) <= opt.tolSec) break;

    if (fMid < 0.0) {
      zLo = zMid;
      fLo = fMid;
    } else {
      zHi = zMid;
      fHi = fMid;
    }
  }

  // Evaluate final z to get y and compute f/g.
  double tFin = 0.0, yFin = 0.0;
  if (!tofForZ(zMid, r1, r2, A, sqrtMu, tFin, yFin)) return res;

  const double f = 1.0 - yFin / r1;
  const double g = A * std::sqrt(yFin / muKm3S2);
  const double gdot = 1.0 - yFin / r2;

  if (std::abs(g) < 1e-12) return res;

  res.v1KmS = (r2Km - r1Km * f) / g;
  res.v2KmS = (r2Km * gdot - r1Km) / g;

  // Transfer angle diagnostic.
  const double a = std::acos(cosD);
  res.transferAngleRad = (sinD >= 0.0) ? a : (2.0 * math::kPi - a);

  res.ok = true;
  res.z = zMid;
  res.iterations = iters;
  return res;
}

} // namespace stellar::sim
