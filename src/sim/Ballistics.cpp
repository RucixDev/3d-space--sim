#include "stellar/sim/Ballistics.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

std::optional<double> solveInterceptTimeSec(const math::Vec3d& shooterPosKm,
                                           const math::Vec3d& shooterVelKmS,
                                           const math::Vec3d& targetPosKm,
                                           const math::Vec3d& targetVelKmS,
                                           double projectileSpeedKmS,
                                           double minTimeSec) {
  const double s = std::max(1.0e-9, projectileSpeedKmS);
  const double tMin = std::max(0.0, minTimeSec);

  // Solve |r + v t| = s t (see header).
  const math::Vec3d r = targetPosKm - shooterPosKm;
  const math::Vec3d v = targetVelKmS - shooterVelKmS;

  const double a = v.lengthSq() - s * s;
  const double b = 2.0 * math::dot(r, v);
  const double c = r.lengthSq();

  // Degenerate: target is basically at the shooter.
  if (c < 1.0e-18) {
    return std::nullopt;
  }

  constexpr double kEps = 1.0e-12;

  // When a ~ 0 the quadratic collapses to a linear form.
  if (std::abs(a) < kEps) {
    if (std::abs(b) < kEps) {
      return std::nullopt;
    }

    const double t = -c / b;
    if (t >= tMin) return t;
    return std::nullopt;
  }

  const double disc = b * b - 4.0 * a * c;
  if (disc < 0.0) return std::nullopt;

  const double sd = std::sqrt(std::max(0.0, disc));
  const double inv2a = 1.0 / (2.0 * a);

  const double t1 = (-b - sd) * inv2a;
  const double t2 = (-b + sd) * inv2a;

  // Pick the smallest positive root.
  double best = 1.0e300;
  if (t1 >= tMin) best = std::min(best, t1);
  if (t2 >= tMin) best = std::min(best, t2);

  if (best < 1.0e200) return best;
  return std::nullopt;
}

std::optional<LeadSolution> solveProjectileLead(const math::Vec3d& shooterPosKm,
                                                const math::Vec3d& shooterVelKmS,
                                                const math::Vec3d& targetPosKm,
                                                const math::Vec3d& targetVelKmS,
                                                double projectileSpeedKmS,
                                                double maxTimeSec,
                                                double minTimeSec) {
  const double maxT = std::max(0.0, maxTimeSec);
  auto tOpt = solveInterceptTimeSec(shooterPosKm, shooterVelKmS, targetPosKm, targetVelKmS,
                                    projectileSpeedKmS, minTimeSec);
  if (!tOpt) return std::nullopt;
  const double t = *tOpt;
  if (t > maxT) return std::nullopt;

  LeadSolution out{};
  out.tSec = t;
  out.leadPointKm = targetPosKm + targetVelKmS * t;
  math::Vec3d d = out.leadPointKm - shooterPosKm;
  if (d.lengthSq() < 1.0e-18) d = {0, 0, 1};
  out.aimDirWorld = d.normalized();
  return out;
}



std::optional<double> solveInterceptTimeSecAccel(const math::Vec3d& shooterPosKm,
                                                const math::Vec3d& shooterVelKmS,
                                                const math::Vec3d& targetPosKm,
                                                const math::Vec3d& targetVelKmS,
                                                const math::Vec3d& targetAccelKmS2,
                                                double projectileSpeedKmS,
                                                double maxTimeSec,
                                                double minTimeSec) {
  const double s = std::max(1.0e-9, projectileSpeedKmS);
  const double tMin = std::max(0.0, minTimeSec);
  const double tMax = std::max(0.0, maxTimeSec);
  if (!(tMax >= tMin)) return std::nullopt;

  const math::Vec3d r0 = targetPosKm - shooterPosKm;
  if (r0.lengthSq() < 1.0e-18) {
    return std::nullopt;
  }

  const math::Vec3d v = targetVelKmS - shooterVelKmS;
  const math::Vec3d a = targetAccelKmS2;

  // If acceleration is negligible, fall back to the closed-form quadratic.
  if (a.lengthSq() < 1.0e-24) {
    auto t = solveInterceptTimeSec(shooterPosKm, shooterVelKmS, targetPosKm, targetVelKmS, s, tMin);
    if (t && *t <= tMax) return t;
    return std::nullopt;
  }

  auto poly = [&](double t) -> double {
    const double tt = t * t;
    const math::Vec3d rt = r0 + v * t + a * (0.5 * tt);
    const double lhs = rt.lengthSq();
    const double rhs = (s * t) * (s * t);
    return lhs - rhs;
  };

  // Search for the earliest root by scanning for the first entry into the "catchable"
  // region (poly <= 0). This avoids needing an explicit quartic solver.
  const int steps = 160;
  const double span = tMax - tMin;
  double prevT = tMin;
  double prevG = poly(prevT);

  constexpr double kZeroTol = 1.0e-10;

  if (std::abs(prevG) <= kZeroTol) return prevT;

  for (int i = 1; i <= steps; ++i) {
    const double t = tMin + span * (double)i / (double)steps;
    const double g = poly(t);

    if (std::abs(g) <= kZeroTol) return t;

    if (prevG > 0.0 && g <= 0.0) {
      // Bracket [prevT, t] contains the first root.
      double lo = prevT;
      double hi = t;
      double gLo = prevG;
      double gHi = g;

      // Bisection for robustness.
      for (int it = 0; it < 64; ++it) {
        const double mid = 0.5 * (lo + hi);
        const double gMid = poly(mid);
        if (gMid > 0.0) {
          lo = mid;
          gLo = gMid;
        } else {
          hi = mid;
          gHi = gMid;
        }
      }

      (void)gLo;
      (void)gHi;
      return hi;
    }

    prevT = t;
    prevG = g;
  }

  return std::nullopt;
}

std::optional<LeadSolution> solveProjectileLeadAccel(const math::Vec3d& shooterPosKm,
                                                     const math::Vec3d& shooterVelKmS,
                                                     const math::Vec3d& targetPosKm,
                                                     const math::Vec3d& targetVelKmS,
                                                     const math::Vec3d& targetAccelKmS2,
                                                     double projectileSpeedKmS,
                                                     double maxTimeSec,
                                                     double minTimeSec) {
  const double maxT = std::max(0.0, maxTimeSec);

  auto tOpt = solveInterceptTimeSecAccel(shooterPosKm, shooterVelKmS, targetPosKm, targetVelKmS,
                                        targetAccelKmS2, projectileSpeedKmS, maxT, minTimeSec);
  if (!tOpt) return std::nullopt;

  const double t = *tOpt;
  if (t > maxT) return std::nullopt;

  LeadSolution out{};
  out.tSec = t;
  out.leadPointKm = targetPosKm + targetVelKmS * t + targetAccelKmS2 * (0.5 * t * t);
  math::Vec3d d = out.leadPointKm - shooterPosKm;
  if (d.lengthSq() < 1.0e-18) d = {0, 0, 1};
  out.aimDirWorld = d.normalized();
  return out;
}

} // namespace stellar::sim
