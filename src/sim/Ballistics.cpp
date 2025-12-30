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

} // namespace stellar::sim
