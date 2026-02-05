#include "stellar/sim/CollisionWarning.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::sim {

double stoppingTimeSec(double speedKmS, double decelKmS2) {
  if (speedKmS <= 0.0) return 0.0;
  if (decelKmS2 <= 1e-12) return std::numeric_limits<double>::infinity();
  return speedKmS / decelKmS2;
}

double stoppingDistanceKm(double speedKmS, double decelKmS2) {
  if (speedKmS <= 0.0) return 0.0;
  if (decelKmS2 <= 1e-12) return std::numeric_limits<double>::infinity();
  // d = v^2 / (2a)
  return (speedKmS * speedKmS) / (2.0 * decelKmS2);
}

static double clamp01(double x) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  return x;
}

CollisionWarningResult computeCollisionWarning(const ProximityFieldKm& field,
                                              const math::Vec3d& posKm,
                                              const math::Vec3d& velKmS,
                                              double maxDecelKmS2,
                                              const CollisionWarningParams& params,
                                              int ignoreId) {
  CollisionWarningResult out{};

  out.speedKmS = velKmS.length();
  out.maxDecelKmS2 = maxDecelKmS2;

  if (out.speedKmS < std::max(0.0, params.minSpeedKmS)) {
    return out;
  }

  const double horizon = std::max(0.0, params.horizonSec);
  if (horizon <= 1e-9 || field.empty()) {
    return out;
  }

  const double padKm = std::max(0.0, params.padKm);

  const ProximityHitKm hit = field.predictLinearImpact(posKm, velKmS, horizon, padKm, ignoreId);
  if (!hit.hit) {
    return out;
  }

  out.hasImpact = true;
  out.obstacleId = hit.id;
  out.ttiSec = hit.tSec;
  out.impactDistKm = hit.tKm;
  out.impactPointKm = hit.pointKm;

  out.stopTimeSec = stoppingTimeSec(out.speedKmS, maxDecelKmS2);
  out.stopDistKm = stoppingDistanceKm(out.speedKmS, maxDecelKmS2);
  out.stopDistWithMarginKm = out.stopDistKm * (1.0 + std::max(0.0, params.stopMarginFactor));

  // Treat numerical weirdness as "no stop model".
  if (!std::isfinite(out.stopDistWithMarginKm)) {
    out.stopDistWithMarginKm = 0.0;
  }

  out.marginKm = out.impactDistKm - out.stopDistWithMarginKm;
  out.canStopBeforeImpact = (out.marginKm >= 0.0);

  // Compute hazard scalar based on time-to-impact and/or stopping distance.
  double hazardTime = 0.0;
  {
    const double danger = std::max(0.0, params.dangerTtiSec);
    const double caution = std::max(danger + 1e-6, params.cautionTtiSec);
    const double t = out.ttiSec;
    // 0 when t >= caution, 1 when t <= danger.
    hazardTime = clamp01((caution - t) / (caution - danger));
  }

  double hazardStop = 0.0;
  if (params.useStopDistance && out.stopDistWithMarginKm > 1e-9) {
    // 0 when impactDist >= stopDistWithMargin, 1 when impactDist <= 0.
    hazardStop = clamp01((out.stopDistWithMarginKm - out.impactDistKm) / out.stopDistWithMarginKm);
  }

  out.hazard01 = std::max(hazardTime, hazardStop);

  // Discrete level.
  const double danger = std::max(0.0, params.dangerTtiSec);
  const double caution = std::max(danger + 1e-6, params.cautionTtiSec);

  if (!out.hasImpact) {
    out.level = CollisionWarningLevel::None;
  } else if ((out.ttiSec <= danger + 1e-6) || (params.useStopDistance && !out.canStopBeforeImpact)) {
    out.level = CollisionWarningLevel::Danger;
  } else if (out.ttiSec <= caution) {
    out.level = CollisionWarningLevel::Caution;
  } else {
    out.level = CollisionWarningLevel::Advisory;
  }

  return out;
}

} // namespace stellar::sim
