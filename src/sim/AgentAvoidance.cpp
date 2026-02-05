#include "stellar/sim/AgentAvoidance.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::sim {

namespace {

static bool isFiniteVec(const math::Vec3d& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

static double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

static double safeAcos(double x) {
  return std::acos(std::clamp(x, -1.0, 1.0));
}

static math::Vec3d slerpUnit(const math::Vec3d& aUnit, const math::Vec3d& bUnit, double t) {
  t = std::clamp(t, 0.0, 1.0);
  const double d = std::clamp(math::dot(aUnit, bUnit), -1.0, 1.0);

  // If the directions are very close, fall back to lerp to avoid precision issues.
  const double theta = safeAcos(d);
  if (theta < 1e-6) {
    const auto v = (aUnit * (1.0 - t)) + (bUnit * t);
    const auto n = v.normalized();
    return (n.lengthSq() > 1e-12) ? n : aUnit;
  }

  const double sinTheta = std::sin(theta);
  const double wA = std::sin((1.0 - t) * theta) / sinTheta;
  const double wB = std::sin(t * theta) / sinTheta;
  const auto v = (aUnit * wA) + (bUnit * wB);
  const auto n = v.normalized();
  return (n.lengthSq() > 1e-12) ? n : aUnit;
}

static math::Vec3d stablePerp(const math::Vec3d& dirUnit) {
  // Pick a world axis that isn't nearly parallel with dirUnit, then cross.
  const math::Vec3d axis = (std::abs(dirUnit.y) < 0.90) ? math::Vec3d{0, 1, 0} : math::Vec3d{1, 0, 0};
  auto v = math::cross(axis, dirUnit);
  const auto n = v.normalized();
  if (n.lengthSq() > 1e-12) return n;
  v = math::cross(math::Vec3d{0, 0, 1}, dirUnit);
  return v.normalized();
}

} // namespace

AgentAvoidanceResult steerAvoidAgents(const math::Vec3d& selfPosKm,
                                     const math::Vec3d& selfVelKmS,
                                     const math::Vec3d& desiredDirUnitWorld,
                                     double desiredSpeedKmS,
                                     std::span<const AgentSphere> neighbors,
                                     core::u64 seed,
                                     const AgentAvoidanceParams& params) {
  AgentAvoidanceResult out{};

  math::Vec3d desiredN = desiredDirUnitWorld.normalized();
  if (desiredN.lengthSq() <= 1e-12) {
    desiredN = {0, 0, 1};
  }

  out.desiredDirUnit = desiredN;
  out.safeDirUnit = desiredN;

  if (!params.enabled || neighbors.empty()) {
    return out;
  }

  const double horizon = std::max(0.0, params.horizonSec);
  if (horizon <= 1e-6) {
    return out;
  }

  // Predict our near-future velocity as a blend between current motion and the
  // commanded direction. This keeps the model stable when accelerating from rest.
  const double blend = clamp01(params.selfVelBlend01);
  const double spDes = std::max(0.0, desiredSpeedKmS);
  const math::Vec3d vCmd = desiredN * spDes;

  math::Vec3d vSelfPred = selfVelKmS * (1.0 - blend) + vCmd * blend;
  if (!isFiniteVec(vSelfPred)) vSelfPred = selfVelKmS;
  if (vSelfPred.lengthSq() <= 1e-12 && spDes > 1e-9) vSelfPred = vCmd;

  const double inflateBase = std::max(0.0, params.selfRadiusKm) + std::max(0.0, params.padKm);
  const double nearMiss = std::max(0.0, params.nearMissExtraKm);
  const double minRelSp = std::max(0.0, params.minRelSpeedKmS);

  math::Vec3d repulse{0, 0, 0};

  core::u64 bestThreatId = 0;
  double bestThreatClear = std::numeric_limits<double>::infinity();
  double bestThreatT = 0.0;
  double bestThreatDist = 0.0;

  for (const auto& n : neighbors) {
    if (n.id == 0) {
      // id=0 is treated as "unknown"; still valid but we can't report it reliably.
    }

    if (!isFiniteVec(n.posKm) || !isFiniteVec(n.velKmS) || !std::isfinite(n.radiusKm) || !std::isfinite(n.hardness01)) {
      continue;
    }

    const double r = std::max(0.0, n.radiusKm);
    if (r <= 0.0) {
      continue;
    }

    const double hard = clamp01(n.hardness01);
    if (hard <= 0.0) {
      continue;
    }

    // Relative motion of the neighbor in our predicted frame.
    const math::Vec3d dp = n.posKm - selfPosKm;
    const math::Vec3d dv = n.velKmS - vSelfPred;

    double tClosest = 0.0;
    const double relSp = dv.length();
    if (relSp >= minRelSp) {
      const double dv2 = math::dot(dv, dv);
      if (dv2 > 1e-12) {
        tClosest = -math::dot(dp, dv) / dv2;
        tClosest = std::clamp(tClosest, 0.0, horizon);
      }
    }

    const math::Vec3d sep = dp + dv * tClosest;
    const double dist = sep.length();

    // Clearance between inflated spheres at closest approach.
    const double rInfl = r + inflateBase;
    const double clear = dist - rInfl;

    if (clear < bestThreatClear) {
      bestThreatClear = clear;
      bestThreatId = n.id;
      bestThreatT = tClosest;
      bestThreatDist = dist;
    }

    if (clear > nearMiss) {
      continue;
    }

    double w = 1.0;
    if (nearMiss > 1e-9) {
      w = clamp01((nearMiss - clear) / nearMiss);
    } else {
      w = (clear <= 0.0) ? 1.0 : 0.0;
    }

    // Prefer resolving sooner threats.
    const double timeW = (horizon > 1e-9) ? clamp01(1.0 - (tClosest / horizon)) : 1.0;

    w = w * w * timeW * hard;
    if (w <= 0.0) {
      continue;
    }

    math::Vec3d away = -sep;
    if (away.lengthSq() <= 1e-12) {
      away = stablePerp(desiredN);
      // Deterministic left/right choice in the perfectly symmetric case.
      const core::u64 h = core::hashCombine(seed, n.id);
      if (h & 1ull) away = away * -1.0;
    }

    const auto awayN = away.normalized();
    if (awayN.lengthSq() > 1e-12) {
      repulse = repulse + awayN * w;
    }
  }

  if (bestThreatId != 0) {
    out.threatId = bestThreatId;
    out.threatTtiSec = bestThreatT;
    out.threatClearanceKm = bestThreatClear;
    out.threatDistAtClosestKm = bestThreatDist;
  }

  const double repLen = repulse.length();
  if (repLen <= 1e-9) {
    return out;
  }

  const math::Vec3d steerDir = repulse / repLen;
  const double steerStrength = std::max(0.0, params.strength) * clamp01(repLen);

  math::Vec3d cand = (desiredN + steerDir * steerStrength).normalized();
  if (cand.lengthSq() <= 1e-12) {
    cand = desiredN;
  }

  const double maxDeg = std::max(0.0, params.maxSteerDeg);
  if (maxDeg > 1e-6) {
    const double maxRad = maxDeg * (3.14159265358979323846 / 180.0);
    const double ang = safeAcos(math::dot(desiredN, cand));
    if (ang > maxRad && ang > 1e-9) {
      const double t = std::clamp(maxRad / ang, 0.0, 1.0);
      cand = slerpUnit(desiredN, cand, t);
    }
  }

  out.safeDirUnit = cand;
  out.steering = (math::dot(desiredN, cand) < 0.99999);
  return out;
}

} // namespace stellar::sim
