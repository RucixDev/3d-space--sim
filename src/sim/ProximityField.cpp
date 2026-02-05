#include "stellar/sim/ProximityField.h"

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

void ProximityFieldKm::clear() {
  obstacles_.clear();
  bvh_ = math::Bvh3d{};
}

void ProximityFieldKm::build(std::vector<SphereObstacleKm> obstacles, std::size_t leafSize) {
  obstacles_ = std::move(obstacles);

  std::vector<math::BvhItem3d> items;
  items.reserve(obstacles_.size());

  for (std::size_t i = 0; i < obstacles_.size(); ++i) {
    auto& o = obstacles_[i];

    if (!isFiniteVec(o.centerKm) || !std::isfinite(o.radiusKm) || !std::isfinite(o.hardness01)) {
      o.radiusKm = 0.0;
      o.hardness01 = 0.0;
      continue;
    }

    o.radiusKm = std::max(0.0, o.radiusKm);
    o.hardness01 = clamp01(o.hardness01);

    if (o.radiusKm <= 0.0) {
      continue;
    }

    const math::Vec3d r{ o.radiusKm, o.radiusKm, o.radiusKm };
    math::Aabb3d aabb;
    aabb.min = o.centerKm - r;
    aabb.max = o.centerKm + r;

    items.push_back(math::BvhItem3d{aabb, (int)i});
  }

  bvh_.build(items, leafSize);
}

const SphereObstacleKm* ProximityFieldKm::obstacle(int id) const {
  if (id < 0) return nullptr;
  const std::size_t idx = (std::size_t)id;
  if (idx >= obstacles_.size()) return nullptr;
  return &obstacles_[idx];
}

bool ProximityFieldKm::raySphereFirstHitKm_(const math::Vec3d& originKm,
                                           const math::Vec3d& dirUnitWorld,
                                           const math::Vec3d& centerKm,
                                           double radiusKm,
                                           double& outTKm) {
  // Algorithm from "Real-Time Collision Detection" (Christer Ericson):
  // - robust for rays starting outside the sphere and for inside hits.
  const math::Vec3d m = originKm - centerKm;
  const double b = math::dot(m, dirUnitWorld);
  const double c = math::dot(m, m) - radiusKm * radiusKm;

  // Ray origin outside sphere and pointing away.
  if (c > 0.0 && b > 0.0) {
    return false;
  }

  const double discr = b * b - c;
  if (discr < 0.0) {
    return false;
  }

  double t = -b - std::sqrt(discr);

  // If inside the sphere, clamp to 0.
  if (t < 0.0) t = 0.0;

  outTKm = t;
  return true;
}

ProximityHitKm ProximityFieldKm::raycastFirstHit(const math::Vec3d& originKm,
                                                const math::Vec3d& dirUnitWorld,
                                                double maxDistKm,
                                                double padKm,
                                                int ignoreId,
                                                int maxSphereTests) const {
  ProximityHitKm hit{};

  if (bvh_.empty()) {
    return hit;
  }

  maxDistKm = std::max(0.0, maxDistKm);
  if (maxDistKm <= 1e-12) {
    return hit;
  }

  const auto dirN = dirUnitWorld.normalized();
  if (dirN.lengthSq() <= 1e-12) {
    return hit;
  }

  padKm = std::max(0.0, padKm);
  if (maxSphereTests < 0) maxSphereTests = 0;

  const math::Vec3d endKm = originKm + dirN * maxDistKm;

  double bestT = maxDistKm + 1.0;
  int bestId = -1;
  double bestR = 0.0;
  double bestHard = 1.0;

  int tests = 0;
  bvh_.querySegment(originKm, endKm, [&](int id) {
    if (id == ignoreId) return true;
    const auto* o = obstacle(id);
    if (!o) return true;
    if (o->radiusKm <= 0.0) return true;

    if (maxSphereTests > 0 && tests >= maxSphereTests) {
      return true;
    }
    ++tests;

    const double r = o->radiusKm + padKm;
    double t;
    if (!raySphereFirstHitKm_(originKm, dirN, o->centerKm, r, t)) {
      return true;
    }

    if (t <= bestT && t <= maxDistKm + 1e-12) {
      bestT = t;
      bestId = id;
      bestR = o->radiusKm;
      bestHard = o->hardness01;
    }

    return true;
  });

  if (bestId >= 0) {
    hit.hit = true;
    hit.id = bestId;
    hit.tKm = bestT;
    hit.tSec = 0.0;
    hit.pointKm = originKm + dirN * bestT;
    hit.obstacleRadiusKm = bestR;
    hit.obstacleHardness01 = bestHard;
  }

  return hit;
}

ProximityHitKm ProximityFieldKm::predictLinearImpact(const math::Vec3d& originKm,
                                                    const math::Vec3d& velKmS,
                                                    double maxTimeSec,
                                                    double padKm,
                                                    int ignoreId,
                                                    int maxSphereTests) const {
  ProximityHitKm hit{};

  maxTimeSec = std::max(0.0, maxTimeSec);
  if (maxTimeSec <= 1e-12) {
    return hit;
  }

  const double sp = velKmS.length();
  if (sp <= 1e-9) {
    return hit;
  }

  const math::Vec3d dirN = velKmS / sp;
  const double maxDistKm = sp * maxTimeSec;

  hit = raycastFirstHit(originKm, dirN, maxDistKm, padKm, ignoreId, maxSphereTests);
  if (hit.hit) {
    hit.tSec = (sp > 1e-12) ? (hit.tKm / sp) : 0.0;
  }

  return hit;
}



TangentBypassResultKm planTangentBypassWaypoint(const ProximityFieldKm& field,
                                                const math::Vec3d& posKm,
                                                const math::Vec3d& goalKm,
                                                const TangentBypassParamsKm& params,
                                                int ignoreId,
                                                int preferredSideSign,
                                                int maxSphereTests) {
  TangentBypassResultKm out{};

  if (!params.enabled || field.empty()) {
    return out;
  }

  const math::Vec3d toGoal = goalKm - posKm;
  if (!isFiniteVec(toGoal)) {
    return out;
  }

  const double distGoal = toGoal.length();
  if (distGoal <= 1e-9) {
    return out;
  }

  const math::Vec3d dir = toGoal / distGoal;

  const double shipR = std::max(0.0, params.shipRadiusKm);
  const double pad = std::max(0.0, params.padKm);
  const double nearMiss = std::max(0.0, params.nearMissExtraKm);
  const double queryPad = shipR + pad + nearMiss;

  // If the direct segment is clear, no bypass is needed.
  const auto directHit = field.raycastFirstHit(posKm, dir, distGoal, queryPad, ignoreId, maxSphereTests);
  if (!directHit.hit) {
    return out;
  }

  const auto* o = field.obstacle(directHit.id);
  if (!o) {
    return out;
  }
  if (!isFiniteVec(o->centerKm) || !std::isfinite(o->radiusKm) || o->radiusKm <= 0.0) {
    return out;
  }

  // Inflate obstacle by the safety bubble used for the direct-hit query.
  const double rInfl = o->radiusKm + queryPad;

  // Compute ray enter/exit against the inflated sphere.
  const auto rayEnterExit = [&](double& outEnterKm, double& outExitKm) -> bool {
    const math::Vec3d m = posKm - o->centerKm;
    const double b = math::dot(m, dir);
    const double c = math::dot(m, m) - rInfl * rInfl;
    const double discr = b * b - c;
    if (discr < 0.0) return false;

    const double s = std::sqrt(std::max(0.0, discr));
    outEnterKm = -b - s;
    outExitKm = -b + s;
    return std::isfinite(outEnterKm) && std::isfinite(outExitKm);
  };

  double tEnterKm = 0.0;
  double tExitKm = 0.0;
  if (!rayEnterExit(tEnterKm, tExitKm)) {
    return out;
  }

  if (tExitKm < 0.0) {
    return out;
  }

  // Determine a stable lateral direction (in the plane normal to dir) to step away from the obstacle.
  // Use the direction from obstacle center -> closest point on the ray if possible.
  const math::Vec3d toC = o->centerKm - posKm;
  const double tca = math::dot(toC, dir);
  const math::Vec3d closest = posKm + dir * tca;
  math::Vec3d d = closest - o->centerKm;

  math::Vec3d sideRef{0, 0, 0};
  if (d.lengthSq() > 1e-12) {
    sideRef = d.normalized();
  } else {
    // Degenerate when the ray passes through the center: pick a deterministic perpendicular.
    sideRef = stablePerp(dir);
  }

  struct Candidate {
    bool valid{false};
    math::Vec3d waypointKm{0, 0, 0};
    double score{0.0};
    int sideSign{0};
  };

  auto evalCandidate = [&](int sideSign) -> Candidate {
    Candidate c{};
    c.sideSign = sideSign;

    const double clearanceKm = std::max(1e-9, rInfl);
    const double ahead = clearanceKm * std::max(0.0, params.aheadClearanceMult) +
                         std::max(0.0, params.aheadExtraKm);

    double alongKm = tExitKm + ahead;
    if (!std::isfinite(alongKm)) {
      return c;
    }

    alongKm = std::clamp(alongKm, 0.0, distGoal);

    const math::Vec3d side = sideRef * (double)sideSign;
    math::Vec3d w = posKm + dir * alongKm + side * clearanceKm;

    // Optional clamp to keep the waypoint from being absurdly far.
    const double maxW = std::max(0.0, params.maxWaypointDistKm);
    if (maxW > 1e-9) {
      const math::Vec3d toW = w - posKm;
      const double dW = toW.length();
      if (dW > maxW && dW > 1e-9) {
        w = posKm + toW * (maxW / dW);
      }
    }

    // Validate: segment pos->waypoint must be clear.
    const math::Vec3d toW = w - posKm;
    const double distW = toW.length();
    if (!std::isfinite(distW) || distW <= 1e-6) {
      return c;
    }

    const math::Vec3d dirW = toW / distW;
    const auto hit1 = field.raycastFirstHit(posKm, dirW, distW, queryPad, ignoreId, maxSphereTests);
    if (hit1.hit && hit1.tKm < distW - 1e-6) {
      return c;
    }

    // Validate: waypoint->goal should also be clear (helps avoid oscillation around the same obstacle).
    const math::Vec3d toG = goalKm - w;
    const double distG = toG.length();
    if (distG > 1e-6) {
      const math::Vec3d dirG = toG / distG;
      const auto hit2 = field.raycastFirstHit(w, dirG, distG, queryPad, ignoreId, maxSphereTests);
      if (hit2.hit && hit2.tKm < distG - 1e-6) {
        return c;
      }
    }

    const double ang = safeAcos(math::dot(dir, dirW));
    c.valid = true;
    c.waypointKm = w;

    // Prefer small steering angles. Add a tiny length term as deterministic tie-breaker.
    c.score = ang + 1e-6 * distW;
    return c;
  };

  // Normalise preferred side sign.
  if (preferredSideSign != -1 && preferredSideSign != 1) {
    preferredSideSign = 0;
  }

  Candidate best{};

  auto consider = [&](int s) {
    Candidate c = evalCandidate(s);
    if (!c.valid) return;
    if (!best.valid || c.score < best.score - 1e-12) {
      best = c;
      return;
    }
    if (std::abs(c.score - best.score) <= 1e-12) {
      // Tie-break deterministically: prefer the preferred side, else +1.
      if (preferredSideSign != 0) {
        if (c.sideSign == preferredSideSign && best.sideSign != preferredSideSign) {
          best = c;
          return;
        }
      } else if (c.sideSign > best.sideSign) {
        best = c;
        return;
      }
    }
  };

  if (preferredSideSign != 0) {
    consider(preferredSideSign);
    if (!best.valid) {
      consider(-preferredSideSign);
    }
  } else {
    consider(1);
    consider(-1);
  }

  if (!best.valid) {
    return out;
  }

  out.used = true;
  out.waypointKm = best.waypointKm;
  out.obstacleId = directHit.id;
  out.sideSign = best.sideSign;

  return out;
}

AvoidanceResultKm steerAvoidObstacles(const ProximityFieldKm& field,
                                     const math::Vec3d& posKm,
                                     const math::Vec3d& velKmS,
                                     const math::Vec3d& desiredDirUnitWorld,
                                     double desiredSpeedKmS,
                                     const AvoidanceParamsKm& params,
                                     int ignoreId) {
  AvoidanceResultKm out{};

  const auto desiredN = desiredDirUnitWorld.normalized();
  if (desiredN.lengthSq() <= 1e-12) {
    out.desiredDirUnit = {0, 0, 1};
    out.safeDirUnit = out.desiredDirUnit;
    return out;
  }

  out.desiredDirUnit = desiredN;
  out.safeDirUnit = desiredN;

  if (!params.enabled || field.empty()) {
    return out;
  }

  const double minSp = std::max(0.0, params.minSpeedForLookaheadKmS);
  const double sp = std::max({velKmS.length(), std::max(0.0, desiredSpeedKmS), minSp});

  const double lookaheadKm = std::max(0.0, params.lookaheadBaseKm) +
                             std::max(0.0, params.lookaheadTimeSec) * sp;
  out.lookaheadKm = lookaheadKm;

  if (lookaheadKm <= 1e-6) {
    return out;
  }

  const double inflateBase = std::max(0.0, params.shipRadiusKm) + std::max(0.0, params.padKm);
  const double nearMiss = std::max(0.0, params.nearMissExtraKm);
  const double queryPad = inflateBase + nearMiss;

  const math::Vec3d a = posKm;
  const math::Vec3d b = posKm + desiredN * lookaheadKm;
  math::Aabb3d q;
  q.min = math::Vec3d{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)} -
          math::Vec3d{queryPad, queryPad, queryPad};
  q.max = math::Vec3d{std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)} +
          math::Vec3d{queryPad, queryPad, queryPad};

  math::Vec3d repulse{0, 0, 0};

  int bestThreatId = -1;
  double bestThreatClear = std::numeric_limits<double>::infinity();
  double bestThreatAlong = 0.0;

  field.queryAabb(q, [&](int id) {
    if (id == ignoreId) return true;

    const auto* o = field.obstacle(id);
    if (!o) return true;
    if (o->radiusKm <= 0.0) return true;

    const double rInfl = o->radiusKm + inflateBase;

    const math::Vec3d toC = o->centerKm - posKm;
    const double along = math::dot(toC, desiredN);
    const double alongClamped = std::clamp(along, 0.0, lookaheadKm);
    const math::Vec3d closest = posKm + desiredN * alongClamped;

    const double distToCenter = (o->centerKm - closest).length();
    const double clearance = distToCenter - rInfl;

    if (clearance < bestThreatClear) {
      bestThreatClear = clearance;
      bestThreatId = id;
      bestThreatAlong = alongClamped;
    }

    // Influence region: within nearMiss of the inflated sphere surface.
    if (clearance > nearMiss) {
      return true;
    }

    const double hard = clamp01(o->hardness01);
    if (hard <= 0.0) {
      return true;
    }

    double w = 1.0;
    if (nearMiss > 1e-9) {
      w = clamp01((nearMiss - clearance) / nearMiss);
    } else {
      w = (clearance <= 0.0) ? 1.0 : 0.0;
    }

    w = w * w * hard;
    if (w <= 0.0) return true;

    math::Vec3d away = closest - o->centerKm;
    if (away.lengthSq() <= 1e-12) {
      away = stablePerp(desiredN);
    }

    const auto awayN = away.normalized();
    if (awayN.lengthSq() > 1e-12) {
      repulse = repulse + awayN * w;
    }

    return true;
  });

  if (bestThreatId >= 0) {
    out.threatId = bestThreatId;
    out.threatAlongKm = bestThreatAlong;
    out.threatClearanceKm = bestThreatClear;
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
