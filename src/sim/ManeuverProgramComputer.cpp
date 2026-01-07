#include "stellar/sim/ManeuverProgramComputer.h"

#include "stellar/math/Math.h"
#include "stellar/sim/OrbitalMechanics.h"

#include <cmath>

namespace stellar::sim {

namespace {

static bool isFinite(const math::Vec3d& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

static bool isFinite(double v) { return std::isfinite(v); }

} // namespace

ManeuverProgramResult planCircularize(const Ship& ship,
                                      double nowTimeDays,
                                      const GravityBody& refBody,
                                      ManeuverProgramKind kind,
                                      double gravityScale,
                                      const ManeuverProgramParams& params) {
  ManeuverProgramResult out{};
  out.refBody = refBody;

  const double muBase = refBody.muKm3S2;
  const double mu = (params.applyGravityScale ? (muBase * gravityScale) : muBase);
  if (!(mu > 0.0) || !isFinite(mu)) {
    out.valid = false;
    out.reason = "Invalid mu";
    return out;
  }

  const math::Vec3d relPos = ship.positionKm() - refBody.posKm;
  const math::Vec3d relVel = ship.velocityKmS() - refBody.velKmS;
  if (!isFinite(relPos) || !isFinite(relVel)) {
    out.valid = false;
    out.reason = "Non-finite ship state";
    return out;
  }

  const auto el = solveClassicalOrbitElements(relPos, relVel, mu);
  if (!el.valid) {
    out.valid = false;
    out.reason = "Invalid orbit (degenerate state vector)";
    return out;
  }
  if (el.type != TwoBodyOrbit::Type::Elliptic) {
    out.valid = false;
    out.reason = "Orbit is not elliptic (circularization not supported)";
    return out;
  }

  const double targetNu = (kind == ManeuverProgramKind::CircularizeAtApoapsis) ? math::kPi : 0.0;
  double timeToNodeSec = timeToTrueAnomalySecElliptic(el, targetNu);
  if (!(timeToNodeSec >= 0.0) || !isFinite(timeToNodeSec)) {
    out.valid = false;
    out.reason = "Failed to compute time-to-node";
    return out;
  }

  // State at target anomaly.
  math::Vec3d posNode{0, 0, 0};
  math::Vec3d velNode{0, 0, 0};
  if (!stateFromClassicalOrbitElementsAtTrueAnomaly(el, targetNu, posNode, velNode)) {
    out.valid = false;
    out.reason = "Failed to evaluate node state";
    return out;
  }

  const double r = posNode.length();
  if (!(r > 1e-9) || !isFinite(r)) {
    out.valid = false;
    out.reason = "Invalid node radius";
    return out;
  }

  const double vNode = velNode.length();
  const double vCirc = std::sqrt(mu / r);
  if (!isFinite(vNode) || !isFinite(vCirc)) {
    out.valid = false;
    out.reason = "Invalid node speed";
    return out;
  }

  // At apoapsis/periapsis the velocity is (approximately) tangential. Use the
  // direction of the Keplerian node velocity to determine burn direction.
  math::Vec3d dvDir{0, 0, 0};
  if (vNode > 1e-12) {
    dvDir = velNode * (1.0 / vNode);
  } else {
    // Degenerate; shouldn't happen for a valid elliptic orbit.
    out.valid = false;
    out.reason = "Zero node velocity";
    return out;
  }

  const double dvSigned = vCirc - vNode;
  math::Vec3d dv = dvDir * dvSigned;
  const double dvMag = dv.length();

  // If we're already circular (within tolerance), return a no-op plan at now.
  if (dvMag < params.minDvKmS) {
    dv = {0, 0, 0};
    timeToNodeSec = 0.0;
  }

  out.valid = true;
  out.timeToNodeSec = timeToNodeSec;
  out.dvKmS = dv.length();
  out.nodeRadiusKm = r;
  out.nodeAltitudeKm = r - refBody.radiusKm;
  out.nodeSpeedKmS = vNode;
  out.targetCircularSpeedKmS = vCirc;

  out.plan.nodeTimeDays = nowTimeDays + timeToNodeSec / 86400.0;
  out.plan.deltaVWorldKmS = dv;

  return out;
}

ManeuverProgramResult planCircularizeDominantBody(const StarSystem& sys,
                                                  double nowTimeDays,
                                                  const Ship& ship,
                                                  const GravityParams& gravityParams,
                                                  ManeuverProgramKind kind,
                                                  const ManeuverProgramParams& params) {
  const auto dom = dominantGravityBody(sys, nowTimeDays, ship.positionKm(), gravityParams);
  if (!dom.valid) {
    ManeuverProgramResult out;
    out.valid = false;
    out.reason = "No gravity bodies enabled";
    return out;
  }
  return planCircularize(ship, nowTimeDays, dom.body, kind, gravityParams.scale, params);
}

} // namespace stellar::sim