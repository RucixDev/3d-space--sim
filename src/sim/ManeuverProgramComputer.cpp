#include "stellar/sim/ManeuverProgramComputer.h"

#include "stellar/math/Math.h"
#include "stellar/sim/OrbitalMechanics.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

namespace {

static bool isFinite(const math::Vec3d& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

static bool isFinite(double v) {
  return std::isfinite(v);
}

static ManeuverProgramResult invalid(const GravityBody& refBody, ManeuverProgramKind program, const char* reason) {
  ManeuverProgramResult out{};
  out.valid = false;
  out.reason = reason ? reason : "";
  out.refBody = refBody;
  out.program = program;
  return out;
}

static double computeMu(const GravityBody& refBody, double gravityScale, const ManeuverProgramParams& params) {
  const double muBase = refBody.muKm3S2;
  return params.applyGravityScale ? (muBase * gravityScale) : muBase;
}

static math::Vec3d safeNormalize(const math::Vec3d& v, const math::Vec3d& fallback) {
  const double lsq = v.lengthSq();
  if (!(lsq > 1.0e-18) || !std::isfinite(lsq)) return fallback;
  return v / std::sqrt(lsq);
}

static math::Vec3d rotateRodrigues(const math::Vec3d& v, const math::Vec3d& axisUnit, double angleRad) {
  const double c = std::cos(angleRad);
  const double s = std::sin(angleRad);
  return v * c + math::cross(axisUnit, v) * s + axisUnit * (math::dot(axisUnit, v) * (1.0 - c));
}

// Align the orbit plane at a specific node state.
// Returns false on degenerate input.
static bool computePlaneAlignDvAtNode(const math::Vec3d& relPosKm,
                                      const math::Vec3d& relVelKmS,
                                      bool forcePrograde,
                                      math::Vec3d& outDvKmS) {
  outDvKmS = {0, 0, 0};

  const math::Vec3d h = math::cross(relPosKm, relVelKmS);
  const double h2 = h.lengthSq();
  if (!(h2 > 1e-18)) return false;

  const math::Vec3d hHat = h * (1.0 / std::sqrt(h2));

  math::Vec3d targetNormal{0, 0, 1};
  if (!forcePrograde && hHat.z < 0.0) targetNormal = {0, 0, -1};

  const double dotN = std::clamp(math::dot(hHat, targetNormal), -1.0, 1.0);
  const double angle = std::acos(dotN);
  if (!(angle >= 0.0) || !std::isfinite(angle)) return false;
  if (angle < 1e-9) {
    outDvKmS = {0, 0, 0};
    return true;
  }

  math::Vec3d rHat = safeNormalize(relPosKm, {1, 0, 0});

  // Choose the sign that best aligns the rotated orbit normal to the target normal.
  const math::Vec3d hPlus = rotateRodrigues(hHat, rHat, angle);
  const math::Vec3d hMinus = rotateRodrigues(hHat, rHat, -angle);
  const double dp = math::dot(hPlus, targetNormal);
  const double dm = math::dot(hMinus, targetNormal);
  const double signedAngle = (dp >= dm) ? angle : -angle;

  const math::Vec3d vDesired = rotateRodrigues(relVelKmS, rHat, signedAngle);
  outDvKmS = vDesired - relVelKmS;
  return true;
}

static double wrap2pi(double a) {
  const double twoPi = 2.0 * math::kPi;
  a = std::fmod(a, twoPi);
  if (a < 0.0) a += twoPi;
  return a;
}

} // namespace

ManeuverProgramResult planCircularize(const Ship& ship,
                                      double nowTimeDays,
                                      const GravityBody& refBody,
                                      ManeuverProgramKind kind,
                                      double gravityScale,
                                      const ManeuverProgramParams& params) {
  ManeuverProgramResult out{};
  out.refBody = refBody;
  out.program = kind;

  const double mu = computeMu(refBody, gravityScale, params);
  if (!(mu > 0.0) || !isFinite(mu)) {
    return invalid(refBody, kind, "Invalid mu");
  }

  const math::Vec3d relPos = ship.positionKm() - refBody.posKm;
  const math::Vec3d relVel = ship.velocityKmS() - refBody.velKmS;
  if (!isFinite(relPos) || !isFinite(relVel)) {
    return invalid(refBody, kind, "Non-finite ship state");
  }

  const auto el = solveClassicalOrbitElements(relPos, relVel, mu);
  if (!el.valid) {
    return invalid(refBody, kind, "Invalid orbit (degenerate state vector)");
  }
  if (el.type != TwoBodyOrbit::Type::Elliptic) {
    return invalid(refBody, kind, "Orbit is not elliptic (circularization not supported)");
  }

  const double targetNu = (kind == ManeuverProgramKind::CircularizeAtApoapsis) ? math::kPi : 0.0;
  double timeToNodeSec = timeToTrueAnomalySecElliptic(el, targetNu);
  if (!(timeToNodeSec >= 0.0) || !isFinite(timeToNodeSec)) {
    return invalid(refBody, kind, "Failed to compute time-to-node");
  }

  // State at target anomaly.
  math::Vec3d posNode{0, 0, 0};
  math::Vec3d velNode{0, 0, 0};
  if (!stateFromClassicalOrbitElementsAtTrueAnomaly(el, targetNu, posNode, velNode)) {
    return invalid(refBody, kind, "Failed to evaluate node state");
  }

  double r = posNode.length();
  if (!(r > 1e-9) || !isFinite(r)) {
    return invalid(refBody, kind, "Invalid node radius");
  }

  double vNode = velNode.length();
  double vCirc = std::sqrt(mu / r);
  if (!isFinite(vNode) || !isFinite(vCirc)) {
    return invalid(refBody, kind, "Invalid node speed");
  }

  // At apoapsis/periapsis the velocity is (approximately) tangential. Use the
  // direction of the Keplerian node velocity to determine burn direction.
  if (!(vNode > 1e-12)) {
    return invalid(refBody, kind, "Zero node velocity");
  }

  const math::Vec3d dvDir = velNode * (1.0 / vNode);
  const double dvSigned = vCirc - vNode;
  math::Vec3d dv = dvDir * dvSigned;

  const double dvMag = dv.length();

  // If we're already circular (within tolerance), return a no-op plan at now.
  // Circular orbits are degenerate (peri/apo not unique), so it is more useful
  // to report an immediate "no maneuver required".
  if (dvMag < params.minDvKmS) {
    dv = {0, 0, 0};
    timeToNodeSec = 0.0;

    // Use the current state for node diagnostics.
    posNode = relPos;
    velNode = relVel;
    r = posNode.length();
    vNode = velNode.length();
    if (r > 1e-9) {
      vCirc = std::sqrt(mu / r);
    }
  }

  out.valid = true;
  out.timeToNodeSec = timeToNodeSec;
  out.dvKmS = dv.length();
  out.nodeRadiusKm = r;
  out.nodeAltitudeKm = r - refBody.radiusKm;
  out.nodeSpeedKmS = vNode;

  out.targetRadiusKm = r;
  out.targetSpeedKmS = vCirc;
  out.targetCircularSpeedKmS = vCirc;

  out.nodeRelPosKm = posNode;
  out.nodeRelVelKmS = velNode;

  out.plan.nodeTimeDays = nowTimeDays + timeToNodeSec / 86400.0;
  out.plan.deltaVWorldKmS = dv;

  return out;
}

ManeuverProgramResult planSetApoapsisAtPeriapsis(const Ship& ship,
                                                 double nowTimeDays,
                                                 const GravityBody& refBody,
                                                 double targetApoapsisKm,
                                                 double gravityScale,
                                                 const ManeuverProgramParams& params) {
  const ManeuverProgramKind program = ManeuverProgramKind::SetApoapsisAtPeriapsis;

  const double mu = computeMu(refBody, gravityScale, params);
  if (!(mu > 0.0) || !isFinite(mu)) {
    return invalid(refBody, program, "Invalid mu");
  }

  const math::Vec3d relPos = ship.positionKm() - refBody.posKm;
  const math::Vec3d relVel = ship.velocityKmS() - refBody.velKmS;
  if (!isFinite(relPos) || !isFinite(relVel)) {
    return invalid(refBody, program, "Non-finite ship state");
  }

  const auto el = solveClassicalOrbitElements(relPos, relVel, mu);
  if (!el.valid) {
    return invalid(refBody, program, "Invalid orbit (degenerate state vector)");
  }
  if (el.type != TwoBodyOrbit::Type::Elliptic) {
    return invalid(refBody, program, "Orbit is not elliptic (apsis change not supported)");
  }

  const double targetNu = 0.0; // periapsis
  double timeToNodeSec = timeToTrueAnomalySecElliptic(el, targetNu);
  if (!(timeToNodeSec >= 0.0) || !isFinite(timeToNodeSec)) {
    return invalid(refBody, program, "Failed to compute time-to-periapsis");
  }

  math::Vec3d posNode{0, 0, 0};
  math::Vec3d velNode{0, 0, 0};
  if (!stateFromClassicalOrbitElementsAtTrueAnomaly(el, targetNu, posNode, velNode)) {
    return invalid(refBody, program, "Failed to evaluate periapsis state");
  }

  const double rNode = posNode.length();
  if (!(rNode > 1e-9) || !isFinite(rNode)) {
    return invalid(refBody, program, "Invalid periapsis radius");
  }

  const double vNode = velNode.length();
  if (!(vNode > 1e-12) || !isFinite(vNode)) {
    return invalid(refBody, program, "Invalid periapsis speed");
  }

  const double minR = std::max(0.0, refBody.radiusKm);
  double raTarget = targetApoapsisKm;
  if (!isFinite(raTarget)) raTarget = rNode;
  raTarget = std::max(raTarget, minR);
  raTarget = std::max(raTarget, rNode); // apo must be >= peri

  const double aNew = 0.5 * (rNode + raTarget);
  if (!(aNew > 1e-9) || !isFinite(aNew)) {
    return invalid(refBody, program, "Invalid target semi-major axis");
  }

  const double vDesired = std::sqrt(mu * (2.0 / rNode - 1.0 / aNew));
  if (!isFinite(vDesired)) {
    return invalid(refBody, program, "Invalid target speed");
  }

  const math::Vec3d dvDir = velNode * (1.0 / vNode);
  math::Vec3d dv = dvDir * (vDesired - vNode);
  if (dv.length() < params.minDvKmS) {
    dv = {0, 0, 0};
  }

  ManeuverProgramResult out{};
  out.valid = true;
  out.program = program;
  out.refBody = refBody;

  out.timeToNodeSec = timeToNodeSec;
  out.dvKmS = dv.length();
  out.nodeRadiusKm = rNode;
  out.nodeAltitudeKm = rNode - refBody.radiusKm;
  out.nodeSpeedKmS = vNode;

  out.targetRadiusKm = raTarget;
  out.targetSpeedKmS = vDesired;

  out.nodeRelPosKm = posNode;
  out.nodeRelVelKmS = velNode;

  out.plan.nodeTimeDays = nowTimeDays + timeToNodeSec / 86400.0;
  out.plan.deltaVWorldKmS = dv;

  return out;
}

ManeuverProgramResult planSetPeriapsisAtApoapsis(const Ship& ship,
                                                 double nowTimeDays,
                                                 const GravityBody& refBody,
                                                 double targetPeriapsisKm,
                                                 double gravityScale,
                                                 const ManeuverProgramParams& params) {
  const ManeuverProgramKind program = ManeuverProgramKind::SetPeriapsisAtApoapsis;

  const double mu = computeMu(refBody, gravityScale, params);
  if (!(mu > 0.0) || !isFinite(mu)) {
    return invalid(refBody, program, "Invalid mu");
  }

  const math::Vec3d relPos = ship.positionKm() - refBody.posKm;
  const math::Vec3d relVel = ship.velocityKmS() - refBody.velKmS;
  if (!isFinite(relPos) || !isFinite(relVel)) {
    return invalid(refBody, program, "Non-finite ship state");
  }

  const auto el = solveClassicalOrbitElements(relPos, relVel, mu);
  if (!el.valid) {
    return invalid(refBody, program, "Invalid orbit (degenerate state vector)");
  }
  if (el.type != TwoBodyOrbit::Type::Elliptic) {
    return invalid(refBody, program, "Orbit is not elliptic (apsis change not supported)");
  }

  const double targetNu = math::kPi; // apoapsis
  double timeToNodeSec = timeToTrueAnomalySecElliptic(el, targetNu);
  if (!(timeToNodeSec >= 0.0) || !isFinite(timeToNodeSec)) {
    return invalid(refBody, program, "Failed to compute time-to-apoapsis");
  }

  math::Vec3d posNode{0, 0, 0};
  math::Vec3d velNode{0, 0, 0};
  if (!stateFromClassicalOrbitElementsAtTrueAnomaly(el, targetNu, posNode, velNode)) {
    return invalid(refBody, program, "Failed to evaluate apoapsis state");
  }

  const double rNode = posNode.length();
  if (!(rNode > 1e-9) || !isFinite(rNode)) {
    return invalid(refBody, program, "Invalid apoapsis radius");
  }

  const double vNode = velNode.length();
  if (!(vNode > 1e-12) || !isFinite(vNode)) {
    return invalid(refBody, program, "Invalid apoapsis speed");
  }

  const double minR = std::max(0.0, refBody.radiusKm);
  double rpTarget = targetPeriapsisKm;
  if (!isFinite(rpTarget)) rpTarget = rNode;

  rpTarget = std::max(rpTarget, minR);
  rpTarget = std::min(rpTarget, rNode); // peri must be <= apo

  const double aNew = 0.5 * (rNode + rpTarget);
  if (!(aNew > 1e-9) || !isFinite(aNew)) {
    return invalid(refBody, program, "Invalid target semi-major axis");
  }

  const double vDesired = std::sqrt(mu * (2.0 / rNode - 1.0 / aNew));
  if (!isFinite(vDesired)) {
    return invalid(refBody, program, "Invalid target speed");
  }

  const math::Vec3d dvDir = velNode * (1.0 / vNode);
  math::Vec3d dv = dvDir * (vDesired - vNode);
  if (dv.length() < params.minDvKmS) {
    dv = {0, 0, 0};
  }

  ManeuverProgramResult out{};
  out.valid = true;
  out.program = program;
  out.refBody = refBody;

  out.timeToNodeSec = timeToNodeSec;
  out.dvKmS = dv.length();
  out.nodeRadiusKm = rNode;
  out.nodeAltitudeKm = rNode - refBody.radiusKm;
  out.nodeSpeedKmS = vNode;

  out.targetRadiusKm = rpTarget;
  out.targetSpeedKmS = vDesired;

  out.nodeRelPosKm = posNode;
  out.nodeRelVelKmS = velNode;

  out.plan.nodeTimeDays = nowTimeDays + timeToNodeSec / 86400.0;
  out.plan.deltaVWorldKmS = dv;

  return out;
}

ManeuverProgramResult planEscapeNow(const Ship& ship,
                                    double nowTimeDays,
                                    const GravityBody& refBody,
                                    double gravityScale,
                                    const ManeuverProgramParams& params) {
  const ManeuverProgramKind program = ManeuverProgramKind::EscapeNow;

  const double mu = computeMu(refBody, gravityScale, params);
  if (!(mu > 0.0) || !isFinite(mu)) {
    return invalid(refBody, program, "Invalid mu");
  }

  const math::Vec3d relPos = ship.positionKm() - refBody.posKm;
  const math::Vec3d relVel = ship.velocityKmS() - refBody.velKmS;
  if (!isFinite(relPos) || !isFinite(relVel)) {
    return invalid(refBody, program, "Non-finite ship state");
  }

  const double r = relPos.length();
  if (!(r > 1e-9) || !isFinite(r)) {
    return invalid(refBody, program, "Invalid radius");
  }

  const double v = relVel.length();
  if (!isFinite(v)) {
    return invalid(refBody, program, "Invalid speed");
  }

  const double vEsc = std::sqrt(2.0 * mu / r);
  if (!isFinite(vEsc)) {
    return invalid(refBody, program, "Invalid escape speed");
  }

  double dvMag = vEsc - v;
  if (dvMag < params.minDvKmS) dvMag = 0.0;

  const math::Vec3d dvDir = (v > 1e-12) ? (relVel * (1.0 / v)) : math::Vec3d{0, 0, 1};
  const math::Vec3d dv = dvDir * std::max(0.0, dvMag);

  ManeuverProgramResult out{};
  out.valid = true;
  out.program = program;
  out.refBody = refBody;

  out.timeToNodeSec = 0.0;
  out.dvKmS = dv.length();
  out.nodeRadiusKm = r;
  out.nodeAltitudeKm = r - refBody.radiusKm;
  out.nodeSpeedKmS = v;

  out.targetRadiusKm = r;
  out.targetSpeedKmS = vEsc;

  out.nodeRelPosKm = relPos;
  out.nodeRelVelKmS = relVel;

  out.plan.nodeTimeDays = nowTimeDays;
  out.plan.deltaVWorldKmS = dv;

  return out;
}

ManeuverProgramResult planAlignPlaneAtNode(const Ship& ship,
                                           double nowTimeDays,
                                           const GravityBody& refBody,
                                           bool ascendingNode,
                                           bool forcePrograde,
                                           double gravityScale,
                                           const ManeuverProgramParams& params) {
  const ManeuverProgramKind program = ascendingNode ?
    ManeuverProgramKind::AlignPlaneAtAscendingNode :
    ManeuverProgramKind::AlignPlaneAtDescendingNode;

  const double mu = computeMu(refBody, gravityScale, params);
  if (!(mu > 0.0) || !isFinite(mu)) {
    return invalid(refBody, program, "Invalid mu");
  }

  const math::Vec3d relPos = ship.positionKm() - refBody.posKm;
  const math::Vec3d relVel = ship.velocityKmS() - refBody.velKmS;
  if (!isFinite(relPos) || !isFinite(relVel)) {
    return invalid(refBody, program, "Non-finite ship state");
  }

  const auto el = solveClassicalOrbitElements(relPos, relVel, mu);
  if (!el.valid) {
    return invalid(refBody, program, "Invalid orbit (degenerate state vector)");
  }
  if (el.type != TwoBodyOrbit::Type::Elliptic) {
    return invalid(refBody, program, "Orbit is not elliptic (plane align not supported)");
  }
  if (el.equatorial) {
    return invalid(refBody, program, "Orbit is equatorial (plane already aligned)");
  }

  // Ascending node: argument of latitude u = ω + ν = 0.
  // Descending node: u = π.
  const double nuNode = ascendingNode ? wrap2pi(-el.argPeriapsisRad)
                                      : wrap2pi(math::kPi - el.argPeriapsisRad);

  const double timeToNodeSec = timeToTrueAnomalySecElliptic(el, nuNode);
  if (!(timeToNodeSec >= 0.0) || !isFinite(timeToNodeSec)) {
    return invalid(refBody, program, "Failed to compute time-to-node");
  }

  math::Vec3d posNode{0, 0, 0};
  math::Vec3d velNode{0, 0, 0};
  if (!stateFromClassicalOrbitElementsAtTrueAnomaly(el, nuNode, posNode, velNode)) {
    return invalid(refBody, program, "Failed to evaluate node state");
  }

  const double rNode = posNode.length();
  if (!(rNode > 1e-9) || !isFinite(rNode)) {
    return invalid(refBody, program, "Invalid node radius");
  }

  const double vNode = velNode.length();
  if (!(vNode > 1e-12) || !isFinite(vNode)) {
    return invalid(refBody, program, "Invalid node speed");
  }

  math::Vec3d dv{0, 0, 0};
  if (!computePlaneAlignDvAtNode(posNode, velNode, forcePrograde, dv)) {
    return invalid(refBody, program, "Failed to compute plane-change dv");
  }

  if (dv.length() < params.minDvKmS) {
    dv = {0, 0, 0};
  }

  ManeuverProgramResult out{};
  out.valid = true;
  out.program = program;
  out.refBody = refBody;

  out.timeToNodeSec = timeToNodeSec;
  out.dvKmS = dv.length();
  out.nodeRadiusKm = rNode;
  out.nodeAltitudeKm = rNode - refBody.radiusKm;
  out.nodeSpeedKmS = vNode;

  out.targetRadiusKm = rNode;
  out.targetSpeedKmS = vNode;

  out.nodeRelPosKm = posNode;
  out.nodeRelVelKmS = velNode;

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
    out.program = kind;
    return out;
  }
  return planCircularize(ship, nowTimeDays, dom.body, kind, gravityParams.scale, params);
}

} // namespace stellar::sim
