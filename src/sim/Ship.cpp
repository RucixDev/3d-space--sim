#include "stellar/sim/Ship.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

Ship::Ship() {
  // Default boost caps preserve legacy behavior (1.8x linear, 1.4x angular).
  // Gameplay code can override these via the boost setters.
  maxLinAccelBoostKmS2_ = maxLinAccelKmS2_ * 1.8;
  maxAngAccelBoostRadS2_ = maxAngAccelRadS2_ * 1.4;
}

double Ship::payloadHandlingScale() const {
  const double payload = (payloadMassKg_ > 0.0) ? payloadMassKg_ : 0.0;
  if (payload <= 1e-12) return 1.0;

  // If the gameplay layer provides an explicit payload capacity, use a simple
  // "cargo fraction" curve so cargo has a noticeable handling impact even when
  // the ship's dry mass is large.
  if (payloadCapacityKg_ > 1e-12) {
    const double frac = payload / payloadCapacityKg_;
    const double t = std::clamp(frac, 0.0, 1.0);

    const double minScale = std::clamp(payloadHandlingMinScale_, 0.05, 1.0);
    double scale = 1.0 - (1.0 - minScale) * t;

    // If overloaded, degrade further but softly (avoid "brick" behavior).
    if (frac > 1.0) {
      scale *= 1.0 / (1.0 + 0.5 * (frac - 1.0));
    }
    return std::clamp(scale, 0.05, 1.0);
  }

  // Fallback: physical scaling vs dry mass.
  const double dry = (massKg_ > 0.0) ? massKg_ : 0.0;
  const double total = dry + payload;
  if (!(total > 1e-12) || !(dry > 1e-12)) return 1.0;

  return std::clamp(dry / total, 0.05, 1.0);
}


void Ship::setMaxLinearAccelKmS2(double a) {
  maxLinAccelKmS2_ = a;
  if (!customLinBoost_) {
    maxLinAccelBoostKmS2_ = maxLinAccelKmS2_ * 1.8;
  }
}

void Ship::setMaxAngularAccelRadS2(double a) {
  maxAngAccelRadS2_ = a;
  if (!customAngBoost_) {
    maxAngAccelBoostRadS2_ = maxAngAccelRadS2_ * 1.4;
  }
}

void Ship::setMaxLinearAccelBoostKmS2(double a) {
  maxLinAccelBoostKmS2_ = a;
  customLinBoost_ = true;
}

void Ship::setMaxAngularAccelBoostRadS2(double a) {
  maxAngAccelBoostRadS2_ = a;
  customAngBoost_ = true;
}

static double decayKeff(double kPerSec, double dt) {
  if (!(kPerSec > 0.0) || !(dt > 0.0)) return 0.0;
  const double x = kPerSec * dt;
  // Use expm1 for numerical stability when x is tiny.
  const double oneMinus = -std::expm1(-x);
  return oneMinus / dt;
}

void Ship::step(double dtSeconds, const ShipInput& input) {
  stepWithExternalForces(dtSeconds, input, {0, 0, 0}, {0, 0, 0});
}

void Ship::stepWithExternalAccel(double dtSeconds,
                                 const ShipInput& input,
                                 const stellar::math::Vec3d& externalAccelWorldKmS2) {
  stepWithExternalForces(dtSeconds, input, externalAccelWorldKmS2, {0, 0, 0});
}

void Ship::stepWithExternalForces(double dtSeconds,
                                  const ShipInput& input,
                                  const stellar::math::Vec3d& externalAccelWorldKmS2,
                                  const stellar::math::Vec3d& externalAngAccelLocalRadS2) {
  if (dtSeconds <= 0.0) return;

  // Clamp user/control inputs defensively. (AI/autopilot may feed slightly out-of-range values.)
  ShipInput in = input;
  in.thrustLocal = stellar::math::clampComponents(in.thrustLocal, -1.0, 1.0);
  in.torqueLocal = stellar::math::clampComponents(in.torqueLocal, -1.0, 1.0);

  // Also clamp vector magnitude so diagonal inputs can't exceed the configured
  // acceleration caps. This keeps maxLinAccel/maxAngAccel semantics intuitive
  // (caps apply to total thruster authority, not per-axis sums).
  in.thrustLocal = stellar::math::clampMagnitude(in.thrustLocal, 1.0);
  in.torqueLocal = stellar::math::clampMagnitude(in.torqueLocal, 1.0);

  // Sub-step integration to keep the simple Euler-ish integrator stable under large dt.
  // This matters when the game uses time acceleration (timeScale) or frames hitch.
  constexpr double kMaxStep = 0.25;   // seconds
  constexpr int kMaxSteps = 4096;     // safety clamp

  int steps = static_cast<int>(std::ceil(dtSeconds / kMaxStep));
  steps = std::clamp(steps, 1, kMaxSteps);
  const double dt = dtSeconds / static_cast<double>(steps);

  for (int si = 0; si < steps; ++si) {
    // --------
    // Linear
    // --------
    const double linCap = (in.boost ? maxLinAccelBoostKmS2_ : maxLinAccelKmS2_) * payloadHandlingScale();

    stellar::math::Vec3d accelWorld = orient_.rotate(in.thrustLocal) * linCap;
    accelWorld += externalAccelWorldKmS2;

    if (in.dampers) {
      // Dampers attempt to kill velocity relative to the configured reference frame.
      //
      // Use a dt-invariant exponential decay in the unsaturated regime:
      //   v(t+dt) = v(t) * exp(-k*dt)
      // and convert it into an acceleration request that can still be capped by thruster authority.
      const stellar::math::Vec3d relVel = velKmS_ - dampingFrameVelKmS_;
      const double kEff = decayKeff(dampingLinear_, dt);
      const stellar::math::Vec3d damp = stellar::math::clampMagnitude(relVel * (-kEff), linCap);
      accelWorld += damp;
    }

    if (in.brake) {
      const double brakeCap = linCap * 2.0;
      const stellar::math::Vec3d relVel = velKmS_ - dampingFrameVelKmS_;
      const double kEff = decayKeff(dampingLinear_ * 6.0, dt);
      const stellar::math::Vec3d brake = stellar::math::clampMagnitude(relVel * (-kEff), brakeCap);
      accelWorld += brake;
    }

    velKmS_ += accelWorld * dt;
    posKm_ += velKmS_ * dt;

    // --------
    // Angular (body-local)
    // --------
    const double angCap = (in.boost ? maxAngAccelBoostRadS2_ : maxAngAccelRadS2_) * payloadHandlingScale();

    stellar::math::Vec3d angAccel = in.torqueLocal * angCap;
    angAccel += externalAngAccelLocalRadS2;

    if (in.dampers) {
      const double kEff = decayKeff(dampingAngular_, dt);
      const stellar::math::Vec3d dampW = stellar::math::clampMagnitude(angVelRadS_ * (-kEff), angCap);
      angAccel += dampW;
    }

    if (in.brake) {
      const double brakeCap = angCap * 2.0;
      const double kEff = decayKeff(dampingAngular_ * 6.0, dt);
      const stellar::math::Vec3d brakeW = stellar::math::clampMagnitude(angVelRadS_ * (-kEff), brakeCap);
      angAccel += brakeW;
    }

    angVelRadS_ += angAccel * dt;
    orient_ = orient_.integrateAngular(angVelRadS_, dt);
  }
}

} // namespace stellar::sim
