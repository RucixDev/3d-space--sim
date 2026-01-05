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

static stellar::math::Vec3d clampMagnitude(const stellar::math::Vec3d& v, double maxLen) {
  const double len = v.length();
  if (len <= maxLen || len <= 1e-12) return v;
  return v * (maxLen / len);
}

static stellar::math::Vec3d clampComponents(const stellar::math::Vec3d& v, double lo, double hi) {
  return {
    std::clamp(v.x, lo, hi),
    std::clamp(v.y, lo, hi),
    std::clamp(v.z, lo, hi),
  };
}

void Ship::step(double dtSeconds, const ShipInput& input) {
  stepWithExternalAccel(dtSeconds, input, {0,0,0});
}

void Ship::stepWithExternalAccel(double dtSeconds,
                                 const ShipInput& input,
                                 const stellar::math::Vec3d& externalAccelWorldKmS2) {
  if (dtSeconds <= 0.0) return;

  // Clamp user/control inputs defensively. (AI/autopilot may feed slightly out-of-range values.)
  ShipInput in = input;
  in.thrustLocal = clampComponents(in.thrustLocal, -1.0, 1.0);
  in.torqueLocal = clampComponents(in.torqueLocal, -1.0, 1.0);

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
    const double linCap = in.boost ? maxLinAccelBoostKmS2_ : maxLinAccelKmS2_;

    stellar::math::Vec3d accelWorld = orient_.rotate(in.thrustLocal) * linCap;
    accelWorld += externalAccelWorldKmS2;

    if (in.dampers) {
      // Dampers attempt to kill velocity (uses thrusters, so cap it).
      const stellar::math::Vec3d relVel = velKmS_ - dampingFrameVelKmS_;
      const stellar::math::Vec3d damp = clampMagnitude(relVel * (-dampingLinear_), linCap);
      accelWorld += damp;
    }

    if (in.brake) {
      const double brakeCap = linCap * 2.0;
      const stellar::math::Vec3d relVel = velKmS_ - dampingFrameVelKmS_;
      const stellar::math::Vec3d brake = clampMagnitude(relVel * (-dampingLinear_ * 6.0), brakeCap);
      accelWorld += brake;
    }

    velKmS_ += accelWorld * dt;
    posKm_ += velKmS_ * dt;

    // --------
    // Angular (body-local)
    // --------
    const double angCap = in.boost ? maxAngAccelBoostRadS2_ : maxAngAccelRadS2_;

    stellar::math::Vec3d angAccel = in.torqueLocal * angCap;

    if (in.dampers) {
      const stellar::math::Vec3d dampW = clampMagnitude(angVelRadS_ * (-dampingAngular_), angCap);
      angAccel += dampW;
    }

    if (in.brake) {
      const double brakeCap = angCap * 2.0;
      const stellar::math::Vec3d brakeW = clampMagnitude(angVelRadS_ * (-dampingAngular_ * 6.0), brakeCap);
      angAccel += brakeW;
    }

    angVelRadS_ += angAccel * dt;
    orient_ = orient_.integrateAngular(angVelRadS_, dt);
  }
}

} // namespace stellar::sim
