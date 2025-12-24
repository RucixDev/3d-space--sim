#include "stellar/sim/Ship.h"

#include "stellar/math/Math.h"

namespace stellar::sim {

Ship::Ship() = default;

static stellar::math::Vec3d clampMagnitude(const stellar::math::Vec3d& v, double maxLen) {
  const double len = v.length();
  if (len <= maxLen || len <= 1e-12) return v;
  return v * (maxLen / len);
}

void Ship::step(double dtSeconds, const ShipInput& input) {
  if (dtSeconds <= 0.0) return;

  // --------
  // Linear
  // --------
  double linCap = maxLinAccelKmS2_;
  if (input.boost) linCap *= 1.8;

  stellar::math::Vec3d accelWorld = orient_.rotate(input.thrustLocal) * linCap;

  if (input.dampers) {
    // Dampers attempt to kill velocity (uses thrusters, so cap it).
    const stellar::math::Vec3d damp = clampMagnitude(velKmS_ * (-dampingLinear_), linCap);
    accelWorld += damp;
  }

  if (input.brake) {
    const double brakeCap = linCap * 2.0;
    const stellar::math::Vec3d brake = clampMagnitude(velKmS_ * (-dampingLinear_ * 6.0), brakeCap);
    accelWorld += brake;
  }

  velKmS_ += accelWorld * dtSeconds;
  posKm_ += velKmS_ * dtSeconds;

  // --------
  // Angular (body-local)
  // --------
  double angCap = maxAngAccelRadS2_;
  if (input.boost) angCap *= 1.4;

  stellar::math::Vec3d angAccel = input.torqueLocal * angCap;

  if (input.dampers) {
    const stellar::math::Vec3d dampW = clampMagnitude(angVelRadS_ * (-dampingAngular_), angCap);
    angAccel += dampW;
  }

  if (input.brake) {
    const double brakeCap = angCap * 2.0;
    const stellar::math::Vec3d brakeW = clampMagnitude(angVelRadS_ * (-dampingAngular_ * 6.0), brakeCap);
    angAccel += brakeW;
  }

  angVelRadS_ += angAccel * dtSeconds;
  orient_ = orient_.integrateAngular(angVelRadS_, dtSeconds);
}

void Ship::stepAngularOnly(double dtSeconds, const ShipInput& input) {
  if (dtSeconds <= 0.0) return;

  // --------
  // Angular (body-local)
  // --------
  double angCap = maxAngAccelRadS2_;
  if (input.boost) angCap *= 1.4;

  stellar::math::Vec3d angAccel = input.torqueLocal * angCap;

  if (input.dampers) {
    const stellar::math::Vec3d dampW = clampMagnitude(angVelRadS_ * (-dampingAngular_), angCap);
    angAccel += dampW;
  }

  if (input.brake) {
    const double brakeCap = angCap * 2.0;
    const stellar::math::Vec3d brakeW = clampMagnitude(angVelRadS_ * (-dampingAngular_ * 6.0), brakeCap);
    angAccel += brakeW;
  }

  angVelRadS_ += angAccel * dtSeconds;
  orient_ = orient_.integrateAngular(angVelRadS_, dtSeconds);
}

} // namespace stellar::sim
