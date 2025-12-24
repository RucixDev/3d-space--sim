#pragma once

#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"

namespace stellar::sim {

// Units for ship dynamics:
// - position: km (kilometers) in local system frame
// - linear velocity: km/s
// - angular velocity: rad/s
struct ShipInput {
  // Translation thruster intent in body-local axes:
  //  X = right, Y = up, Z = forward
  // Range typically [-1,1].
  math::Vec3d thrustLocal{0,0,0};

  // Rotation thruster intent in body-local axes (pitch/yaw/roll):
  //  X = pitch (nose up positive), Y = yaw (turn right positive), Z = roll (right wing down positive)
  // Range typically [-1,1].
  math::Vec3d torqueLocal{0,0,0};

  bool boost{false};
  bool brake{false};
  bool dampers{true};
};

class Ship {
public:
  Ship();

  // State
  math::Vec3d positionKm() const { return posKm_; }
  math::Vec3d velocityKmS() const { return velKmS_; }
  math::Quatd orientation() const { return orient_; }
  math::Vec3d angularVelocityRadS() const { return angVelRadS_; }

  void setPositionKm(const math::Vec3d& p) { posKm_ = p; }
  void setVelocityKmS(const math::Vec3d& v) { velKmS_ = v; }
  void setOrientation(const math::Quatd& q) { orient_ = q.normalized(); }
  void setAngularVelocityRadS(const math::Vec3d& w) { angVelRadS_ = w; }

  // Params
  double massKg() const { return massKg_; }
  math::Vec3d inertiaDiagKgKm2() const { return inertiaDiagKgKm2_; }

  void setMassKg(double m) { massKg_ = m; }
  void setInertiaDiagKgKm2(const math::Vec3d& i) { inertiaDiagKgKm2_ = i; }

  // Thruster limits (force in kN translated into km/s^2 via mass)
  void setMaxLinearAccelKmS2(double a) { maxLinAccelKmS2_ = a; }
  void setMaxAngularAccelRadS2(double a) { maxAngAccelRadS2_ = a; }

  double maxLinearAccelKmS2() const { return maxLinAccelKmS2_; }
  double maxAngularAccelRadS2() const { return maxAngAccelRadS2_; }

  void setDampingLinear(double d) { dampingLinear_ = d; }     // per-second
  void setDampingAngular(double d) { dampingAngular_ = d; }   // per-second

  double dampingLinear() const { return dampingLinear_; }
  double dampingAngular() const { return dampingAngular_; }

  // Update physics
  void step(double dtSeconds, const ShipInput& input);

  // Convenience vectors in world space
  math::Vec3d forward() const { return orient_.rotate({0,0,1}); }
  math::Vec3d right() const { return orient_.rotate({1,0,0}); }
  math::Vec3d up() const { return orient_.rotate({0,1,0}); }

private:
  // Mass properties
  double massKg_{1.0e5}; // 100t
  // For simplicity: diagonal inertia in (kg*km^2) (km to keep magnitudes reasonable).
  math::Vec3d inertiaDiagKgKm2_{ 2.0e4, 2.0e4, 3.0e4 };

  // Limits
  double maxLinAccelKmS2_{0.05};     // km/s^2 (~50 m/s^2) (arcade-ish but usable)
  double maxAngAccelRadS2_{0.8};     // rad/s^2

  // Dampers (simple exponential decay)
  double dampingLinear_{0.15};       // 1/s
  double dampingAngular_{0.25};      // 1/s

  // State
  math::Vec3d posKm_{0,0,0};
  math::Vec3d velKmS_{0,0,0};
  math::Quatd orient_{1,0,0,0};
  math::Vec3d angVelRadS_{0,0,0};
};

} // namespace stellar::sim
