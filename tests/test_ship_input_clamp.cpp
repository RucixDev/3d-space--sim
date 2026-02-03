#include "test_harness.h"

#include "stellar/sim/Ship.h"

#include <cmath>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

static bool approxVec(const math::Vec3d& a, const math::Vec3d& b, double eps = 1e-9) {
  return approx(a.x, b.x, eps) && approx(a.y, b.y, eps) && approx(a.z, b.z, eps);
}

int test_ship_input_clamp() {
  int failures = 0;

  // We interpret maxLinAccel/maxAngAccel as *total* thruster authority caps.
  // Therefore, diagonal control inputs must be clamped to unit magnitude.

  // ---- Linear thrustLocal magnitude clamp ----
  {
    sim::Ship s;
    s.setOrientation({1,0,0,0});
    s.setVelocityKmS({0,0,0});
    s.setMaxLinearAccelKmS2(0.10);

    sim::ShipInput in{};
    in.dampers = false;
    in.brake = false;
    in.boost = false;
    in.thrustLocal = {1.0, 1.0, 1.0};

    s.step(1.0, in);

    const double inv = 1.0 / std::sqrt(3.0);
    const math::Vec3d expected{0.10 * inv, 0.10 * inv, 0.10 * inv};

    CHECK(approxVec(s.velocityKmS(), expected, 1e-9));
    CHECK(approx(s.velocityKmS().length(), 0.10, 1e-9));
  }

  // ---- Linear clamp after component clamp (e.g. >1 inputs) ----
  {
    sim::Ship s;
    s.setOrientation({1,0,0,0});
    s.setVelocityKmS({0,0,0});
    s.setMaxLinearAccelKmS2(0.20);

    sim::ShipInput in{};
    in.dampers = false;
    in.brake = false;
    in.boost = false;
    in.thrustLocal = {2.0, 2.0, 0.0}; // components will clamp to {1,1,0} then normalize

    s.step(1.0, in);

    const double inv = 1.0 / std::sqrt(2.0);
    const math::Vec3d expected{0.20 * inv, 0.20 * inv, 0.0};

    CHECK(approxVec(s.velocityKmS(), expected, 1e-9));
    CHECK(approx(s.velocityKmS().length(), 0.20, 1e-9));
  }

  // ---- Angular torqueLocal magnitude clamp ----
  {
    sim::Ship s;
    s.setOrientation({1,0,0,0});
    s.setAngularVelocityRadS({0,0,0});
    s.setMaxAngularAccelRadS2(0.60);

    sim::ShipInput in{};
    in.dampers = false;
    in.brake = false;
    in.boost = false;
    in.torqueLocal = {1.0, 1.0, 1.0};

    s.step(1.0, in);

    const double inv = 1.0 / std::sqrt(3.0);
    const math::Vec3d expected{0.60 * inv, 0.60 * inv, 0.60 * inv};

    CHECK(approxVec(s.angularVelocityRadS(), expected, 1e-9));
    CHECK(approx(s.angularVelocityRadS().length(), 0.60, 1e-9));
  }

  return failures;
}
