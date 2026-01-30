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

static sim::Ship makeTestShip() {
  sim::Ship ship;
  ship.setPositionKm({0.0, 0.0, 0.0});
  ship.setVelocityKmS({0.0, 0.0, 0.0});
  ship.setOrientation({1, 0, 0, 0});
  ship.setAngularVelocityRadS({0.0, 0.0, 0.0});

  // Very high caps so damping runs in the unsaturated regime.
  ship.setMaxLinearAccelKmS2(1.0e6);
  ship.setMaxAngularAccelRadS2(1.0e6);

  // Default-ish damping rates (per second).
  ship.setDampingLinear(0.7);
  ship.setDampingAngular(1.1);

  ship.setDampingFrameVelocityKmS({0.0, 0.0, 0.0});
  return ship;
}

int test_ship_dampers() {
  int failures = 0;

  sim::ShipInput in{};
  in.thrustLocal = {0.0, 0.0, 0.0};
  in.torqueLocal = {0.0, 0.0, 0.0};
  in.dampers = true;
  in.brake = false;
  in.boost = false;

  // ---- Linear dampers should be dt-invariant (one big step == many small steps). ----
  {
    sim::Ship a = makeTestShip();
    a.setVelocityKmS({10.0, -5.0, 2.0});
    a.step(0.1, in);
    const auto vA = a.velocityKmS();

    sim::Ship b = makeTestShip();
    b.setVelocityKmS({10.0, -5.0, 2.0});
    for (int i = 0; i < 10; ++i) {
      b.step(0.01, in);
    }
    const auto vB = b.velocityKmS();

    CHECK(approxVec(vA, vB, 1e-9));
  }

  // ---- Angular dampers should be dt-invariant. ----
  {
    sim::Ship a = makeTestShip();
    a.setAngularVelocityRadS({1.0, 2.0, -0.5});
    a.step(0.1, in);
    const auto wA = a.angularVelocityRadS();

    sim::Ship b = makeTestShip();
    b.setAngularVelocityRadS({1.0, 2.0, -0.5});
    for (int i = 0; i < 10; ++i) {
      b.step(0.01, in);
    }
    const auto wB = b.angularVelocityRadS();

    CHECK(approxVec(wA, wB, 1e-9));
  }

  // ---- Dampers should decay velocity relative to the damping frame. ----
  {
    sim::Ship s = makeTestShip();
    s.setDampingLinear(0.5);
    s.setDampingFrameVelocityKmS({1.0, 0.0, 0.0});
    s.setVelocityKmS({11.0, 0.0, 0.0});

    const double dt = 0.1;
    s.step(dt, in);

    const double expectedRel = 10.0 * std::exp(-0.5 * dt);
    const math::Vec3d expected{1.0 + expectedRel, 0.0, 0.0};

    CHECK(approxVec(s.velocityKmS(), expected, 1e-9));
  }

  // ---- Brake-only damping should also be dt-invariant. ----
  {
    sim::ShipInput br{};
    br.thrustLocal = {0.0, 0.0, 0.0};
    br.torqueLocal = {0.0, 0.0, 0.0};
    br.dampers = false;
    br.brake = true;
    br.boost = false;

    sim::Ship a = makeTestShip();
    a.setDampingLinear(0.4);
    a.setVelocityKmS({5.0, 1.0, 0.0});
    a.step(0.1, br);
    const auto vA = a.velocityKmS();

    sim::Ship b = makeTestShip();
    b.setDampingLinear(0.4);
    b.setVelocityKmS({5.0, 1.0, 0.0});
    for (int i = 0; i < 10; ++i) {
      b.step(0.01, br);
    }
    const auto vB = b.velocityKmS();

    CHECK(approxVec(vA, vB, 1e-9));
  }

  return failures;
}
