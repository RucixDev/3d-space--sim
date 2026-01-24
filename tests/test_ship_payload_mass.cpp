#include "stellar/sim/Ship.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_ship_payload_mass() {
  int failures = 0;

  // No payload => full authority.
  {
    sim::Ship s;
    s.setMassKg(1000.0);
    s.setMaxLinearAccelKmS2(1.0);
    s.setMaxAngularAccelRadS2(2.0);

    s.setPayloadCapacityKg(100.0);
    s.setPayloadHandlingMinScale(0.5);
    s.setPayloadMassKg(0.0);

    sim::ShipInput in{};
    in.thrustLocal = {0, 0, 1};
    in.torqueLocal = {1, 0, 0};
    in.dampers = false;

    s.step(1.0, in);

    CHECK(approx(s.velocityKmS().z, 1.0, 1e-9));
    CHECK(approx(s.angularVelocityRadS().x, 2.0, 1e-9));
  }

  // Full payload => scaled to minScale at capacity.
  {
    sim::Ship s;
    s.setMassKg(1000.0);
    s.setMaxLinearAccelKmS2(1.0);
    s.setMaxAngularAccelRadS2(2.0);

    s.setPayloadCapacityKg(100.0);
    s.setPayloadHandlingMinScale(0.5);
    s.setPayloadMassKg(100.0); // 100% full

    sim::ShipInput in{};
    in.thrustLocal = {0, 0, 1};
    in.torqueLocal = {1, 0, 0};
    in.dampers = false;

    s.step(1.0, in);

    CHECK(approx(s.payloadHandlingScale(), 0.5, 1e-12));
    CHECK(approx(s.velocityKmS().z, 0.5, 1e-9));
    CHECK(approx(s.angularVelocityRadS().x, 1.0, 1e-9));
  }

  // Overloaded payload => soft degradation below minScale but never to zero.
  {
    sim::Ship s;
    s.setMassKg(1000.0);
    s.setMaxLinearAccelKmS2(1.0);

    s.setPayloadCapacityKg(100.0);
    s.setPayloadHandlingMinScale(0.5);
    s.setPayloadMassKg(200.0); // 200% full

    sim::ShipInput in{};
    in.thrustLocal = {0, 0, 1};
    in.dampers = false;

    s.step(1.0, in);

    // Expected: minScale / (1 + 0.5*(frac-1)) where frac=2 => 0.5 / 1.5 = 1/3.
    CHECK(std::fabs(s.velocityKmS().z - (1.0 / 3.0)) < 1e-6);
    CHECK(s.payloadHandlingScale() < 0.5 + 1e-12);
    CHECK(s.payloadHandlingScale() > 0.05 - 1e-12);
  }

  return failures;
}
