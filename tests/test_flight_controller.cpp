#include "stellar/sim/FlightController.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool inRange(double v, double lo, double hi, double eps = 1e-12) {
  return v >= (lo - eps) && v <= (hi + eps);
}

int test_flight_controller() {
  int fails = 0;

  // --- Basic chase: converge toward a point and roughly face it. ---
  {
    sim::Ship s;
    s.setPositionKm({0, 0, 0});
    s.setVelocityKmS({0, 0, 0});
    s.setOrientation({1, 0, 0, 0});
    s.setAngularVelocityRadS({0, 0, 0});

    const math::Vec3d tgt{0, 0, 5.0};

    sim::FlightControlParams fp{};
    fp.maxSpeedKmS = 0.60;
    fp.speedGain = 0.25;   // small-distance test: reach maxSpeed quickly
    fp.velGain = 2.0;
    fp.desiredDistKm = 0.0;
    fp.allowBoost = false;
    fp.dampers = true;

    sim::AttitudeControlParams ap{};
    ap.faceGain = 2.0;
    ap.alignUp = false;

    const double dt = 0.10;
    for (int i = 0; i < 260; ++i) {
      const auto out = sim::chaseTarget(s, tgt, {0, 0, 0}, fp, ap);
      s.step(dt, out.input);
    }

    const double dist = (tgt - s.positionKm()).length();
    if (dist > 0.35) {
      std::cerr << "[test_flight_controller] chase did not converge. dist=" << dist << "\n";
      ++fails;
    }

    const math::Vec3d to = (tgt - s.positionKm());
    if (to.lengthSq() > 1e-12) {
      const double aim = math::dot(s.forward().normalized(), to.normalized());
      if (aim < 0.70) {
        std::cerr << "[test_flight_controller] chase did not face target. aim=" << aim << "\n";
        ++fails;
      }
    }
  }

  // --- Standoff: try to maintain ~desiredDistKm. ---
  {
    sim::Ship s;
    s.setPositionKm({0, 0, 0});
    s.setVelocityKmS({0, 0, 0});
    s.setOrientation({1, 0, 0, 0});
    s.setAngularVelocityRadS({0, 0, 0});

    const math::Vec3d tgt{0, 0, 10.0};

    sim::FlightControlParams fp{};
    fp.maxSpeedKmS = 0.75;
    fp.speedGain = 0.20;
    fp.velGain = 2.0;
    fp.desiredDistKm = 2.0;
    fp.allowBoost = false;
    fp.dampers = true;

    sim::AttitudeControlParams ap{};
    ap.faceGain = 1.7;

    const double dt = 0.10;
    for (int i = 0; i < 380; ++i) {
      const auto out = sim::chaseTarget(s, tgt, {0, 0, 0}, fp, ap);
      s.step(dt, out.input);
    }

    const double dist = (tgt - s.positionKm()).length();
    if (!inRange(dist, 0.9, 3.3)) {
      std::cerr << "[test_flight_controller] standoff out of range. dist=" << dist << "\n";
      ++fails;
    }
  }

  // --- Output safety: inputs should stay in normalized ranges (even under large dv). ---
  {
    sim::Ship s;
    s.setPositionKm({0, 0, 0});
    s.setVelocityKmS({8.0, -3.0, 2.0});

    sim::FlightControlParams fp{};
    fp.maxSpeedKmS = 1.0;
    fp.speedGain = 0.10;
    fp.velGain = 8.0;
    fp.desiredDistKm = 0.0;
    fp.allowBoost = false;

    sim::AttitudeControlParams ap{};
    ap.faceGain = 2.5;

    const auto out = sim::chaseTarget(s, {0, 0, 50.0}, {0, 0, 0}, fp, ap);

    if (!inRange(out.input.thrustLocal.x, -1.0, 1.0) ||
        !inRange(out.input.thrustLocal.y, -1.0, 1.0) ||
        !inRange(out.input.thrustLocal.z, -1.0, 1.0)) {
      std::cerr << "[test_flight_controller] thrustLocal out of range.\n";
      ++fails;
    }
    if (!inRange(out.input.torqueLocal.x, -1.0, 1.0) ||
        !inRange(out.input.torqueLocal.y, -1.0, 1.0) ||
        !inRange(out.input.torqueLocal.z, -1.0, 1.0)) {
      std::cerr << "[test_flight_controller] torqueLocal out of range.\n";
      ++fails;
    }
    if (out.input.boost) {
      std::cerr << "[test_flight_controller] unexpected boost when allowBoost=false.\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_flight_controller] PASS\n";
  }
  return fails;
}
