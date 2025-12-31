#include "stellar/sim/FlightController.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool near(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

int test_intercept_course() {
  int fails = 0;

  // Scenario: a target is moving laterally. Pure pursuit will "aim at" the target's
  // current position, while intercept-course guidance should lead it.
  {
    sim::Ship s;
    s.setPositionKm({0, 0, 0});
    s.setVelocityKmS({0, 0, 0});
    s.setOrientation({1, 0, 0, 0});
    s.setAngularVelocityRadS({0, 0, 0});

    const math::Vec3d targetPos{0, 0, 20.0};
    const math::Vec3d targetVel{0.30, 0.0, 0.0};

    sim::FlightControlParams fp{};
    fp.maxSpeedKmS = 0.50;   // desired closing speed
    fp.speedGain = 1.0;      // saturate quickly
    fp.velGain = 2.0;
    fp.desiredDistKm = 0.0;
    fp.allowBoost = false;

    sim::AttitudeControlParams ap{};
    ap.faceGain = 2.0;

    // Expected intercept lead direction for this setup:
    // Solve |(v*t, 0, d)| = s*t  => d^2 = (s^2 - v^2) t^2
    // d=20km, v=0.30km/s, s=0.50km/s => t=50s
    // lead vector = (v*t, 0, d) = (15,0,20) -> normalized (0.6,0,0.8)
    const math::Vec3d expectedLead{0.6, 0.0, 0.8};

    sim::InterceptCourseParams ic{};
    ic.enabled = true;
    ic.maxLeadTimeSec = 120.0;
    ic.minSpeedKmS = 0.05;
    ic.useMaxSpeedForSolve = true;

    const auto out = sim::chaseTargetIntercept(s, targetPos, targetVel, fp, ap, ic);

    const math::Vec3d rel = out.desiredVelKmS - targetVel;
    math::Vec3d relN = rel;
    if (relN.lengthSq() > 1e-18) relN = relN.normalized();

    if (!near(relN.x, expectedLead.x, 0.03) ||
        !near(relN.y, expectedLead.y, 0.03) ||
        !near(relN.z, expectedLead.z, 0.03)) {
      std::cerr << "[test_intercept_course] lead direction mismatch. got=("
                << relN.x << "," << relN.y << "," << relN.z << ") expected=("
                << expectedLead.x << "," << expectedLead.y << "," << expectedLead.z << ")\n";
      ++fails;
    }
  }

  // Gating: if maxLeadTimeSec is too small, the lead solve should be rejected and
  // the controller should fall back to pure pursuit (aim at current position).
  {
    sim::Ship s;
    s.setPositionKm({0, 0, 0});
    s.setVelocityKmS({0, 0, 0});
    s.setOrientation({1, 0, 0, 0});

    const math::Vec3d targetPos{0, 0, 20.0};
    const math::Vec3d targetVel{0.30, 0.0, 0.0};

    sim::FlightControlParams fp{};
    fp.maxSpeedKmS = 0.50;
    fp.speedGain = 1.0;
    fp.velGain = 2.0;
    fp.desiredDistKm = 0.0;
    fp.allowBoost = false;

    sim::AttitudeControlParams ap{};
    ap.faceGain = 2.0;

    sim::InterceptCourseParams ic{};
    ic.enabled = true;
    ic.maxLeadTimeSec = 10.0; // solution requires ~50s, so this should reject
    ic.minSpeedKmS = 0.05;
    ic.useMaxSpeedForSolve = true;

    const auto out = sim::chaseTargetIntercept(s, targetPos, targetVel, fp, ap, ic);
    const math::Vec3d rel = out.desiredVelKmS - targetVel;
    math::Vec3d relN = rel;
    if (relN.lengthSq() > 1e-18) relN = relN.normalized();

    // Pure pursuit aims straight down +Z here.
    if (std::abs(relN.x) > 0.08 || std::abs(relN.y) > 0.08 || relN.z < 0.90) {
      std::cerr << "[test_intercept_course] expected fallback to pure pursuit. got=("
                << relN.x << "," << relN.y << "," << relN.z << ")\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_intercept_course] PASS\n";
  }
  return fails;
}
