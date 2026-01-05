#include "stellar/sim/ManeuverComputer.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) { return std::abs(a - b) <= eps; }

int test_maneuver_computer() {
  int fails = 0;

  // --- Basic immediate burn reaches requested delta-v (continuous approximation) ---
  {
    stellar::sim::Ship ship;
    ship.setPositionKm({0,0,0});
    ship.setVelocityKmS({0,0,0});
    ship.setOrientation({1,0,0,0});

    stellar::sim::ManeuverComputer mc;
    stellar::sim::ManeuverPlan plan;
    double nowDays = 0.0;

    plan.nodeTimeDays = nowDays; // start now
    plan.deltaVWorldKmS = {0, 0, 0.5}; // 500 m/s prograde along +Z

    mc.engage(ship, plan);

    stellar::sim::ManeuverComputerParams p;
    p.alignToleranceDeg = 25.0;         // forgiving for the unit test
    p.faceGain = 3.0;
    p.dvToleranceKmS = 0.002;           // 2 m/s
    p.disableDampersDuringBurn = true;
    p.allowBoost = false;

    const double dt = 0.1;
    int steps = 0;
    while (mc.active() && steps < 2000) {
      auto out = mc.update(ship, nowDays, dt, p);
      ship.step(dt, out.input);
      nowDays += dt / 86400.0;
      ++steps;
    }

    if (mc.phase() != stellar::sim::ManeuverComputerPhase::Complete) {
      std::cerr << "[test_maneuver_computer] expected Complete, got phase=" << (int)mc.phase() << "\n";
      ++fails;
    }

    const auto v = ship.velocityKmS();
    if (!near(v.z, 0.5, 0.006) || std::abs(v.x) > 0.006 || std::abs(v.y) > 0.006) {
      std::cerr << "[test_maneuver_computer] dv reached (" << v.x << "," << v.y << "," << v.z
                << ") km/s; expected ~ (0,0,0.5)\n";
      ++fails;
    }
  }

  // --- Burn start timing: should start ~half burn duration early ---
  {
    stellar::sim::Ship ship;
    ship.setPositionKm({0,0,0});
    ship.setVelocityKmS({0,0,0});
    ship.setOrientation({1,0,0,0});

    stellar::sim::ManeuverComputer mc;
    stellar::sim::ManeuverPlan plan;
    double nowDays = 0.0;

    plan.nodeTimeDays = nowDays + 20.0 / 86400.0; // 20s from now
    plan.deltaVWorldKmS = {0, 0, 0.5};            // 10s burn at 0.05 km/s^2
    mc.engage(ship, plan);

    stellar::sim::ManeuverComputerParams p;
    p.alignToleranceDeg = 45.0;
    p.faceGain = 2.5;
    p.dvToleranceKmS = 0.002;
    p.disableDampersDuringBurn = true;
    p.allowBoost = false;
    p.extraLeadTimeSec = 0.0;
    p.abortAfterMissedSec = 0.0;

    const double dt = 1.0;
    bool startedAt15 = false;

    for (int t = 0; t <= 16; ++t) {
      auto out = mc.update(ship, nowDays, dt, p);

      if (t == 14) {
        if (out.phase != stellar::sim::ManeuverComputerPhase::Orient || std::abs(out.input.thrustLocal.z) > 1e-6) {
          std::cerr << "[test_maneuver_computer] expected no burn at t=14s\n";
          ++fails;
          break;
        }
      }

      if (t == 15) {
        startedAt15 = (out.phase == stellar::sim::ManeuverComputerPhase::Burn && out.input.thrustLocal.z > 1e-6);
      }

      ship.step(dt, out.input);
      nowDays += dt / 86400.0;
    }

    if (!startedAt15) {
      std::cerr << "[test_maneuver_computer] expected burn to start at ~15s\n";
      ++fails;
    }
  }

  // --- Abort behavior when node is badly missed ---
  {
    stellar::sim::Ship ship;
    ship.setVelocityKmS({0,0,0});

    stellar::sim::ManeuverComputer mc;
    stellar::sim::ManeuverPlan plan;
    const double nowDays = 10.0;

    plan.nodeTimeDays = nowDays - 60.0 / 86400.0; // missed by 60s
    plan.deltaVWorldKmS = {0, 0, 0.2};

    mc.engage(ship, plan);

    stellar::sim::ManeuverComputerParams p;
    p.abortAfterMissedSec = 45.0;

    auto out = mc.update(ship, nowDays, 0.1, p);
    if (out.phase != stellar::sim::ManeuverComputerPhase::Aborted) {
      std::cerr << "[test_maneuver_computer] expected Aborted when missed by 60s\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_maneuver_computer] pass\n";
  return fails;
}
