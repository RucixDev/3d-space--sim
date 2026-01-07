#include "stellar/sim/ManeuverProgramComputer.h"

#include "stellar/sim/Gravity.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

int test_maneuver_program_computer() {
  int fails = 0;

  // --- Circularize at apoapsis from periapsis state (Earth-like mu) ---
  {
    const double mu = stellar::sim::muFromEarthMass(1.0);

    const double rp = 7000.0;
    const double ra = 12000.0;
    const double a = (rp + ra) * 0.5;

    const double vPeri = std::sqrt(mu * (2.0 / rp - 1.0 / a));
    const double vApo = std::sqrt(mu * (2.0 / ra - 1.0 / a));
    const double vCircApo = std::sqrt(mu / ra);
    const double dvExpected = vCircApo - vApo;

    const double period = 2.0 * 3.14159265358979323846 * std::sqrt((a * a * a) / mu);
    const double tApoExpected = period * 0.5;

    stellar::sim::Ship ship;
    ship.setPositionKm({rp, 0, 0});
    ship.setVelocityKmS({0, vPeri, 0});

    stellar::sim::GravityBody earth;
    earth.kind = stellar::sim::GravityBody::Kind::Planet;
    earth.id = 0;
    earth.name = "Earth";
    earth.posKm = {0, 0, 0};
    earth.velKmS = {0, 0, 0};
    earth.muKm3S2 = mu;
    earth.radiusKm = stellar::sim::radiusKmFromEarthRadius(1.0);

    const double nowDays = 0.0;
    const auto res = stellar::sim::planCircularize(
      ship,
      nowDays,
      earth,
      stellar::sim::ManeuverProgramKind::CircularizeAtApoapsis);

    if (!res.valid) {
      std::cerr << "[test_maneuver_program_computer] expected valid plan, got: " << res.reason << "\n";
      ++fails;
    } else {
      if (!near(res.timeToNodeSec, tApoExpected, 1e-1)) {
        std::cerr << "[test_maneuver_program_computer] timeToNodeSec mismatch: got=" << res.timeToNodeSec
                  << " expected=" << tApoExpected << "\n";
        ++fails;
      }

      if (!near(res.dvKmS, std::abs(dvExpected), 1e-3)) {
        std::cerr << "[test_maneuver_program_computer] dv mismatch: got=" << res.dvKmS
                  << " expected=" << std::abs(dvExpected) << "\n";
        ++fails;
      }

      // At apoapsis the velocity is along -Y for this initial condition.
      if (res.plan.deltaVWorldKmS.y >= 0.0) {
        std::cerr << "[test_maneuver_program_computer] expected prograde dv along -Y at apoapsis\n";
        ++fails;
      }
      if (std::abs(res.plan.deltaVWorldKmS.x) > 1e-3 || std::abs(res.plan.deltaVWorldKmS.z) > 1e-3) {
        std::cerr << "[test_maneuver_program_computer] expected mostly tangential dv, got ("
                  << res.plan.deltaVWorldKmS.x << "," << res.plan.deltaVWorldKmS.y << ","
                  << res.plan.deltaVWorldKmS.z << ")\n";
        ++fails;
      }
    }
  }

  // --- Already circular: should return a ~zero dv plan at now ---
  {
    const double mu = stellar::sim::muFromEarthMass(1.0);
    const double r = 7000.0;
    const double vCirc = std::sqrt(mu / r);

    stellar::sim::Ship ship;
    ship.setPositionKm({r, 0, 0});
    ship.setVelocityKmS({0, vCirc, 0});

    stellar::sim::GravityBody earth;
    earth.kind = stellar::sim::GravityBody::Kind::Planet;
    earth.id = 0;
    earth.name = "Earth";
    earth.posKm = {0, 0, 0};
    earth.velKmS = {0, 0, 0};
    earth.muKm3S2 = mu;

    const double nowDays = 10.0;
    const auto res = stellar::sim::planCircularize(
      ship,
      nowDays,
      earth,
      stellar::sim::ManeuverProgramKind::CircularizeAtPeriapsis);

    if (!res.valid) {
      std::cerr << "[test_maneuver_program_computer] expected valid circular plan, got: " << res.reason << "\n";
      ++fails;
    } else {
      if (res.dvKmS > 1e-6) {
        std::cerr << "[test_maneuver_program_computer] expected ~0 dv, got dv=" << res.dvKmS << "\n";
        ++fails;
      }
      if (!near(res.timeToNodeSec, 0.0, 1e-12)) {
        std::cerr << "[test_maneuver_program_computer] expected immediate node, got t=" << res.timeToNodeSec << "\n";
        ++fails;
      }
      if (!near(res.plan.nodeTimeDays, nowDays, 1e-12)) {
        std::cerr << "[test_maneuver_program_computer] expected nodeTimeDays==now\n";
        ++fails;
      }
    }
  }

  if (fails == 0) std::cout << "[test_maneuver_program_computer] pass\n";
  return fails;
}
