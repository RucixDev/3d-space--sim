#include "stellar/sim/ManeuverProgramComputer.h"

#include "stellar/sim/Gravity.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

int test_maneuver_program_computer() {
  int fails = 0;

  const double mu = stellar::sim::muFromEarthMass(1.0);
  const double earthRadius = stellar::sim::radiusKmFromEarthRadius(1.0);

  auto makeEarth = [&](double muVal) {
    stellar::sim::GravityBody earth;
    earth.kind = stellar::sim::GravityBody::Kind::Planet;
    earth.id = 0;
    earth.name = "Earth";
    earth.posKm = {0, 0, 0};
    earth.velKmS = {0, 0, 0};
    earth.muKm3S2 = muVal;
    earth.radiusKm = earthRadius;
    return earth;
  };

  // --- Circularize at apoapsis from periapsis state (Earth-like mu) ---
  {
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

    const auto earth = makeEarth(mu);

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
    const double r = 7000.0;
    const double vCirc = std::sqrt(mu / r);

    stellar::sim::Ship ship;
    ship.setPositionKm({r, 0, 0});
    ship.setVelocityKmS({0, vCirc, 0});

    const auto earth = makeEarth(mu);

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

  // --- Set apoapsis (burn @ periapsis) ---
  {
    const double rp = 7000.0;
    const double ra = 12000.0;
    const double a = (rp + ra) * 0.5;

    const double vPeri = std::sqrt(mu * (2.0 / rp - 1.0 / a));

    const double raTarget = 15000.0;
    const double aNew = 0.5 * (rp + raTarget);
    const double vPeriNew = std::sqrt(mu * (2.0 / rp - 1.0 / aNew));
    const double dvExpected = vPeriNew - vPeri;

    stellar::sim::Ship ship;
    ship.setPositionKm({rp, 0, 0});
    ship.setVelocityKmS({0, vPeri, 0});

    const auto earth = makeEarth(mu);
    const double nowDays = 0.0;

    const auto res = stellar::sim::planSetApoapsisAtPeriapsis(ship, nowDays, earth, raTarget);

    if (!res.valid) {
      std::cerr << "[test_maneuver_program_computer] expected valid set-apo plan, got: " << res.reason << "\n";
      ++fails;
    } else {
      if (!near(res.timeToNodeSec, 0.0, 1e-9)) {
        std::cerr << "[test_maneuver_program_computer] expected periapsis node now, got t=" << res.timeToNodeSec << "\n";
        ++fails;
      }

      if (!near(res.targetRadiusKm, raTarget, 1e-6)) {
        std::cerr << "[test_maneuver_program_computer] target apo radius mismatch: got=" << res.targetRadiusKm
                  << " expected=" << raTarget << "\n";
        ++fails;
      }

      if (!near(res.dvKmS, std::abs(dvExpected), 1e-3)) {
        std::cerr << "[test_maneuver_program_computer] set-apo dv mismatch: got=" << res.dvKmS
                  << " expected=" << std::abs(dvExpected) << "\n";
        ++fails;
      }

      // Burn should be prograde at periapsis (+Y).
      if (res.plan.deltaVWorldKmS.y <= 0.0) {
        std::cerr << "[test_maneuver_program_computer] expected +Y dv at periapsis\n";
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

  // --- Set periapsis (burn @ apoapsis) ---
  {
    const double rp = 7000.0;
    const double ra = 12000.0;
    const double a = (rp + ra) * 0.5;

    const double vApo = std::sqrt(mu * (2.0 / ra - 1.0 / a));

    const double rpTarget = 8000.0;
    const double aNew = 0.5 * (ra + rpTarget);
    const double vApoNew = std::sqrt(mu * (2.0 / ra - 1.0 / aNew));
    const double dvExpected = vApoNew - vApo;

    stellar::sim::Ship ship;
    // Apoapsis state for the same orbit: x=-ra, v along -Y.
    ship.setPositionKm({-ra, 0, 0});
    ship.setVelocityKmS({0, -vApo, 0});

    const auto earth = makeEarth(mu);
    const double nowDays = 0.0;

    const auto res = stellar::sim::planSetPeriapsisAtApoapsis(ship, nowDays, earth, rpTarget);

    if (!res.valid) {
      std::cerr << "[test_maneuver_program_computer] expected valid set-peri plan, got: " << res.reason << "\n";
      ++fails;
    } else {
      if (!near(res.timeToNodeSec, 0.0, 1e-9)) {
        std::cerr << "[test_maneuver_program_computer] expected apoapsis node now, got t=" << res.timeToNodeSec << "\n";
        ++fails;
      }

      if (!near(res.targetRadiusKm, rpTarget, 1e-6)) {
        std::cerr << "[test_maneuver_program_computer] target peri radius mismatch: got=" << res.targetRadiusKm
                  << " expected=" << rpTarget << "\n";
        ++fails;
      }

      if (!near(res.dvKmS, std::abs(dvExpected), 1e-3)) {
        std::cerr << "[test_maneuver_program_computer] set-peri dv mismatch: got=" << res.dvKmS
                  << " expected=" << std::abs(dvExpected) << "\n";
        ++fails;
      }

      // Burn should be prograde at apoapsis (velocity is -Y, so dv should be -Y).
      if (res.plan.deltaVWorldKmS.y >= 0.0) {
        std::cerr << "[test_maneuver_program_computer] expected -Y dv at apoapsis\n";
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

  // --- Escape now: prograde burn to v_escape at current radius ---
  {
    const double r = 7000.0;
    const double vCirc = std::sqrt(mu / r);
    const double vEsc = std::sqrt(2.0 * mu / r);
    const double dvExpected = vEsc - vCirc;

    stellar::sim::Ship ship;
    ship.setPositionKm({r, 0, 0});
    ship.setVelocityKmS({0, vCirc, 0});

    const auto earth = makeEarth(mu);
    const double nowDays = 123.0;

    const auto res = stellar::sim::planEscapeNow(ship, nowDays, earth);

    if (!res.valid) {
      std::cerr << "[test_maneuver_program_computer] expected valid escape plan, got: " << res.reason << "\n";
      ++fails;
    } else {
      if (!near(res.timeToNodeSec, 0.0, 1e-12)) {
        std::cerr << "[test_maneuver_program_computer] expected escape node now, got t=" << res.timeToNodeSec << "\n";
        ++fails;
      }
      if (!near(res.dvKmS, dvExpected, 1e-3)) {
        std::cerr << "[test_maneuver_program_computer] escape dv mismatch: got=" << res.dvKmS
                  << " expected=" << dvExpected << "\n";
        ++fails;
      }
      if (res.plan.deltaVWorldKmS.y <= 0.0) {
        std::cerr << "[test_maneuver_program_computer] expected +Y dv for escape (prograde)\n";
        ++fails;
      }
    }
  }

  // --- Plane align: inclination change at ascending node ---
  {
    const double r = 7000.0;
    const double vCirc = std::sqrt(mu / r);
    const double incDeg = 30.0;
    const double incRad = incDeg * 3.14159265358979323846 / 180.0;

    // Ascending node on +X, inclined about +X.
    stellar::sim::Ship ship;
    ship.setPositionKm({r, 0, 0});
    ship.setVelocityKmS({0, vCirc * std::cos(incRad), vCirc * std::sin(incRad)});

    const auto earth = makeEarth(mu);
    const double nowDays = 0.0;

    const auto res = stellar::sim::planAlignPlaneAtAscendingNode(ship, nowDays, earth, /*forcePrograde=*/true);

    if (!res.valid) {
      std::cerr << "[test_maneuver_program_computer] expected valid plane-align plan, got: " << res.reason << "\n";
      ++fails;
    } else {
      if (!near(res.timeToNodeSec, 0.0, 1e-6)) {
        std::cerr << "[test_maneuver_program_computer] expected AN node now, got t=" << res.timeToNodeSec << "\n";
        ++fails;
      }

      // Expected: rotate velocity around +X by -inc to remove +Z component.
      const double dvY = vCirc - vCirc * std::cos(incRad);
      const double dvZ = -vCirc * std::sin(incRad);
      const double dvExpected = std::sqrt(dvY * dvY + dvZ * dvZ);

      if (!near(res.dvKmS, dvExpected, 1e-3)) {
        std::cerr << "[test_maneuver_program_computer] plane-align dv mismatch: got=" << res.dvKmS
                  << " expected=" << dvExpected << "\n";
        ++fails;
      }

      if (std::abs(res.plan.deltaVWorldKmS.x) > 1e-3) {
        std::cerr << "[test_maneuver_program_computer] expected dv.x ~ 0, got " << res.plan.deltaVWorldKmS.x << "\n";
        ++fails;
      }
      if (res.plan.deltaVWorldKmS.y <= 0.0 || res.plan.deltaVWorldKmS.z >= 0.0) {
        std::cerr << "[test_maneuver_program_computer] expected +Y and -Z dv components, got ("
                  << res.plan.deltaVWorldKmS.x << "," << res.plan.deltaVWorldKmS.y << "," << res.plan.deltaVWorldKmS.z
                  << ")\n";
        ++fails;
      }

      if (!near(res.plan.deltaVWorldKmS.y, dvY, 1e-3) || !near(res.plan.deltaVWorldKmS.z, dvZ, 1e-3)) {
        std::cerr << "[test_maneuver_program_computer] plane-align dv components mismatch: got ("
                  << res.plan.deltaVWorldKmS.y << "," << res.plan.deltaVWorldKmS.z
                  << ") expected (" << dvY << "," << dvZ << ")\n";
        ++fails;
      }
    }
  }

  if (fails == 0) std::cout << "[test_maneuver_program_computer] pass\n";
  return fails;
}
