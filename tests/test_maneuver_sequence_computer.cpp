#include "stellar/sim/ManeuverSequenceComputer.h"

#include <cmath>
#include <iostream>

static bool approxEq(double a, double b, double eps) { return std::abs(a - b) <= eps; }

int test_maneuver_sequence_computer() {
  int fails = 0;

  // --- Two-node sequence reaches the sum of delta-v vectors ---
  {
    stellar::sim::Ship ship;
    ship.setPositionKm({0, 0, 0});
    ship.setVelocityKmS({0, 0, 0});
    ship.setOrientation({1, 0, 0, 0});

    stellar::sim::ManeuverSequenceComputer seq;

    double nowDays = 0.0;

    std::vector<stellar::sim::ManeuverPlan> plans;
    plans.push_back({nowDays, {0, 0, 0.2}});                 // 200 m/s prograde
    plans.push_back({nowDays + 20.0 / 86400.0, {0.1, 0, 0}}); // 100 m/s to +X

    seq.engage(ship, plans, /*sortByTime=*/true);

    stellar::sim::ManeuverComputerParams p;
    p.alignToleranceDeg = 35.0; // forgiving to keep test fast
    p.faceGain = 3.2;
    p.dvToleranceKmS = 0.002;
    p.disableDampersDuringBurn = true;
    p.allowBoost = false;
    p.abortAfterMissedSec = 0.0;

    const double dt = 0.1;
    bool sawAdvance = false;
    int steps = 0;

    while (seq.active() && steps < 6000) {
      auto out = seq.update(ship, nowDays, dt, p);
      sawAdvance = sawAdvance || out.advanced;
      ship.step(dt, out.input);
      nowDays += dt / 86400.0;
      ++steps;
    }

    if (seq.phase() != stellar::sim::ManeuverSequencePhase::Complete) {
      std::cerr << "[test_maneuver_sequence_computer] expected Complete, got phase="
                << (int)seq.phase() << "\n";
      ++fails;
    }

    if (!sawAdvance) {
      std::cerr << "[test_maneuver_sequence_computer] expected sequence to advance to node 2\n";
      ++fails;
    }

    const auto v = ship.velocityKmS();
    if (!approxEq(v.x, 0.1, 0.012) || !approxEq(v.z, 0.2, 0.012) || std::abs(v.y) > 0.012) {
      std::cerr << "[test_maneuver_sequence_computer] dv reached (" << v.x << "," << v.y << "," << v.z
                << ") km/s; expected ~ (0.1,0,0.2)\n";
      ++fails;
    }
  }

  // --- Abort propagates from underlying node computer ---
  {
    stellar::sim::Ship ship;
    ship.setVelocityKmS({0, 0, 0});

    stellar::sim::ManeuverSequenceComputer seq;

    const double nowDays = 10.0;

    std::vector<stellar::sim::ManeuverPlan> plans;
    plans.push_back({nowDays - 60.0 / 86400.0, {0, 0, 0.2}}); // missed by 60s
    plans.push_back({nowDays + 10.0 / 86400.0, {0, 0, 0.1}});

    seq.engage(ship, plans, /*sortByTime=*/true);

    stellar::sim::ManeuverComputerParams p;
    p.abortAfterMissedSec = 45.0;

    auto out = seq.update(ship, nowDays, 0.1, p);
    (void)out;

    if (seq.phase() != stellar::sim::ManeuverSequencePhase::Aborted) {
      std::cerr << "[test_maneuver_sequence_computer] expected Aborted when first node missed by 60s\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_maneuver_sequence_computer] pass\n";
  return fails;
}
