#include "stellar/sim/TrajectoryAnalysis.h"

#include "stellar/sim/Units.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

int test_trajectory_analysis() {
  int fails = 0;

  // --- Closest approach to a (nearly) stationary planet ---
  {
    stellar::sim::StarSystem sys;
    sys.stub.id = 1;
    sys.stub.name = "TestSystem";

    stellar::sim::Planet p;
    p.name = "TestPlanet";
    p.massEarth = 1.0;
    p.radiusEarth = 0.01; // ~63.71 km
    p.orbit.eccentricity = 0.0;
    p.orbit.semiMajorAxisAU = 10000.0 / stellar::sim::kAU_KM; // 10,000 km from star
    p.orbit.meanAnomalyAtEpochRad = 0.0;
    p.orbit.epochDays = 0.0;
    p.orbit.periodDays = 1e9; // effectively static for a few seconds
    sys.planets.push_back(p);

    std::vector<stellar::sim::TrajectorySample> samples;
    samples.reserve(11);

    // Ship passes by at a closest distance of 500 km at t=5s.
    for (int i = 0; i <= 10; ++i) {
      stellar::sim::TrajectorySample s;
      s.tSec = (double)i;
      s.posKm = {10500.0, -1000.0 + 200.0 * (double)i, 0.0};
      s.velKmS = {0.0, 200.0, 0.0};
      samples.push_back(s);
    }

    stellar::sim::GravityParams gp;
    gp.includeStar = false;
    gp.includePlanets = true;

    const auto res = stellar::sim::analyzeTrajectory(sys, /*startTimeDays=*/0.0, samples, gp);

    bool found = false;
    for (const auto& a : res.closestApproachByBody) {
      if (!a.valid) continue;
      if (a.body.kind != stellar::sim::GravityBody::Kind::Planet) continue;
      if (a.body.id != 0) continue;
      found = true;

      if (!near(a.tSec, 5.0, 1e-6)) {
        std::cerr << "[test_trajectory_analysis] expected t=5s, got t=" << a.tSec << "\n";
        ++fails;
      }
      if (!near(a.distanceKm, 500.0, 1e-6)) {
        std::cerr << "[test_trajectory_analysis] expected d=500km, got d=" << a.distanceKm << "\n";
        ++fails;
      }
      if (!(a.altitudeKm > 0.0)) {
        std::cerr << "[test_trajectory_analysis] expected positive altitude, got alt=" << a.altitudeKm << "\n";
        ++fails;
      }
      break;
    }

    if (!found) {
      std::cerr << "[test_trajectory_analysis] failed to find planet approach\n";
      ++fails;
    }
  }

  // --- Impact detection: intersecting the planet radius ---
  {
    stellar::sim::StarSystem sys;
    sys.stub.id = 2;
    sys.stub.name = "ImpactSystem";

    stellar::sim::Planet p;
    p.name = "ImpactPlanet";
    p.massEarth = 1.0;
    p.radiusEarth = 1000.0 / 6371.0; // radius = 1000 km
    p.orbit.eccentricity = 0.0;
    p.orbit.semiMajorAxisAU = 10000.0 / stellar::sim::kAU_KM;
    p.orbit.meanAnomalyAtEpochRad = 0.0;
    p.orbit.epochDays = 0.0;
    p.orbit.periodDays = 1e9;
    sys.planets.push_back(p);

    std::vector<stellar::sim::TrajectorySample> samples;
    samples.reserve(11);

    // Pass through the planet center at t=5s (x matches planet, y crosses 0).
    // Entry at y=-1000 => t=2.5, exit at y=+1000 => t=7.5.
    for (int i = 0; i <= 10; ++i) {
      stellar::sim::TrajectorySample s;
      s.tSec = (double)i;
      s.posKm = {10000.0, -2000.0 + 400.0 * (double)i, 0.0};
      s.velKmS = {0.0, 400.0, 0.0};
      samples.push_back(s);
    }

    stellar::sim::GravityParams gp;
    gp.includeStar = false;
    gp.includePlanets = true;

    const auto res = stellar::sim::analyzeTrajectory(sys, /*startTimeDays=*/0.0, samples, gp);
    if (!res.firstImpact.valid) {
      std::cerr << "[test_trajectory_analysis] expected an impact\n";
      ++fails;
    } else {
      if (res.firstImpact.body.kind != stellar::sim::GravityBody::Kind::Planet || res.firstImpact.body.id != 0) {
        std::cerr << "[test_trajectory_analysis] expected impact with planet 0\n";
        ++fails;
      }

      if (!near(res.firstImpact.tSec, 2.5, 1e-6)) {
        std::cerr << "[test_trajectory_analysis] expected impact at t=2.5s, got t=" << res.firstImpact.tSec << "\n";
        ++fails;
      }
    }
  }

  if (fails == 0) std::cout << "[test_trajectory_analysis] pass\n";
  return fails;
}
