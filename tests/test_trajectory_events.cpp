#include "stellar/sim/TrajectoryEvents.h"

#include "stellar/sim/Celestial.h"
#include "stellar/sim/Units.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

int test_trajectory_events() {
  int fails = 0;

  // --- Dominant gravity transitions: star -> planet along a straight pass ---
  {
    stellar::sim::StarSystem sys;
    sys.stub.id = 42;
    sys.stub.name = "DomChangeSys";

    sys.star.massSol = 1e-6;
    sys.star.radiusSol = 1e-6;

    stellar::sim::Planet p;
    p.name = "PlanetX";
    p.massEarth = 1.0;
    p.radiusEarth = 0.01;
    p.orbit.eccentricity = 0.0;
    const double D = 1.0e6; // km
    p.orbit.semiMajorAxisAU = D / stellar::sim::kAU_KM;
    p.orbit.meanAnomalyAtEpochRad = 0.0;
    p.orbit.epochDays = 0.0;
    p.orbit.periodDays = 1e9; // effectively static
    sys.planets.push_back(p);

    const double muStar = stellar::sim::muStarKm3S2(sys.star);
    const double muPlanet = stellar::sim::muPlanetKm3S2(sys.planets[0]);
    if (!(muStar > 0.0 && muPlanet > 0.0)) {
      std::cerr << "[test_trajectory_events] invalid mu\n";
      ++fails;
    }

    // For two point masses on a line: mu_s/x^2 = mu_p/(D-x)^2
    // => x = D / (1 + sqrt(mu_p/mu_s))
    const double xTransition = D / (1.0 + std::sqrt(muPlanet / muStar));

    const double x0 = 0.1 * D;
    const double x1 = 0.9 * D;
    const double T = 1000.0;
    const double v = (x1 - x0) / T;

    std::vector<stellar::sim::TrajectorySample> samples;
    for (int i = 0; i <= 10; ++i) {
      const double t = 100.0 * (double)i;
      const double x = x0 + v * t;

      stellar::sim::TrajectorySample s;
      s.tSec = t;
      s.posKm = {x, 0.0, 0.0};
      s.velKmS = {v, 0.0, 0.0};
      samples.push_back(s);
    }

    stellar::sim::GravityParams gp;
    gp.includeStar = true;
    gp.includePlanets = true;

    stellar::sim::DominantBodyTransitionParams dp;
    dp.refineDepth = 18;

    const auto events = stellar::sim::detectDominantBodyTransitions(sys, /*startTimeDays=*/0.0, samples, gp, dp);

    if (events.empty()) {
      std::cerr << "[test_trajectory_events] expected >=1 transition\n";
      ++fails;
    } else {
      const auto& e = events.front();

      if (e.from.kind != stellar::sim::GravityBody::Kind::Star) {
        std::cerr << "[test_trajectory_events] expected from=Star\n";
        ++fails;
      }
      if (e.to.kind != stellar::sim::GravityBody::Kind::Planet || e.to.id != 0) {
        std::cerr << "[test_trajectory_events] expected to=Planet[0]\n";
        ++fails;
      }

      const double expectedTSec = (xTransition - x0) / v;
      if (!near(e.tSec, expectedTSec, 1.0)) {
        std::cerr << "[test_trajectory_events] expected t~" << expectedTSec << "s, got t=" << e.tSec << "s\n";
        ++fails;
      }
    }
  }

  // --- No included bodies => no transitions ---
  {
    stellar::sim::StarSystem sys;
    sys.stub.id = 7;
    sys.stub.name = "NoBodies";

    std::vector<stellar::sim::TrajectorySample> samples;
    {
      stellar::sim::TrajectorySample a;
      a.tSec = 0.0;
      a.posKm = {0.0, 0.0, 0.0};
      a.velKmS = {0.0, 0.0, 0.0};
      samples.push_back(a);
    }
    {
      stellar::sim::TrajectorySample b;
      b.tSec = 1.0;
      b.posKm = {1.0, 0.0, 0.0};
      b.velKmS = {1.0, 0.0, 0.0};
      samples.push_back(b);
    }

    stellar::sim::GravityParams gp;
    gp.includeStar = false;
    gp.includePlanets = false;

    const auto events = stellar::sim::detectDominantBodyTransitions(sys, /*startTimeDays=*/0.0, samples, gp);
    if (!events.empty()) {
      std::cerr << "[test_trajectory_events] expected no events when no bodies are included\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_trajectory_events] pass\n";
  return fails;
}
