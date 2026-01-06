#include "stellar/sim/TrajectoryPredictor.h"

#include "stellar/math/Math.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps = 1e-6) { return std::abs(a - b) <= eps; }

int test_trajectory_predictor() {
  int fails = 0;

  // --- No gravity: straight line ---
  {
    stellar::sim::StarSystem sys;
    sys.stub.id = 1;

    stellar::sim::TrajectoryPredictParams p;
    p.horizonSec = 10.0;
    p.stepSec = 1.0;
    p.includeGravity = false;

    const stellar::math::Vec3d x0{0, 0, 0};
    const stellar::math::Vec3d v0{1.0, 2.0, -0.5};

    const auto samples = stellar::sim::predictTrajectoryRK4(sys, /*startTimeDays=*/0.0, x0, v0, p);
    if (samples.size() < 2) {
      std::cerr << "[test_trajectory_predictor] expected >= 2 samples\n";
      ++fails;
    } else {
      for (const auto& s : samples) {
        const stellar::math::Vec3d expected = x0 + v0 * s.tSec;
        const double err = (s.posKm - expected).length();
        if (err > 1e-9) {
          std::cerr << "[test_trajectory_predictor] straight-line position error=" << err
                    << " at t=" << s.tSec << "\n";
          ++fails;
          break;
        }
      }
    }
  }

  // --- Maneuver node: instantaneous delta-v applied at exact time ---
  {
    stellar::sim::StarSystem sys;
    sys.stub.id = 2;

    stellar::sim::TrajectoryPredictParams p;
    p.horizonSec = 12.0;
    p.stepSec = 2.0;
    p.includeGravity = false;

    stellar::sim::ManeuverNode node;
    node.timeSec = 5.0;
    node.deltaVKmS = {0.25, -0.5, 0.0};

    const stellar::math::Vec3d x0{0, 0, 0};
    const stellar::math::Vec3d v0{1.0, 0.0, 0.0};

    const auto samples = stellar::sim::predictTrajectoryRK4(sys, 0.0, x0, v0, p, &node);
    bool found = false;
    for (const auto& s : samples) {
      if (std::abs(s.tSec - node.timeSec) < 1e-6) {
        found = true;
        const stellar::math::Vec3d expectedVel = v0 + node.deltaVKmS;
        if ((s.velKmS - expectedVel).length() > 1e-9) {
          std::cerr << "[test_trajectory_predictor] node vel mismatch at t=" << s.tSec << "\n";
          ++fails;
        }
        break;
      }
    }
    if (!found) {
      std::cerr << "[test_trajectory_predictor] expected an exact node sample at t=" << node.timeSec << "\n";
      ++fails;
    }
  }

  // --- Gravity: near-circular orbit stays near constant radius (RK4 sanity) ---
  {
    stellar::sim::StarSystem sys;
    sys.stub.id = 3;

    // Configure the system "star" to have Earth's mass & radius.
    // This keeps the test self-contained while exercising systemGravityAccelKmS2().
    const double earthMassSol = stellar::sim::kMassEarthKg / stellar::sim::kMassSunKg;
    const double earthRadiusSol = stellar::sim::kRadiusEarthKm / stellar::sim::kRadiusSunKm;
    sys.star.massSol = earthMassSol;
    sys.star.radiusSol = earthRadiusSol;

    stellar::sim::TrajectoryPredictParams p;
    p.horizonSec = 180.0; // 3 minutes
    p.stepSec = 1.0;
    p.includeGravity = true;
    p.gravity.includeStar = true;
    p.gravity.includePlanets = false;
    p.gravity.scale = 1.0;
    p.gravity.softeningRadiusScale = 1.0;

    const double mu = stellar::sim::muFromEarthMass(1.0);
    const double rKm = 7000.0;
    const double vKmS = std::sqrt(mu / rKm);

    const stellar::math::Vec3d x0{rKm, 0, 0};
    const stellar::math::Vec3d v0{0, vKmS, 0};

    const auto samples = stellar::sim::predictTrajectoryRK4(sys, 0.0, x0, v0, p);
    double maxDev = 0.0;
    for (const auto& s : samples) {
      const double r = s.posKm.length();
      maxDev = std::max(maxDev, std::abs(r - rKm));
    }

    // RK4 with dt=1s should keep deviation small over a few minutes.
    if (maxDev > 5.0) {
      std::cerr << "[test_trajectory_predictor] circular radius drift too high: maxDev=" << maxDev << " km\n";
      ++fails;
    }
  }

  // --- Adaptive RK45: no gravity stays a perfect straight line (within float noise) ---
  {
    stellar::sim::StarSystem sys;
    sys.stub.id = 4;

    stellar::sim::TrajectoryPredictParams p;
    p.horizonSec = 10.0;
    p.stepSec = 1.0;        // output samples each second
    p.includeGravity = false;

    // Tight tolerances (should accept large steps because the solution is exactly linear).
    p.minStepSec = 1e-3;
    p.maxStepSec = 1.0;
    p.relTol = 1e-10;
    p.absTolPosKm = 1e-12;
    p.absTolVelKmS = 1e-12;

    const stellar::math::Vec3d x0{0, 0, 0};
    const stellar::math::Vec3d v0{1.0, 2.0, -0.5};

    const auto samples = stellar::sim::predictTrajectoryRK45Adaptive(sys, 0.0, x0, v0, p);

    // Expect samples at t=0..10 inclusive.
    if (samples.size() != 11) {
      std::cerr << "[test_trajectory_predictor] RK45 expected 11 samples, got " << samples.size() << "\n";
      ++fails;
    }

    for (const auto& s : samples) {
      const stellar::math::Vec3d expected = x0 + v0 * s.tSec;
      const double err = (s.posKm - expected).length();
      if (err > 1e-8) {
        std::cerr << "[test_trajectory_predictor] RK45 straight-line position error=" << err
                  << " at t=" << s.tSec << "\n";
        ++fails;
        break;
      }
    }
  }

  // --- Adaptive RK45: maneuver node applied at exact time ---
  {
    stellar::sim::StarSystem sys;
    sys.stub.id = 5;

    stellar::sim::TrajectoryPredictParams p;
    p.horizonSec = 12.0;
    p.stepSec = 2.0;        // output samples at 0,2,4,6,... plus node sample at 5
    p.includeGravity = false;

    p.minStepSec = 1e-3;
    p.maxStepSec = 2.0;
    p.relTol = 1e-10;
    p.absTolPosKm = 1e-12;
    p.absTolVelKmS = 1e-12;

    stellar::sim::ManeuverNode node;
    node.timeSec = 5.0;
    node.deltaVKmS = {0.25, -0.5, 0.0};

    const stellar::math::Vec3d x0{0, 0, 0};
    const stellar::math::Vec3d v0{1.0, 0.0, 0.0};

    const auto samples = stellar::sim::predictTrajectoryRK45Adaptive(sys, 0.0, x0, v0, p, &node);

    bool found = false;
    for (const auto& s : samples) {
      if (std::abs(s.tSec - node.timeSec) < 1e-6) {
        found = true;
        const stellar::math::Vec3d expectedVel = v0 + node.deltaVKmS;
        if ((s.velKmS - expectedVel).length() > 1e-9) {
          std::cerr << "[test_trajectory_predictor] RK45 node vel mismatch at t=" << s.tSec << "\n";
          ++fails;
        }
        break;
      }
    }

    if (!found) {
      std::cerr << "[test_trajectory_predictor] RK45 expected an exact node sample at t=" << node.timeSec << "\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_trajectory_predictor] pass\n";
  return fails;
}
