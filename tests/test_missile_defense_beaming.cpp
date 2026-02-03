#include "stellar/sim/MissileDefense.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_missile_defense_beaming() {
  int fails = 0;

  // Scenario: inbound radar missile with a doppler notch.
  //
  // We construct a clean geometry where the classic closest-approach jink produces
  // a purely Y/Z-plane direction (x ~= 0), and enabling the radar-beaming bias
  // introduces a positive X component.
  {
    sim::Missile m{};
    m.seeker = sim::MissileSeekerType::Radar;
    m.radarDopplerNotchKmS = 0.75;

    // Target at origin.
    const math::Vec3d tgtPos{0, 0, 0};

    // Place the missile off-axis so LOS has no X component.
    m.posKm = {0, -100, 100};

    // Keep relative velocity purely along -Z so the closest-approach miss vector
    // is exactly -Y (-> jink direction +Y before LOS projection).
    const math::Vec3d tgtVel{10, 5, 0};
    m.velKmS = {10, 5, -30};

    const core::u64 seed = 0x1234abcdull;

    sim::MissileEvasionParams pJink{};
    pJink.enforceLateralToLos = true;
    pJink.enableRadarBeaming = false;

    const sim::MissileEvasionPlan jink = sim::planMissileEvasion(m, tgtPos, tgtVel, seed, pJink);
    if (!jink.valid) {
      std::cerr << "[test_missile_defense_beaming] expected jink plan to be valid.\n";
      ++fails;
    } else {
      // By construction, the pure jink direction should have essentially no X component.
      if (std::fabs(jink.dirWorld.x) > 1e-6) {
        std::cerr << "[test_missile_defense_beaming] expected jink x ~ 0. got x=" << jink.dirWorld.x << "\n";
        ++fails;
      }
      const double len = jink.dirWorld.length();
      if (!approx(len, 1.0, 1e-6)) {
        std::cerr << "[test_missile_defense_beaming] expected normalized jink dir. len=" << len << "\n";
        ++fails;
      }
    }

    sim::MissileEvasionParams pBeam{};
    pBeam.enforceLateralToLos = true;
    pBeam.enableRadarBeaming = true;
    pBeam.radarBeamBlend = 1.0;               // strong bias
    pBeam.radarBeamEngageNotchMultiple = 0.0; // always engage (test only)

    const sim::MissileEvasionPlan beam = sim::planMissileEvasion(m, tgtPos, tgtVel, seed, pBeam);
    if (!beam.valid) {
      std::cerr << "[test_missile_defense_beaming] expected beaming plan to be valid.\n";
      ++fails;
    } else {
      // Beaming bias should introduce an X component in this setup.
      if (!(beam.dirWorld.x > 0.02)) {
        std::cerr << "[test_missile_defense_beaming] expected beaming to introduce +x. got x=" << beam.dirWorld.x << "\n";
        ++fails;
      }
      const double len = beam.dirWorld.length();
      if (!approx(len, 1.0, 1e-6)) {
        std::cerr << "[test_missile_defense_beaming] expected normalized beaming dir. len=" << len << "\n";
        ++fails;
      }
    }
  }

  return fails;
}
