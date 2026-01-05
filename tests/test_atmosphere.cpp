#include "stellar/sim/Atmosphere.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps = 1e-6) { return std::abs(a - b) <= eps; }

int test_atmosphere() {
  int fails = 0;

  // A stationary "Earth-ish" ocean planet at the origin (semiMajorAxisAU=0).
  stellar::sim::StarSystem sys;
  sys.star.massSol = 1.0;
  sys.star.radiusSol = 1.0;
  sys.star.cls = stellar::sim::StarClass::G;

  {
    stellar::sim::Planet p;
    p.name = "TestEarth";
    p.type = stellar::sim::PlanetType::Ocean;
    p.radiusEarth = 1.0;
    p.massEarth = 1.0;
    p.orbit.semiMajorAxisAU = 0.0;
    p.orbit.eccentricity = 0.0;
    p.orbit.periodDays = 365.0;
    p.orbit.meanAnomalyAtEpochRad = 0.0;
    sys.planets.push_back(p);
  }

  const double timeDays = 0.0;
  const double radiusKm = 6371.0;

  stellar::sim::AtmosphereParams params;
  params.densityScale = 1.0;
  params.dragCd = 0.9;
  params.referenceAreaM2 = 80.0;
  params.maxDecelKmS2 = 10.0; // don't clamp in test
  params.applyHeating = true;
  params.heatPerKPaSec = 0.05;
  params.heatingScale = 1.0;

  // --- In-atmosphere sample should produce drag opposing relative velocity ---
  {
    const stellar::math::Vec3d posKm{radiusKm + 10.0, 0.0, 0.0};
    const stellar::math::Vec3d velKmS{0.0, 1.0, 0.0};
    const double massKg = 100000.0;

    const auto s = stellar::sim::sampleSystemAtmosphere(sys, timeDays, posKm, velKmS, massKg, params);
    if (!s.inAtmosphere) {
      std::cerr << "[test_atmosphere] expected inAtmosphere=true\n";
      ++fails;
    } else {
      // Drag should oppose +Y motion -> negative Y acceleration.
      if (!(s.dragAccelKmS2.y < 0.0)) {
        std::cerr << "[test_atmosphere] expected dragAccel.y < 0, got " << s.dragAccelKmS2.y << "\n";
        ++fails;
      }

      // Dynamic pressure should be positive.
      if (!(s.dynamicPressurePa > 0.0)) {
        std::cerr << "[test_atmosphere] expected q>0\n";
        ++fails;
      }

      // Heating should be non-negative.
      if (!(s.heatingHeatPerSec >= 0.0)) {
        std::cerr << "[test_atmosphere] expected heating >= 0\n";
        ++fails;
      }
    }
  }

  // --- High altitude should be vacuum ---
  {
    const stellar::math::Vec3d posKm{radiusKm + 1000.0, 0.0, 0.0};
    const stellar::math::Vec3d velKmS{0.0, 1.0, 0.0};
    const double massKg = 100000.0;

    const auto s = stellar::sim::sampleSystemAtmosphere(sys, timeDays, posKm, velKmS, massKg, params);
    if (s.inAtmosphere) {
      std::cerr << "[test_atmosphere] expected vacuum at high altitude, got inAtmosphere=true\n";
      ++fails;
    }
    if (!near(s.dynamicPressurePa, 0.0, 1e-9) && s.dynamicPressurePa > 0.0) {
      std::cerr << "[test_atmosphere] expected q ~= 0 in vacuum, got " << s.dynamicPressurePa << "\n";
      ++fails;
    }
  }

  // --- Rocky airless worlds: most should return no atmosphere ---
  {
    stellar::sim::Planet rocky;
    rocky.name = "TinyRock";
    rocky.type = stellar::sim::PlanetType::Rocky;
    rocky.radiusEarth = 0.35;
    rocky.massEarth = 0.15;
    const auto m = stellar::sim::atmosphereModelForPlanet(rocky);
    if (m.hasAtmosphere) {
      std::cerr << "[test_atmosphere] expected tiny rocky world to be airless\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_atmosphere] pass\n";
  return fails;
}
