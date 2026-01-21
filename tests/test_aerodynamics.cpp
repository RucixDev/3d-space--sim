#include "test_harness.h"

#include "stellar/sim/Aerodynamics.h"
#include "stellar/sim/Atmosphere.h"

#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"

#include <cmath>

using namespace stellar;

int test_aerodynamics() {
  int failures = 0;

  sim::AerodynamicsParams p{};
  p.enabled = true;
  p.minDynamicPressurePa = 25.0;
  p.wingAreaM2 = 120.0;
  p.liftSlopePerRad = 5.0;
  p.clMax = 1.25;
  p.stallAngleDeg = 18.0;
  p.controlSurfaces = false;

  sim::AtmosphereSample atmo{};
  atmo.inAtmosphere = true;
  atmo.dynamicPressurePa = 5000.0;
  atmo.relVelKmS = {0, 0, 10};

  // Neutral AoA -> near-zero lift.
  {
    const auto s = sim::computeAerodynamics(atmo, math::Quatd::identity(), {0, 0, 0}, 10000.0, p);
    CHECK(s.active);
    CHECK(std::abs(s.alphaRad) < 1e-6);
    CHECK(std::abs(s.cl) < 1e-6);
    CHECK(s.liftAccelKmS2.length() < 1e-9);
    CHECK(std::isfinite(s.extraDragAccelKmS2.x));
    CHECK(std::isfinite(s.extraDragAccelKmS2.y));
    CHECK(std::isfinite(s.extraDragAccelKmS2.z));
  }

  // Positive AoA -> positive lift (up).
  {
    atmo.relVelKmS = {0, -1, 10}; // downward component => wind comes from below => alpha > 0
    const auto s = sim::computeAerodynamics(atmo, math::Quatd::identity(), {0, 0, 0}, 10000.0, p);
    CHECK(s.active);
    CHECK(s.alphaRad > 0.0);
    CHECK(s.cl > 0.0);
    CHECK(s.liftAccelKmS2.y > 0.0);
    // Extra drag should oppose forward motion (roughly -Z in world).
    CHECK(s.extraDragAccelKmS2.z < 0.0);
  }

  // Below min dynamic pressure -> inactive.
  {
    atmo.dynamicPressurePa = 1.0;
    const auto s = sim::computeAerodynamics(atmo, math::Quatd::identity(), {0, 0, 0}, 10000.0, p);
    CHECK(!s.active);
  }

  return failures;
}
