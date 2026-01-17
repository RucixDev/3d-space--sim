#include "stellar/proc/GalaxyMorphology.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

int test_galaxy_morphology() {
  int failures = 0;

  proc::GalaxyParams gp{};
  gp.radiusLy = 1000.0;
  gp.thicknessLy = 1000.0;

  // ---- Bar mask: should be elongated along its major axis. ----
  gp.barStrength = 1.0;
  gp.barAngleDeg = 0.0; // major axis along +X
  gp.barLengthLy = 1000.0;
  gp.barWidthLy = 200.0;
  gp.barPower = 2.0;

  const double barMajor = proc::galaxyBarMask01(gp, 750.0, 0.0);
  const double barMinor = proc::galaxyBarMask01(gp, 0.0, 750.0);
  CHECK(barMajor > 0.0);
  CHECK(barMinor == 0.0);

  // ---- Ring mask: peaks around ringRadius. ----
  gp.ringStrength = 1.0;
  gp.ringRadiusLy = 500.0;
  gp.ringWidthLy = 50.0;
  gp.ringPower = 1.0;

  const double ringPeak = proc::galaxyRingMask01(gp, 500.0);
  const double ringOff = proc::galaxyRingMask01(gp, 700.0);
  CHECK(ringPeak > 0.95);
  CHECK(ringOff < ringPeak);

  // ---- Flare: thickness should increase with radius. ----
  gp.flareStrength = 1.0;
  gp.flarePower = 1.0;

  const double half0 = proc::galaxyThicknessHalfLy(gp, 0.0);
  const double halfR = proc::galaxyThicknessHalfLy(gp, gp.radiusLy);
  CHECK(halfR > half0);
  CHECK(std::abs(halfR - half0 * 2.0) < 1.0e-9);

  // ---- Warp: sinusoidal midplane offset, amplitude rises with radius. ----
  gp.warpAmplitudeLy = 200.0;
  gp.warpStartRadiusLy = 0.0;
  gp.warpPower = 1.0;
  gp.warpLobes = 2;
  gp.warpPhaseDeg = 90.0; // theta=0 => +sin(pi/2)
  gp.warpNoiseStrength = 0.0;
  gp.warpNoiseFreq = 0.0;

  const core::u64 seed = 1234;
  const double w0 = proc::galaxyWarpZLy(seed, gp, 500.0, 0.0);
  const double w90 = proc::galaxyWarpZLy(seed, gp, 0.0, 500.0);
  CHECK(std::abs(w0) <= gp.warpAmplitudeLy + 1.0e-6);
  CHECK(std::abs(w90) <= gp.warpAmplitudeLy + 1.0e-6);
  CHECK(w0 > 0.0);
  CHECK(w90 < 0.0);

  // ---- Sample struct sanity ----
  const auto s = proc::sampleGalaxyMorphology(seed, gp, math::Vec3d{500.0, 0.0, 0.0});
  CHECK(s.densityMul >= 1.0);
  CHECK(s.bar01 >= 0.0 && s.bar01 <= 1.0);
  CHECK(s.ring01 >= 0.0 && s.ring01 <= 1.0);
  CHECK(s.thicknessHalfLy > 0.0);

  return failures ? 1 : 0;
}
