#include "stellar/math/Fov.h"

#include "test_harness.h"

#include <cmath>

// Small sanity checks for FOV conversions and distance-dependent projection helpers.

int test_fov_math() {
  int failures = 0;
  using stellar::math::degToRad;
  using stellar::math::radToDeg;

  // --- Vertical <-> horizontal conversions (aspect = 16:9) ---
  {
    const double aspect = 16.0 / 9.0;
    const double fovYDeg = 90.0;

    const double fovXRad = stellar::math::fov::verticalToHorizontalRad(degToRad(fovYDeg), aspect);
    const double fovYRad2 = stellar::math::fov::horizontalToVerticalRad(fovXRad, aspect);

    // Round-trip stability.
    CHECK(std::abs(radToDeg(fovYRad2) - fovYDeg) < 1.0e-10);

    // Known value: fovX = 2*atan(tan(45deg)*16/9) ~= 120.51deg
    CHECK(std::abs(radToDeg(fovXRad) - 120.510) < 0.02);
  }

  // --- Dolly zoom (keep projected size constant) ---
  {
    const double refDist = 10.0;
    const double refFov = 60.0;
    const double dist = 20.0;

    const double fovDeg = stellar::math::fov::dollyZoomFovDeg(refFov, refDist, dist);

    // tan(fov/2) should scale by (refDist/dist)
    const double lhs = std::tan(degToRad(fovDeg) * 0.5);
    const double rhs = std::tan(degToRad(refFov) * 0.5) * (refDist / dist);
    CHECK(std::abs(lhs - rhs) < 1.0e-12);
  }

  // --- Angular diameter ---
  {
    // radius=1 at distance=2 -> 2*asin(0.5) = 60deg.
    const double angDeg = stellar::math::fov::angularDiameterDeg(1.0, 2.0);
    CHECK(std::abs(angDeg - 60.0) < 1.0e-10);
  }

  // --- Physical lens equivalence ---
  {
    // 24mm sensor height, 50mm focal length.
    const double fovRad = stellar::math::fov::fovRadFromFocalLengthMm(50.0, 24.0);
    const double fovDeg = radToDeg(fovRad);

    // 2*atan(24/(2*50)) ~= 26.99deg.
    CHECK(std::abs(fovDeg - 26.99) < 0.2);

    // Invert: focal ~= 50mm
    const double focal = stellar::math::fov::focalLengthMmFromFovRad(fovRad, 24.0);
    CHECK(std::abs(focal - 50.0) < 1.0e-6);
  }

  return failures;
}
