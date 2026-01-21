#pragma once

#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

// Utilities for field-of-view (FOV) conversions and distance-related projection math.
//
// These helpers are intentionally dependency-free so they can be used from:
// - render code (projection setup)
// - camera rigs / cinematic tools
// - gameplay systems (sensor cones, scan ranges, etc.)
// - unit tests
//
// Conventions:
// - aspect = width / height
// - angles are in radians unless the function name includes "Deg"

namespace stellar::math::fov {

namespace {
inline double clampPositive(double v, double fallback) {
  return (v > 1e-12) ? v : fallback;
}
} // namespace

// Convert vertical FOV -> horizontal FOV using:
//   tan(fovX/2) = aspect * tan(fovY/2)
inline double verticalToHorizontalRad(double fovYRad, double aspect) {
  const double a = clampPositive(aspect, 1.0);
  const double t = std::tan(fovYRad * 0.5);
  return 2.0 * std::atan(t * a);
}

// Convert horizontal FOV -> vertical FOV using:
//   tan(fovY/2) = tan(fovX/2) / aspect
inline double horizontalToVerticalRad(double fovXRad, double aspect) {
  const double a = clampPositive(aspect, 1.0);
  const double t = std::tan(fovXRad * 0.5);
  return 2.0 * std::atan(t / a);
}

// Physical camera equivalence helpers.
//
// For a pinhole camera:
//   fov = 2 * atan(sensorExtent / (2 * focalLength))
inline double fovRadFromFocalLengthMm(double focalLengthMm, double sensorExtentMm) {
  const double f = clampPositive(focalLengthMm, 1.0);
  const double s = std::max(1e-12, sensorExtentMm);
  return 2.0 * std::atan(s / (2.0 * f));
}

inline double focalLengthMmFromFovRad(double fovRad, double sensorExtentMm) {
  const double s = std::max(1e-12, sensorExtentMm);
  const double t = std::tan(fovRad * 0.5);
  if (t <= 1e-12) return 1.0e12; // extremely narrow
  return s / (2.0 * t);
}

// Height of the view frustum slice at a given distance for a vertical FOV.
inline double viewHeightAtDistance(double distance, double fovYRad) {
  const double d = std::max(0.0, distance);
  return 2.0 * d * std::tan(fovYRad * 0.5);
}

// Approximate world-units per pixel at the center of the screen.
inline double unitsPerPixelAtDistance(double distance, double fovYRad, double viewportHeightPx) {
  const double h = clampPositive(viewportHeightPx, 1.0);
  return viewHeightAtDistance(distance, fovYRad) / h;
}

// "Dolly zoom" FOV helper.
//
// Keeps a target at the look-at depth at constant screen scale while the camera distance changes.
// In a perspective projection, NDC scale is proportional to (1 / tan(fov/2)) / distance.
// Holding this constant gives tan(fov/2) ∝ 1/distance.
inline double dollyZoomFovRad(double referenceFovYRad, double referenceDistance, double distance) {
  const double d0 = std::max(1e-12, referenceDistance);
  const double d1 = std::max(1e-12, distance);
  const double t0 = std::tan(referenceFovYRad * 0.5);
  const double t1 = t0 * (d0 / d1);
  return 2.0 * std::atan(t1);
}

// Exact angular diameter (radians) of a sphere/disc with radius `radius` viewed from
// distance `distanceToCenter` to its center:
//   diameter = 2 * asin(radius / distanceToCenter)
inline double angularDiameterRad(double radius, double distanceToCenter) {
  const double d = std::max(0.0, distanceToCenter);
  if (d <= 1e-12) return kPi;
  const double x = std::clamp(radius / d, 0.0, 1.0);
  return 2.0 * std::asin(x);
}

inline double verticalToHorizontalDeg(double fovYDeg, double aspect) {
  return radToDeg(verticalToHorizontalRad(degToRad(fovYDeg), aspect));
}

inline double horizontalToVerticalDeg(double fovXDeg, double aspect) {
  return radToDeg(horizontalToVerticalRad(degToRad(fovXDeg), aspect));
}

inline double dollyZoomFovDeg(double referenceFovYDeg, double referenceDistance, double distance) {
  return radToDeg(dollyZoomFovRad(degToRad(referenceFovYDeg), referenceDistance, distance));
}

inline double angularDiameterDeg(double radius, double distanceToCenter) {
  return radToDeg(angularDiameterRad(radius, distanceToCenter));
}

} // namespace stellar::math::fov
