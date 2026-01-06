#pragma once

#include "stellar/math/Vec3.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Two-body orbital mechanics (state vector -> classical scalars)
// -----------------------------------------------------------------------------
//
// This is intentionally lightweight: enough for debugging UI and simple gameplay
// features (periapsis/apoapsis readouts, safety checks, etc.).
//
// Units:
//   - position: km
//   - velocity: km/s
//   - mu: km^3/s^2

struct TwoBodyOrbit {
  enum class Type : unsigned char {
    Invalid = 0,
    Elliptic,
    Parabolic,
    Hyperbolic,
  };

  double muKm3S2{0.0};

  Type type{Type::Invalid};

  double rKm{0.0};
  double vKmS{0.0};

  // Specific orbital energy (km^2/s^2).
  double specificEnergyKm2S2{0.0};

  // Semi-major axis (km). Negative indicates a hyperbolic trajectory.
  double semiMajorAxisKm{0.0};

  // Eccentricity magnitude.
  double eccentricity{0.0};

  // Periapsis / apoapsis distances from the central body center (km).
  // For unbound trajectories, apoapsisKm = 0.
  double periapsisKm{0.0};
  double apoapsisKm{0.0};

  // Orbital period (seconds) for bound ellipses. 0 for unbound.
  double periodSec{0.0};

  // Phase / timing.
  // For elliptic orbits:
  //  - trueAnomalyRad is normalized to [0, 2pi)
  //  - meanAnomalyRad is normalized to [0, 2pi)
  //  - timeToPeriapsisSec/timeToApoapsisSec refer to the *next* occurrence
  // For hyperbolic orbits:
  //  - trueAnomalyRad is in (-pi, pi)
  //  - timeSincePeriapsisSec is signed (negative before periapsis)
  //  - timeToPeriapsisSec is max(0, -timeSincePeriapsisSec)
  double trueAnomalyRad{0.0};
  double eccentricAnomalyRad{0.0};
  double meanAnomalyRad{0.0};
  double meanMotionRadPerSec{0.0};
  double timeSincePeriapsisSec{0.0};
  double timeToPeriapsisSec{0.0};
  double timeToApoapsisSec{0.0};

  // Velocity decomposition at the current state (body-centric two-body frame).
  double radialSpeedKmS{0.0};
  double tangentialSpeedKmS{0.0};

  // Debug vectors
  math::Vec3d angularMomentumKm2S{0, 0, 0};
  double angularMomentumMagKm2S{0.0};
  math::Vec3d eccentricityVec{0, 0, 0};
};

// Compute basic two-body orbital parameters from a state vector.
TwoBodyOrbit solveTwoBodyOrbit(const math::Vec3d& relPosKm,
                              const math::Vec3d& relVelKmS,
                              double muKm3S2);

} // namespace stellar::sim
