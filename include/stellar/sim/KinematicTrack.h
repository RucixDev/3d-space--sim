#pragma once

#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// KinematicTrack — tiny constant-velocity tracker (alpha-beta filter)
// -----------------------------------------------------------------------------
//
// Purpose:
//   Radar/UI often wants a plausible, deterministic estimate of a contact's
//   position/velocity that can coast briefly when measurements drop out.
//
// This is not a full Kalman filter. It is an alpha-beta (g-h) filter:
//   predict:   p += v * dt
//   residual:  r = z - p
//   correct:   p += alpha * r
//              v += (beta / dt) * r
//
// We also maintain a simple scalar uncertainty (sigmaKm) that grows over time
// and shrinks on measurement updates. Callers can use it to draw an
// "uncertainty ring" or to modulate UI behavior.
//
// Design constraints:
//   - deterministic and header-only
//   - no dynamic allocations
//   - stable under varying dt (guards for tiny dt)

struct KinematicTrackParams {
  // Gains are interpolated by measurement strength in [0,1].
  // Larger alpha makes position snap harder to measurements.
  double minAlpha{0.20};
  double maxAlpha{0.85};

  // Larger beta makes velocity converge faster (but can be noisier).
  double minBeta{0.002};
  double maxBeta{0.14};

  // Clamp the estimated velocity magnitude (safety against pathological updates).
  double maxSpeedKmS{25000.0};

  // Uncertainty model (scalar 1-sigma in km).
  double sigmaInitKm{8.0};
  double sigmaMinKm{0.03};
  double sigmaMaxKm{800000.0};

  // Process noise: sigma grows by this amount every second.
  double sigmaGrowthKmPerSec{0.35};

  // Measurement update: sigma *= lerp(1.0, sigmaShrinkFactor, strength01).
  // Lower values shrink uncertainty faster for strong measurements.
  double sigmaShrinkFactor{0.55};
};

struct KinematicTrack3d {
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
  double sigmaKm{0.0};

  // Seconds since last measurement update.
  double ageSinceMeasSec{0.0};

  bool initialized{false};
};

struct KinematicTrackResult {
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
  double sigmaKm{0.0};
  double ageSinceMeasSec{0.0};

  bool hasMeasurement{false};
};

inline double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

inline double lerp(double a, double b, double t) {
  return a + (b - a) * t;
}

inline KinematicTrackResult updateKinematicTrack(KinematicTrack3d& track,
                                                 double dtSec,
                                                 bool hasMeasurement,
                                                 const math::Vec3d& measPosKm,
                                                 double measStrength01,
                                                 const KinematicTrackParams& params = {}) {
  if (!std::isfinite(dtSec) || dtSec < 0.0) dtSec = 0.0;

  measStrength01 = clamp01(measStrength01);

  // Predict.
  if (track.initialized && dtSec > 0.0) {
    track.posKm += track.velKmS * dtSec;
  }

  // Uncertainty always grows a bit (process noise).
  if (track.initialized && params.sigmaGrowthKmPerSec > 0.0 && dtSec > 0.0) {
    track.sigmaKm += params.sigmaGrowthKmPerSec * dtSec;
  }

  if (track.initialized && dtSec > 0.0) {
    track.ageSinceMeasSec += dtSec;
  } else if (!track.initialized && dtSec > 0.0) {
    // Keep the age counter meaningful even before init (used by UI).
    track.ageSinceMeasSec += dtSec;
  }

  // Correct.
  if (hasMeasurement) {
    if (!track.initialized) {
      track.initialized = true;
      track.posKm = measPosKm;
      track.velKmS = {0, 0, 0};
      track.sigmaKm = std::clamp(params.sigmaInitKm, params.sigmaMinKm, params.sigmaMaxKm);
      track.ageSinceMeasSec = 0.0;
    } else {
      const double alpha = std::clamp(lerp(params.minAlpha, params.maxAlpha, measStrength01), 0.0, 1.0);
      const double beta = std::clamp(lerp(params.minBeta, params.maxBeta, measStrength01), 0.0, 1.0);

      const math::Vec3d residual = measPosKm - track.posKm;
      track.posKm += residual * alpha;

      // Only update velocity when dt is sane.
      if (dtSec > 1.0e-6) {
        track.velKmS += residual * (beta / dtSec);
        if (params.maxSpeedKmS > 0.0 && std::isfinite(params.maxSpeedKmS)) {
          track.velKmS = math::clampMagnitude(track.velKmS, params.maxSpeedKmS);
        }
      }

      // Shrink uncertainty in proportion to measurement strength.
      const double shrink = std::clamp(lerp(1.0, params.sigmaShrinkFactor, measStrength01), 0.0, 1.0);
      track.sigmaKm *= shrink;
      track.sigmaKm = std::clamp(track.sigmaKm, params.sigmaMinKm, params.sigmaMaxKm);

      track.ageSinceMeasSec = 0.0;
    }
  }

  // Final clamps.
  if (track.initialized) {
    track.sigmaKm = std::clamp(track.sigmaKm, params.sigmaMinKm, params.sigmaMaxKm);
  }

  KinematicTrackResult out{};
  out.posKm = track.posKm;
  out.velKmS = track.velKmS;
  out.sigmaKm = track.sigmaKm;
  out.ageSinceMeasSec = track.ageSinceMeasSec;
  out.hasMeasurement = hasMeasurement;
  return out;
}

} // namespace stellar::sim
