#pragma once

#include "stellar/core/Types.h"

#include <algorithm>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// SensorModel — lightweight, headless sensor / radar detection helpers
// -----------------------------------------------------------------------------
//
// The game has an in-flight "radar" HUD overlay. Historically it behaved like a
// perfect omniscient list of objects inside a fixed range. This module provides
// a tiny, deterministic signal-strength model so UI/AI can reason about:
//   - range falloff (inverse-square)
//   - occlusion attenuation
//   - smoothing / hysteresis for identification
//
// It intentionally avoids dependencies on renderer code, SDL, or global state.

struct SensorParams {
  // Range at which sensorPower=1 and signature=1 yields strength=0.5.
  double halfRangeKm{200000.0};

  // Attenuation applied to SNR when the target is occluded by a large body.
  // 0 => fully blocked, 1 => no attenuation.
  double occlusionAtten{0.08};

  // Prevent singularities at extremely close range.
  double minDistKm{1.0};
};

// Compute an instantaneous detection strength in [0,1].
//
// The model treats (sensorPower * signature) as a normalized "gain" term and
// applies an inverse-square falloff relative to params.halfRangeKm, mapping an
// SNR-like quantity into [0,1] with: strength = snr / (1 + snr).
double computeSensorStrength01(double distKm,
                               double signature,
                               double sensorPower,
                               bool occluded,
                               const SensorParams& params = {});

struct SensorTrackParams {
  // Response half-life for rising/falling strength (seconds).
  double riseHalfLifeSec{0.30};
  double fallHalfLifeSec{0.60};

  // Visibility thresholds.
  double ghostThreshold{0.18};        // show a blip when strength >= ghostThreshold
  double identifyThreshold{0.65};     // lock/identify when strength >= identifyThreshold

  // Hysteresis: once identified, keep identified until falling below this threshold.
  double maintainIdentifyThreshold{0.45};
};

struct SensorTrack {
  // Filtered strength in [0,1].
  double strength01{0.0};

  // Sticky identification with hysteresis.
  bool identified{false};
};

struct SensorTrackResult {
  double strength01{0.0};
  bool visible{false};
  bool identified{false};
};

// Update a persistent sensor track using exponential smoothing and
// identification hysteresis.
SensorTrackResult updateSensorTrack(SensorTrack& track,
                                    double dtSec,
                                    double measuredStrength01,
                                    const SensorTrackParams& params = {});

} // namespace stellar::sim
