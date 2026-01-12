#pragma once

#include <algorithm>
#include <cmath>

namespace stellar::render {

// A tiny, GL-agnostic auto-exposure helper for HDR post-processing.
//
// Feed this with an average scene luminance (linear) to compute a smooth
// exposure multiplier. A typical pipeline:
//  1) measure average luminance from the HDR scene (often log-average)
//  2) compute target exposure that maps avg luminance to `key`
//  3) smooth exposure over time with separate brighten/darken speeds
struct AutoExposureConfig {
  // Middle-grey target (often ~0.18 in linear space).
  float key{0.18f};

  // Clamp range for the exposure multiplier.
  float minExposure{0.25f};
  float maxExposure{4.0f};

  // Adaptation speeds (1/sec). Higher = quicker adaptation.
  float speedUp{2.5f};
  float speedDown{1.0f};
};

struct AutoExposureState {
  // Current exposure multiplier.
  float exposure{1.0f};

  // Most recent measurements (for UI/debug).
  float avgLuminance{1.0f};
  float avgLogLuminance{0.0f};

  bool valid{false};
};

// Compute the ideal exposure multiplier given an average luminance.
inline float computeTargetExposure(float avgLuminance, const AutoExposureConfig& cfg) {
  const float lum = std::max(1e-6f, avgLuminance);
  const float target = cfg.key / lum;
  return std::clamp(target, cfg.minExposure, cfg.maxExposure);
}

// Update `state` given a measured average luminance and a timestep.
// Returns the updated exposure multiplier.
float stepAutoExposure(AutoExposureState& state,
                       float avgLuminance,
                       float avgLogLuminance,
                       float dtSeconds,
                       const AutoExposureConfig& cfg);

} // namespace stellar::render
