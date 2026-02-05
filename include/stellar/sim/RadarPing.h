#pragma once

#include "stellar/core/Types.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// RadarPing — deterministic ping sweep helpers
// -----------------------------------------------------------------------------
//
// This header provides small, deterministic helpers used by the HUD radar ping.
// A ping is represented by a (start time, duration) and a sweep fraction in [0,1]
// that can either be outbound-only or out-and-back (triangle wave).
//
// The core idea is to gate the ping "power boost" by a thin ring around the
// current sweep radius rather than applying a blanket multiplier.
//
// All functions are header-only to keep integration simple.

struct RadarPingParams {
  // If true, the sweep goes out then back (triangle wave: 0→1→0).
  // If false, it only goes outbound (0→1).
  bool outAndBack{true};

  // Ring thickness as a fraction of the radar range [0,1].
  // Example: 0.08 means an 8% range-thick band receives the strongest boost.
  double ringThicknessFrac{0.08};

  // Edge feathering as a fraction of ring half-width [0,1].
  // 0 = sharp edges, 1 = very soft edges.
  double ringFeatherFrac{0.85};
};

inline double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

inline double smoothstep01(double edge0, double edge1, double x) {
  if (edge1 <= edge0) return (x < edge0) ? 0.0 : 1.0;
  double t = (x - edge0) / (edge1 - edge0);
  t = clamp01(t);
  return t * t * (3.0 - 2.0 * t);
}

inline double pingTime01(double nowSec, double startSec, double durationSec) {
  if (!(durationSec > 1.0e-9)) return 1.0;
  const double t01 = (nowSec - startSec) / durationSec;
  return clamp01(t01);
}

inline double pingFrac01(double t01, bool outAndBack) {
  t01 = clamp01(t01);
  if (!outAndBack) return t01;

  // Triangle wave on [0,1]: f(t) = 1 - |2t - 1|
  return 1.0 - std::abs(2.0 * t01 - 1.0);
}

inline double pingFrac01(double nowSec,
                         double startSec,
                         double durationSec,
                         const RadarPingParams& params = {}) {
  const double t01 = pingTime01(nowSec, startSec, durationSec);
  return pingFrac01(t01, params.outAndBack);
}

inline double pingRadiusKm(double pingFrac, double radarRangeKm) {
  if (radarRangeKm <= 0.0) return 0.0;
  return radarRangeKm * clamp01(pingFrac);
}

// Compute a [0,1] boost factor for a given target distance based on the current
// ping sweep ring.
inline double pingRingBoost01(double distKm,
                              double pingFrac,
                              double radarRangeKm,
                              const RadarPingParams& params = {}) {
  if (!(radarRangeKm > 0.0)) return 0.0;
  distKm = std::max(0.0, distKm);

  double thicknessFrac = std::clamp(params.ringThicknessFrac, 0.0, 1.0);
  if (thicknessFrac < 1.0e-6) thicknessFrac = 1.0e-6;

  const double halfWidthKm = 0.5 * thicknessFrac * radarRangeKm;
  const double featherFrac = std::clamp(params.ringFeatherFrac, 0.0, 1.0);
  const double featherKm = halfWidthKm * featherFrac;

  const double radiusKm = pingRadiusKm(pingFrac, radarRangeKm);
  const double dKm = std::abs(distKm - radiusKm);

  const double edge0 = std::max(0.0, halfWidthKm - featherKm);
  const double edge1 = halfWidthKm + featherKm;

  double w = 1.0 - smoothstep01(edge0, edge1, dKm);

  // Mild sharpening so the ring reads as a "return" rather than a wide blanket.
  w = w * w;

  return clamp01(w);
}

inline double pingRingBoost01(double distKm,
                              double nowSec,
                              double startSec,
                              double durationSec,
                              double radarRangeKm,
                              const RadarPingParams& params = {}) {
  const double frac = pingFrac01(nowSec, startSec, durationSec, params);
  return pingRingBoost01(distKm, frac, radarRangeKm, params);
}

} // namespace stellar::sim
