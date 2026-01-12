#include "stellar/render/AutoExposure.h"

#include <algorithm>
#include <cmath>

namespace stellar::render {

static float clampDt(float dt) {
  if (!std::isfinite(dt)) return 0.0f;
  // Clamp to avoid huge time jumps (pause/resume) causing instant exposure snaps.
  return std::clamp(dt, 0.0f, 0.25f);
}

float stepAutoExposure(AutoExposureState& state,
                       float avgLuminance,
                       float avgLogLuminance,
                       float dtSeconds,
                       const AutoExposureConfig& cfg) {
  // Keep last measurements for debug/UI.
  if (std::isfinite(avgLuminance)) state.avgLuminance = std::max(1e-6f, avgLuminance);
  if (std::isfinite(avgLogLuminance)) state.avgLogLuminance = avgLogLuminance;

  const float dt = clampDt(dtSeconds);
  const float target = computeTargetExposure(state.avgLuminance, cfg);

  if (!state.valid || !std::isfinite(state.exposure) || state.exposure <= 0.0f) {
    state.exposure = target;
    state.valid = true;
    return state.exposure;
  }

  const float speed = (target > state.exposure) ? cfg.speedUp : cfg.speedDown;
  const float s = std::max(0.0f, speed);

  // Exponential smoothing: alpha = 1 - exp(-dt * speed).
  const float alpha = 1.0f - std::exp(-dt * s);
  state.exposure = state.exposure + (target - state.exposure) * alpha;
  state.exposure = std::clamp(state.exposure, cfg.minExposure, cfg.maxExposure);

  return state.exposure;
}

} // namespace stellar::render
