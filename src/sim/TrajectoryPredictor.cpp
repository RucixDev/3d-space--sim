#include "stellar/sim/TrajectoryPredictor.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static math::Vec3d accelAt(const StarSystem& sys,
                           double startTimeDays,
                           double tSec,
                           const math::Vec3d& posKm,
                           const TrajectoryPredictParams& params) {
  if (!params.includeGravity) return {0, 0, 0};
  const double tDays = startTimeDays + (tSec / 86400.0);
  return systemGravityAccelKmS2(sys, tDays, posKm, params.gravity);
}

static void rk4Step(const StarSystem& sys,
                    double startTimeDays,
                    double tSec,
                    double dt,
                    math::Vec3d& posKm,
                    math::Vec3d& velKmS,
                    const TrajectoryPredictParams& params) {
  // State: x = pos, v = vel
  // dx/dt = v
  // dv/dt = a(t, x)
  const math::Vec3d a1 = accelAt(sys, startTimeDays, tSec, posKm, params);
  const math::Vec3d k1x = velKmS;
  const math::Vec3d k1v = a1;

  const math::Vec3d x2 = posKm + k1x * (dt * 0.5);
  const math::Vec3d v2 = velKmS + k1v * (dt * 0.5);
  const math::Vec3d a2 = accelAt(sys, startTimeDays, tSec + dt * 0.5, x2, params);
  const math::Vec3d k2x = v2;
  const math::Vec3d k2v = a2;

  const math::Vec3d x3 = posKm + k2x * (dt * 0.5);
  const math::Vec3d v3 = velKmS + k2v * (dt * 0.5);
  const math::Vec3d a3 = accelAt(sys, startTimeDays, tSec + dt * 0.5, x3, params);
  const math::Vec3d k3x = v3;
  const math::Vec3d k3v = a3;

  const math::Vec3d x4 = posKm + k3x * dt;
  const math::Vec3d v4 = velKmS + k3v * dt;
  const math::Vec3d a4 = accelAt(sys, startTimeDays, tSec + dt, x4, params);
  const math::Vec3d k4x = v4;
  const math::Vec3d k4v = a4;

  posKm += (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
  velKmS += (dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
}

std::vector<TrajectorySample> predictTrajectoryRK4(const StarSystem& sys,
                                                   double startTimeDays,
                                                   const math::Vec3d& startPosKm,
                                                   const math::Vec3d& startVelKmS,
                                                   const TrajectoryPredictParams& params,
                                                   const ManeuverNode* node) {
  std::vector<TrajectorySample> out;

  const double horizon = std::max(0.0, params.horizonSec);
  const double step = std::max(1e-6, params.stepSec);
  const int maxSamples = std::max(2, params.maxSamples);

  // Worst case: 1 sample per step + initial sample.
  const int reserve = std::min(maxSamples, 2 + (int)std::ceil(horizon / step));
  out.reserve((std::size_t)reserve);

  double t = 0.0;
  math::Vec3d pos = startPosKm;
  math::Vec3d vel = startVelKmS;

  out.push_back({t, pos, vel});

  const bool wantNode = (node != nullptr) && (node->timeSec >= 0.0) && (node->timeSec <= horizon);
  bool nodeApplied = false;

  while (t + 1e-9 < horizon && (int)out.size() < maxSamples) {
    double dt = std::min(step, horizon - t);
    if (dt <= 0.0) break;

    // If a maneuver node falls inside this step, split the step so the burn occurs at
    // the exact requested time.
    if (wantNode && !nodeApplied) {
      const double tn = node->timeSec;
      if (tn >= t - 1e-9 && tn <= t + dt + 1e-9) {
        const double dt1 = std::clamp(tn - t, 0.0, dt);
        const double dt2 = dt - dt1;

        if (dt1 > 1e-9) {
          rk4Step(sys, startTimeDays, t, dt1, pos, vel, params);
          t += dt1;
        } else {
          // No pre-burn integration needed; we're essentially at the node already.
          t = tn;
        }

        // Apply burn.
        vel += node->deltaVKmS;
        nodeApplied = true;

        // Emit a sample at the node time (post-burn), so render code can draw a visible kink.
        out.push_back({t, pos, vel});
        if ((int)out.size() >= maxSamples) break;

        if (dt2 > 1e-9 && t + 1e-9 < horizon) {
          rk4Step(sys, startTimeDays, t, dt2, pos, vel, params);
          t += dt2;
          out.push_back({t, pos, vel});
        }
        continue;
      }
    }

    rk4Step(sys, startTimeDays, t, dt, pos, vel, params);
    t += dt;
    out.push_back({t, pos, vel});
  }

  // Ensure last sample is exactly at horizon (within numerical tolerances).
  if (!out.empty()) {
    out.back().tSec = std::min(out.back().tSec, horizon);
  }

  return out;
}

} // namespace stellar::sim
