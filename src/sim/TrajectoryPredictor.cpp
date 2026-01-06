#include "stellar/sim/TrajectoryPredictor.h"

#include <algorithm>
#include <cmath>
#include <limits>

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

namespace {

struct State {
  math::Vec3d pos{0, 0, 0};
  math::Vec3d vel{0, 0, 0};
};

inline State operator+(const State& a, const State& b) { return {a.pos + b.pos, a.vel + b.vel}; }
inline State operator-(const State& a, const State& b) { return {a.pos - b.pos, a.vel - b.vel}; }
inline State operator*(const State& a, double s) { return {a.pos * s, a.vel * s}; }
inline State operator*(double s, const State& a) { return a * s; }

static State derivAt(const StarSystem& sys,
                     double startTimeDays,
                     double tSec,
                     const State& y,
                     const TrajectoryPredictParams& params) {
  return {y.vel, accelAt(sys, startTimeDays, tSec, y.pos, params)};
}

struct DP45StepResult {
  State y5;   // 5th order solution (used for propagation)
  State err;  // y5 - y4 (embedded error estimate)
};

static DP45StepResult dormandPrince45Step(const StarSystem& sys,
                                         double startTimeDays,
                                         double tSec,
                                         double dt,
                                         const State& y,
                                         const TrajectoryPredictParams& params) {
  // Dormand–Prince RK5(4)7 coefficients.
  // Butcher tableau reference:
  //   https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
  constexpr double c2 = 1.0 / 5.0;
  constexpr double c3 = 3.0 / 10.0;
  constexpr double c4 = 4.0 / 5.0;
  constexpr double c5 = 8.0 / 9.0;
  constexpr double c6 = 1.0;
  constexpr double c7 = 1.0;

  constexpr double a21 = 1.0 / 5.0;

  constexpr double a31 = 3.0 / 40.0;
  constexpr double a32 = 9.0 / 40.0;

  constexpr double a41 = 44.0 / 45.0;
  constexpr double a42 = -56.0 / 15.0;
  constexpr double a43 = 32.0 / 9.0;

  constexpr double a51 = 19372.0 / 6561.0;
  constexpr double a52 = -25360.0 / 2187.0;
  constexpr double a53 = 64448.0 / 6561.0;
  constexpr double a54 = -212.0 / 729.0;

  constexpr double a61 = 9017.0 / 3168.0;
  constexpr double a62 = -355.0 / 33.0;
  constexpr double a63 = 46732.0 / 5247.0;
  constexpr double a64 = 49.0 / 176.0;
  constexpr double a65 = -5103.0 / 18656.0;

  constexpr double a71 = 35.0 / 384.0;
  constexpr double a72 = 0.0;
  constexpr double a73 = 500.0 / 1113.0;
  constexpr double a74 = 125.0 / 192.0;
  constexpr double a75 = -2187.0 / 6784.0;
  constexpr double a76 = 11.0 / 84.0;

  // 5th order weights (b)
  constexpr double b1 = 35.0 / 384.0;
  constexpr double b2 = 0.0;
  constexpr double b3 = 500.0 / 1113.0;
  constexpr double b4 = 125.0 / 192.0;
  constexpr double b5 = -2187.0 / 6784.0;
  constexpr double b6 = 11.0 / 84.0;
  constexpr double b7 = 0.0;

  // 4th order weights (b*)
  constexpr double bs1 = 5179.0 / 57600.0;
  constexpr double bs2 = 0.0;
  constexpr double bs3 = 7571.0 / 16695.0;
  constexpr double bs4 = 393.0 / 640.0;
  constexpr double bs5 = -92097.0 / 339200.0;
  constexpr double bs6 = 187.0 / 2100.0;
  constexpr double bs7 = 1.0 / 40.0;

  const State k1 = derivAt(sys, startTimeDays, tSec, y, params);

  const State y2 = y + (dt * a21) * k1;
  const State k2 = derivAt(sys, startTimeDays, tSec + c2 * dt, y2, params);

  const State y3 = y + (dt * a31) * k1 + (dt * a32) * k2;
  const State k3 = derivAt(sys, startTimeDays, tSec + c3 * dt, y3, params);

  const State y4 = y + (dt * a41) * k1 + (dt * a42) * k2 + (dt * a43) * k3;
  const State k4 = derivAt(sys, startTimeDays, tSec + c4 * dt, y4, params);

  const State y5 = y + (dt * a51) * k1 + (dt * a52) * k2 + (dt * a53) * k3 + (dt * a54) * k4;
  const State k5 = derivAt(sys, startTimeDays, tSec + c5 * dt, y5, params);

  const State y6 = y + (dt * a61) * k1 + (dt * a62) * k2 + (dt * a63) * k3 + (dt * a64) * k4 +
                   (dt * a65) * k5;
  const State k6 = derivAt(sys, startTimeDays, tSec + c6 * dt, y6, params);

  const State y7 = y + (dt * a71) * k1 + (dt * a72) * k2 + (dt * a73) * k3 + (dt * a74) * k4 +
                   (dt * a75) * k5 + (dt * a76) * k6;
  const State k7 = derivAt(sys, startTimeDays, tSec + c7 * dt, y7, params);

  const State sol5 = y + dt * (b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6 + b7 * k7);
  const State sol4 = y + dt * (bs1 * k1 + bs2 * k2 + bs3 * k3 + bs4 * k4 + bs5 * k5 + bs6 * k6 + bs7 * k7);

  return {sol5, sol5 - sol4};
}

static double errorNorm(const State& y,
                        const State& yNew,
                        const State& err,
                        const TrajectoryPredictParams& params) {
  // A simple scalar error norm: max(position_err/pos_scale, velocity_err/vel_scale).
  // This keeps the behavior deterministic and avoids per-component branching.
  const double posScale = std::max(1e-12, params.absTolPosKm + params.relTol * std::max(y.pos.length(), yNew.pos.length()));
  const double velScale = std::max(1e-12, params.absTolVelKmS + params.relTol * std::max(y.vel.length(), yNew.vel.length()));

  const double ePos = err.pos.length() / posScale;
  const double eVel = err.vel.length() / velScale;
  return std::max(ePos, eVel);
}

} // namespace

std::vector<TrajectorySample> predictTrajectoryRK45Adaptive(const StarSystem& sys,
                                                           double startTimeDays,
                                                           const math::Vec3d& startPosKm,
                                                           const math::Vec3d& startVelKmS,
                                                           const TrajectoryPredictParams& params,
                                                           const ManeuverNode* node) {
  std::vector<TrajectorySample> out;

  const double horizon = std::max(0.0, params.horizonSec);
  const double sampleStep = std::max(1e-6, params.stepSec);
  const int maxSamples = std::max(2, params.maxSamples);

  const double minStep = std::max(1e-9, params.minStepSec);
  const double maxStep = std::max(minStep, params.maxStepSec);

  // Rough reserve: 1 per output sample + initial + optional node.
  const int reserve = std::min(maxSamples, 2 + (int)std::ceil(horizon / sampleStep) + 4);
  out.reserve((std::size_t)reserve);

  double t = 0.0;
  State y{startPosKm, startVelKmS};

  out.push_back({t, y.pos, y.vel});
  if (horizon <= 0.0 || (int)out.size() >= maxSamples) return out;

  const bool wantNode = (node != nullptr) && (node->timeSec >= 0.0) && (node->timeSec <= horizon);
  const double tNode = wantNode ? node->timeSec : -1.0;
  bool nodeApplied = false;

  // Output samples are emitted at this cadence (plus a sample at the node time).
  double nextSampleT = std::min(sampleStep, horizon);

  // Internal step guess (will be adapted).
  double dt = std::clamp(sampleStep, minStep, maxStep);

  // Standard adaptive RK safety factors.
  constexpr double safety = 0.9;
  constexpr double minFactor = 0.2;
  constexpr double maxFactor = 5.0;
  constexpr double exponent = 0.2; // 1/5 for RK5(4)

  int guard = 0;
  const int maxGuard = 2'000'000;

  auto integrateTo = [&](double targetT) {
    while (t + 1e-12 < targetT && guard < maxGuard) {
      const double remaining = targetT - t;
      if (remaining <= 0.0) {
        t = targetT;
        break;
      }

      // Proposed step.
      double h = std::min(dt, remaining);

      // If we're not in the final sliver of the segment, respect minStep.
      if (remaining > minStep) {
        h = std::max(h, minStep);
      }

      h = std::clamp(h, 0.0, maxStep);
      if (h <= 0.0) {
        t = targetT;
        break;
      }

      const DP45StepResult sr = dormandPrince45Step(sys, startTimeDays, t, h, y, params);
      const State yNew = sr.y5;
      const double errN = errorNorm(y, yNew, sr.err, params);

      const bool accept = (errN <= 1.0) || (h <= minStep) || (remaining <= minStep);

      if (accept) {
        y = yNew;
        t += h;

        // Update dt for the next attempt.
        double factor = maxFactor;
        if (errN > 0.0) {
          factor = safety * std::pow(1.0 / errN, exponent);
          factor = std::clamp(factor, minFactor, maxFactor);
        }
        dt = std::clamp(h * factor, minStep, maxStep);
      } else {
        // Reject: reduce step.
        double factor = safety * std::pow(1.0 / errN, exponent);
        factor = std::clamp(factor, minFactor, 1.0);
        dt = std::clamp(h * factor, minStep, maxStep);
      }

      ++guard;
    }

    // Snap to the target if we're extremely close (helps keep exact sample times).
    if (std::abs(t - targetT) < 1e-9) t = targetT;
  };

  while (t + 1e-9 < horizon && (int)out.size() < maxSamples && guard < maxGuard) {
    // Determine the next event time: either the next output sample, or the maneuver node.
    double eventT = nextSampleT;
    bool doSample = true;
    bool doNode = false;

    if (wantNode && !nodeApplied) {
      if (tNode <= eventT + 1e-9) {
        eventT = std::clamp(tNode, t, horizon);
        doNode = true;
        doSample = (std::abs(eventT - nextSampleT) < 1e-6);
      }
    }

    eventT = std::clamp(eventT, t, horizon);
    integrateTo(eventT);

    if (t + 1e-9 < eventT) {
      // Failed to advance (guard hit). Bail out.
      break;
    }

    bool pushed = false;

    if (doNode && wantNode && !nodeApplied) {
      y.vel += node->deltaVKmS;
      nodeApplied = true;
      out.push_back({t, y.pos, y.vel});
      pushed = true;
      if ((int)out.size() >= maxSamples) break;
    }

    if (doSample) {
      if (!pushed) {
        out.push_back({t, y.pos, y.vel});
        if ((int)out.size() >= maxSamples) break;
      }

      // Advance to next sample.
      if (nextSampleT + 1e-9 < horizon) {
        nextSampleT = std::min(nextSampleT + sampleStep, horizon);
      } else {
        nextSampleT = horizon;
      }

      // If numerical drift caused nextSampleT to slip behind, bump it forward.
      if (nextSampleT <= t + 1e-12) {
        nextSampleT = std::min(t + sampleStep, horizon);
      }
    }
  }

  if (!out.empty()) {
    out.back().tSec = std::min(out.back().tSec, horizon);
  }

  return out;
}

} // namespace stellar::sim
