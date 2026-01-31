#include "stellar/sim/TrajectoryPredictor.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::sim {
namespace {

static constexpr double kSecondsPerDay = 86400.0;

inline double clampD(double v, double lo, double hi) {
  if (!std::isfinite(v)) return lo;
  return std::clamp(v, lo, hi);
}

struct BurnProfile {
  bool active{false};
  double startSec{0.0};
  double endSec{0.0};
  double nodeTimeSec{0.0}; // requested maneuver node time (for sampling)
  math::Vec3d accelWorldKmS2{0, 0, 0};
};

static BurnProfile buildBurnProfile(const ManeuverNode* node,
                                   double horizonSec,
                                   math::Vec3d& ioStartVelKmS) {
  BurnProfile out{};
  if (!node) return out;

  const double accel = node->accelKmS2;
  if (!(accel > 0.0) || !std::isfinite(accel)) return out;

  const double dvMag = node->deltaVKmS.length();
  if (!(dvMag > 1e-12)) return out;

  double duration = dvMag / accel;
  if (!(duration > 1e-9) || !std::isfinite(duration)) return out;

  out.nodeTimeSec = node->timeSec;
  out.accelWorldKmS2 = node->deltaVKmS * (accel / dvMag);

  // Match ManeuverComputer semantics: start time uses a centered-burn lead
  // (duration/2) plus an extra lead time.
  out.startSec = node->timeSec - (duration * 0.5 + node->extraLeadTimeSec);
  out.endSec = out.startSec + duration;

  // If the burn is entirely before the start of the prediction horizon, best-effort
  // apply the full delta-v to the initial velocity.
  if (out.endSec <= 0.0) {
    ioStartVelKmS += node->deltaVKmS;
    return BurnProfile{};
  }

  // If the burn began before t=0, apply the already-completed fraction to the
  // initial velocity and clamp the remaining burn to start at t=0.
  if (out.startSec < 0.0) {
    const double fracDone = std::clamp((-out.startSec) / duration, 0.0, 1.0);
    ioStartVelKmS += node->deltaVKmS * fracDone;

    const double remainMag = dvMag * (1.0 - fracDone);
    if (!(remainMag > 1e-12)) return BurnProfile{};

    duration = remainMag / accel;
    out.startSec = 0.0;
    out.endSec = out.startSec + duration;
  }

  if (out.startSec >= horizonSec) return BurnProfile{};

  out.startSec = clampD(out.startSec, 0.0, horizonSec);
  out.endSec = clampD(out.endSec, 0.0, horizonSec);
  if (out.endSec <= out.startSec + 1e-9) return BurnProfile{};

  out.active = true;
  return out;
}

static math::Vec3d thrustForInterval(const BurnProfile& burn, double t0Sec, double t1Sec) {
  if (!burn.active) return {0, 0, 0};
  const double mid = 0.5 * (t0Sec + t1Sec);
  if (mid >= burn.startSec && mid <= burn.endSec) return burn.accelWorldKmS2;
  return {0, 0, 0};
}

static math::Vec3d accelAt(const StarSystem& sys,
                           double startTimeDays,
                           double tSec,
                           const math::Vec3d& posKm,
                           const TrajectoryPredictParams& params,
                           const math::Vec3d& thrustAccelKmS2) {
  math::Vec3d a = thrustAccelKmS2;
  if (!params.includeGravity) return a;
  const double tDays = startTimeDays + (tSec / kSecondsPerDay);
  a += systemGravityAccelKmS2(sys, tDays, posKm, params.gravity);
  return a;
}

static void rk4Step(const StarSystem& sys,
                    double startTimeDays,
                    double tSec,
                    double dt,
                    math::Vec3d& posKm,
                    math::Vec3d& velKmS,
                    const TrajectoryPredictParams& params,
                    const math::Vec3d& thrustAccelKmS2) {
  // State: x = pos, v = vel
  // dx/dt = v
  // dv/dt = a(t, x)
  const math::Vec3d a1 = accelAt(sys, startTimeDays, tSec, posKm, params, thrustAccelKmS2);
  const math::Vec3d k1x = velKmS;
  const math::Vec3d k1v = a1;

  const math::Vec3d x2 = posKm + k1x * (dt * 0.5);
  const math::Vec3d v2 = velKmS + k1v * (dt * 0.5);
  const math::Vec3d a2 = accelAt(sys, startTimeDays, tSec + dt * 0.5, x2, params, thrustAccelKmS2);
  const math::Vec3d k2x = v2;
  const math::Vec3d k2v = a2;

  const math::Vec3d x3 = posKm + k2x * (dt * 0.5);
  const math::Vec3d v3 = velKmS + k2v * (dt * 0.5);
  const math::Vec3d a3 = accelAt(sys, startTimeDays, tSec + dt * 0.5, x3, params, thrustAccelKmS2);
  const math::Vec3d k3x = v3;
  const math::Vec3d k3v = a3;

  const math::Vec3d x4 = posKm + k3x * dt;
  const math::Vec3d v4 = velKmS + k3v * dt;
  const math::Vec3d a4 = accelAt(sys, startTimeDays, tSec + dt, x4, params, thrustAccelKmS2);
  const math::Vec3d k4x = v4;
  const math::Vec3d k4v = a4;

  posKm += (dt / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x);
  velKmS += (dt / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
}

} // namespace

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

  // Worst case: 1 sample per step + event samples.
  const int reserve = std::min(maxSamples, 2 + (int)std::ceil(horizon / step) + 8);
  out.reserve((std::size_t)reserve);

  double t = 0.0;
  math::Vec3d pos = startPosKm;
  math::Vec3d vel = startVelKmS;

  const bool haveNode = (node != nullptr) && (node->timeSec >= 0.0) && (node->timeSec <= horizon);
  const bool finiteBurn = haveNode && (node->accelKmS2 > 0.0) && (node->deltaVKmS.lengthSq() > 1e-24);
  const bool impulsive = haveNode && !finiteBurn;

  BurnProfile burn{};
  if (finiteBurn) {
    burn = buildBurnProfile(node, horizon, vel);
  }

  out.push_back({t, pos, vel});

  bool nodeApplied = false;
  if (impulsive && std::abs(node->timeSec) <= 1e-12) {
    vel += node->deltaVKmS;
    nodeApplied = true;
    out.push_back({t, pos, vel});
    if ((int)out.size() >= maxSamples) return out;
  }

  while (t + 1e-9 < horizon && (int)out.size() < maxSamples) {
    double nextEvent = horizon;

    // Impulsive node: always land exactly on node->timeSec.
    if (impulsive && !nodeApplied) {
      const double tn = clampD(node->timeSec, 0.0, horizon);
      if (tn > t + 1e-12) nextEvent = std::min(nextEvent, tn);
    }

    // Burn boundaries + optional sampling at node time (for UI kink markers).
    if (burn.active) {
      if (burn.startSec > t + 1e-12) nextEvent = std::min(nextEvent, burn.startSec);
      if (burn.endSec > t + 1e-12) nextEvent = std::min(nextEvent, burn.endSec);

      const double tNode = clampD(burn.nodeTimeSec, 0.0, horizon);
      if (tNode > t + 1e-12 && tNode < burn.endSec - 1e-12) {
        nextEvent = std::min(nextEvent, tNode);
      }
    } else if (finiteBurn) {
      // The burn may have been clamped away, but still sample the node time if requested.
      const double tNode = clampD(node->timeSec, 0.0, horizon);
      if (tNode > t + 1e-12) nextEvent = std::min(nextEvent, tNode);
    }

    double dt = step;
    dt = std::min(dt, horizon - t);
    dt = std::min(dt, nextEvent - t);
    if (dt <= 1e-12) {
      // We're effectively at an event time already. If this is the impulsive node, apply.
      if (impulsive && !nodeApplied && std::abs(t - node->timeSec) < 1e-9) {
        vel += node->deltaVKmS;
        nodeApplied = true;
        out.push_back({t, pos, vel});
      }
      break;
    }

    const math::Vec3d thrust = thrustForInterval(burn, t, t + dt);
    rk4Step(sys, startTimeDays, t, dt, pos, vel, params, thrust);

    t += dt;

    if (impulsive && !nodeApplied && std::abs(t - node->timeSec) < 1e-6) {
      vel += node->deltaVKmS;
      nodeApplied = true;
    }

    out.push_back({t, pos, vel});
  }

  // Ensure last sample is clamped to horizon.
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
                     const TrajectoryPredictParams& params,
                     const math::Vec3d& thrustAccelKmS2) {
  return {y.vel, accelAt(sys, startTimeDays, tSec, y.pos, params, thrustAccelKmS2)};
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
                                         const TrajectoryPredictParams& params,
                                         const math::Vec3d& thrustAccelKmS2) {
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

  const State k1 = derivAt(sys, startTimeDays, tSec, y, params, thrustAccelKmS2);

  const State y2 = y + (dt * a21) * k1;
  const State k2 = derivAt(sys, startTimeDays, tSec + c2 * dt, y2, params, thrustAccelKmS2);

  const State y3 = y + (dt * a31) * k1 + (dt * a32) * k2;
  const State k3 = derivAt(sys, startTimeDays, tSec + c3 * dt, y3, params, thrustAccelKmS2);

  const State y4 = y + (dt * a41) * k1 + (dt * a42) * k2 + (dt * a43) * k3;
  const State k4 = derivAt(sys, startTimeDays, tSec + c4 * dt, y4, params, thrustAccelKmS2);

  const State y5 = y + (dt * a51) * k1 + (dt * a52) * k2 + (dt * a53) * k3 + (dt * a54) * k4;
  const State k5 = derivAt(sys, startTimeDays, tSec + c5 * dt, y5, params, thrustAccelKmS2);

  const State y6 = y + (dt * a61) * k1 + (dt * a62) * k2 + (dt * a63) * k3 + (dt * a64) * k4 +
                   (dt * a65) * k5;
  const State k6 = derivAt(sys, startTimeDays, tSec + c6 * dt, y6, params, thrustAccelKmS2);

  const State y7 = y + (dt * a71) * k1 + (dt * a72) * k2 + (dt * a73) * k3 + (dt * a74) * k4 +
                   (dt * a75) * k5 + (dt * a76) * k6;
  const State k7 = derivAt(sys, startTimeDays, tSec + c7 * dt, y7, params, thrustAccelKmS2);

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

  // Rough reserve: 1 per output sample + initial + event samples.
  const int reserve = std::min(maxSamples, 2 + (int)std::ceil(horizon / sampleStep) + 16);
  out.reserve((std::size_t)reserve);

  double t = 0.0;
  State y{startPosKm, startVelKmS};

  const bool haveNode = (node != nullptr) && (node->timeSec >= 0.0) && (node->timeSec <= horizon);
  const bool finiteBurn = haveNode && (node->accelKmS2 > 0.0) && (node->deltaVKmS.lengthSq() > 1e-24);
  const bool impulsive = haveNode && !finiteBurn;

  BurnProfile burn{};
  if (finiteBurn) {
    burn = buildBurnProfile(node, horizon, y.vel);
  }

  bool nodeApplied = false;

  out.push_back({t, y.pos, y.vel});

  // Handle an impulsive node at t=0 by emitting a post-burn sample (matches RK4 behavior).
  if (impulsive && std::abs(node->timeSec) <= 1e-12) {
    y.vel += node->deltaVKmS;
    nodeApplied = true;
    out.push_back({t, y.pos, y.vel});
  }

  if (horizon <= 0.0 || (int)out.size() >= maxSamples) return out;

  const double tNodeImp = impulsive ? clampD(node->timeSec, 0.0, horizon) : std::numeric_limits<double>::infinity();
  const double tNodeMark = (haveNode && !impulsive) ? clampD(node->timeSec, 0.0, horizon) : std::numeric_limits<double>::infinity();

  const double tBurnStart = burn.active ? burn.startSec : std::numeric_limits<double>::infinity();
  const double tBurnEnd = burn.active ? burn.endSec : std::numeric_limits<double>::infinity();

  // Output samples are emitted at this cadence (plus burn boundary samples and node samples).
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

  auto integrateTo = [&](double targetT, const math::Vec3d& thrustAccelKmS2) {
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

      const DP45StepResult sr = dormandPrince45Step(sys, startTimeDays, t, h, y, params, thrustAccelKmS2);
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

  enum class EventKind : int {
    Sample = 0,
    NodeImpulsive = 1,
    BurnStart = 2,
    BurnEnd = 3,
    NodeMark = 4,
  };

  while (t + 1e-9 < horizon && (int)out.size() < maxSamples && guard < maxGuard) {
    // Determine the next event time: sample cadence + burn boundaries + optional node sample.
    double eventT = nextSampleT;
    EventKind kind = EventKind::Sample;

    auto take = [&](double candT, EventKind candKind, bool preferOnTie) {
      if (!(candT > t + 1e-12)) return;
      if (candT < eventT - 1e-9 || (std::abs(candT - eventT) < 1e-9 && preferOnTie)) {
        eventT = candT;
        kind = candKind;
      }
    };

    // Earlier boundaries should win over later sample cadence.
    take(tBurnStart, EventKind::BurnStart, false);
    take(tBurnEnd, EventKind::BurnEnd, false);
    take(tNodeMark, EventKind::NodeMark, false);

    // Impulsive node should win ties so the sample is post-burn.
    if (impulsive && !nodeApplied) {
      take(tNodeImp, EventKind::NodeImpulsive, true);
    }

    eventT = clampD(eventT, t, horizon);

    const math::Vec3d thrust = thrustForInterval(burn, t, eventT);
    integrateTo(eventT, thrust);

    if (t + 1e-9 < eventT) {
      // Failed to advance (guard hit). Bail out.
      break;
    }

    const bool doSample = (std::abs(eventT - nextSampleT) < 1e-6);

    bool pushed = false;

    if (kind == EventKind::NodeImpulsive && impulsive && !nodeApplied) {
      y.vel += node->deltaVKmS;
      nodeApplied = true;
      out.push_back({t, y.pos, y.vel});
      pushed = true;
      if ((int)out.size() >= maxSamples) break;
    }

    if (!pushed) {
      out.push_back({t, y.pos, y.vel});
      pushed = true;
      if ((int)out.size() >= maxSamples) break;
    }

    if (doSample) {
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
    } else {
      // If eventT equals the horizon, stop.
      if (eventT >= horizon - 1e-12) break;
    }
  }

  if (!out.empty()) {
    out.back().tSec = std::min(out.back().tSec, horizon);
  }

  return out;
}

} // namespace stellar::sim
