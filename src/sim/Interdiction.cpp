#include "stellar/sim/Interdiction.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static math::Vec3d safeNormalize(const math::Vec3d& v, const math::Vec3d& fallback) {
  const double lsq = v.lengthSq();
  if (lsq < 1.0e-18) return fallback;
  return v / std::sqrt(lsq);
}

static math::Vec3d randomUnit(core::SplitMix64& rng) {
  // Rejection sample-ish without a loop: sample in cube, normalize.
  // This is good enough for gameplay randomness.
  math::Vec3d v{rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)};
  if (v.lengthSq() < 1.0e-12) v = {1.0, 0.0, 0.0};
  return v.normalized();
}

static math::Vec3d rotateRodrigues(const math::Vec3d& v, const math::Vec3d& axisUnit, double angleRad) {
  const double c = std::cos(angleRad);
  const double s = std::sin(angleRad);
  // v' = v*c + (axis x v)*s + axis*(axis·v)*(1-c)
  return v * c + math::cross(axisUnit, v) * s + axisUnit * (math::dot(axisUnit, v) * (1.0 - c));
}

double interdictionCloseness01(double pirateDistanceKm, double maxRangeKm) {
  const double r = std::max(1.0, maxRangeKm);
  const double c = 1.0 - (pirateDistanceKm / r);
  return std::clamp(c, 0.0, 1.0);
}

static double cargoRisk01(double cargoValueCr, double riskValueCr) {
  const double denom = std::max(1.0, riskValueCr);
  return std::clamp(cargoValueCr / denom, 0.0, 2.0);
}

double interdictionTriggerChance(double dtSec,
                                 double pirateCloseness01,
                                 double cargoValueCr,
                                 double pirateStrength,
                                 const InterdictionTriggerParams& p) {
  const double dt = std::max(0.0, dtSec);
  const double close = std::clamp(pirateCloseness01, 0.0, 1.0);
  const double strength = std::clamp(pirateStrength, 0.0, 3.0);
  const double risk = cargoRisk01(std::max(0.0, cargoValueCr), p.cargoRiskValueCr);

  const double closeTerm = p.closenessRatePerSec * std::pow(close, std::max(0.1, p.closenessPower));
  const double cargoTerm = p.cargoRatePerSec * risk;

  // Stronger pirates apply more pressure and can initiate more interdictions.
  const double strengthMul = 0.75 + 0.55 * std::clamp(strength, 0.0, 2.0);
  const double rate = std::max(0.0, (p.baseRatePerSec + closeTerm + cargoTerm) * strengthMul);

  // Poisson process -> per-step probability.
  const double chance = 1.0 - std::exp(-rate * dt);
  return std::clamp(chance, 0.0, p.maxChancePerStep);
}

bool rollInterdictionStart(core::SplitMix64& rng,
                           double dtSec,
                           double pirateCloseness01,
                           double cargoValueCr,
                           double pirateStrength,
                           const InterdictionTriggerParams& params) {
  const double p = interdictionTriggerChance(dtSec, pirateCloseness01, cargoValueCr, pirateStrength, params);
  return rng.chance(p);
}

static math::Vec3d makeEscapeDir(core::SplitMix64& rng,
                                const math::Vec3d& shipForwardWorld,
                                const math::Vec3d& shipUpWorld,
                                const InterdictionParams& params) {
  const math::Vec3d fwd = safeNormalize(shipForwardWorld, {0.0, 0.0, 1.0});
  math::Vec3d up = safeNormalize(shipUpWorld, {0.0, 1.0, 0.0});

  // Build an orthonormal basis (right, upOrtho, fwd).
  math::Vec3d right = math::cross(up, fwd);
  right = safeNormalize(right, {1.0, 0.0, 0.0});
  math::Vec3d upOrtho = math::cross(fwd, right);
  upOrtho = safeNormalize(upOrtho, {0.0, 1.0, 0.0});

  const double minCos = std::clamp(params.minInitialCos, -1.0, 1.0);
  const double maxCos = std::clamp(params.maxInitialCos, -1.0, 1.0);
  const double t = rng.range(std::min(minCos, maxCos), std::max(minCos, maxCos));

  const double theta = rng.range(0.0, 2.0 * M_PI);
  const math::Vec3d lateral = (right * std::cos(theta) + upOrtho * std::sin(theta)).normalized();

  // Construct a unit vector with dot(fwd, dir)=t.
  const double latMag = std::sqrt(std::max(0.0, 1.0 - t * t));
  math::Vec3d dir = (fwd * t + lateral * latMag).normalized();

  // Ensure it's in front (positive dot) so the HUD indicator stays sensible.
  if (math::dot(dir, fwd) < 0.05) {
    dir = (dir + fwd * 1.25).normalized();
  }
  return dir;
}

InterdictionState beginInterdiction(core::u64 universeSeed,
                                    core::u64 pirateId,
                                    double timeDays,
                                    const math::Vec3d& shipForwardWorld,
                                    const math::Vec3d& shipUpWorld,
                                    double pirateCloseness01,
                                    double pirateStrength,
                                    const InterdictionParams& params) {
  InterdictionState s;
  s.phase = InterdictionPhase::Warning;
  s.warningRemainingSec = std::max(0.0, params.warningSec);
  s.activeRemainingSec = 0.0;
  s.pirateId = pirateId;

  // Deterministic seed: stable within a given day for the same pirate.
  const core::u64 day = (core::u64)std::max(0.0, std::floor(timeDays));
  const core::u64 seed = core::hashCombine(core::hashCombine(universeSeed, pirateId), day);
  s.rng.reseed(core::hashCombine(seed, core::fnv1a64("interdiction")));

  // The initial meter is shaped by pirate strength and closeness.
  const double close = std::clamp(pirateCloseness01, 0.0, 1.0);
  const double strength = std::clamp(pirateStrength, 0.0, 3.0);
  const double meter = params.startMeter + (strength - 1.0) * 0.06 + close * 0.05;
  s.meter = std::clamp(meter, 0.05, 0.95);

  s.escapeDir = makeEscapeDir(s.rng, shipForwardWorld, shipUpWorld, params);
  return s;
}

static void driftEscapeDir(InterdictionState& s,
                           double dtSec,
                           const math::Vec3d& shipForwardWorld,
                           double pirateStrength,
                           const InterdictionParams& params) {
  const double dt = std::max(0.0, dtSec);
  if (dt <= 0.0) return;

  const double strength = std::clamp(pirateStrength, 0.0, 3.0);
  const double base = params.escapeDriftRadPerSec * (0.75 + 0.65 * std::clamp(strength, 0.0, 2.0));
  const double jitter = params.escapeJitterRadPerSec * std::clamp(strength, 0.0, 2.0);
  double angle = (base + s.rng.range(-jitter, jitter)) * dt;
  angle = std::clamp(angle, -params.escapeMaxStepRad, params.escapeMaxStepRad);

  if (std::abs(angle) < 1e-9) return;

  const math::Vec3d dir = safeNormalize(s.escapeDir, {0.0, 0.0, 1.0});
  // Choose a perpendicular rotation axis.
  math::Vec3d axis = math::cross(dir, randomUnit(s.rng));
  axis = safeNormalize(axis, math::cross(dir, {1.0, 0.0, 0.0}));
  axis = safeNormalize(axis, math::cross(dir, {0.0, 1.0, 0.0}));

  math::Vec3d next = rotateRodrigues(dir, axis, angle).normalized();

  // Keep it generally in front of the player so the minigame stays readable.
  const math::Vec3d fwd = safeNormalize(shipForwardWorld, {0.0, 0.0, 1.0});
  if (math::dot(next, fwd) < params.minForwardCos) {
    next = (next + fwd * 1.5).normalized();
  }
  s.escapeDir = next;
}

InterdictionStepOutput stepInterdiction(InterdictionState& s,
                                        double dtSec,
                                        const math::Vec3d& shipForwardWorld,
                                        double pirateCloseness01,
                                        double pirateStrength,
                                        bool submit,
                                        const InterdictionParams& params) {
  InterdictionStepOutput out;
  out.phase = s.phase;
  out.warningRemainingSec = s.warningRemainingSec;
  out.activeRemainingSec = s.activeRemainingSec;
  out.meter = s.meter;
  out.escapeDir = s.escapeDir;

  if (s.phase == InterdictionPhase::None) return out;

  const double dt = std::max(0.0, dtSec);

  if (submit) {
    out.submitted = true;
    out.endedThisFrame = true;
    s = InterdictionState{};
    out.phase = InterdictionPhase::None;
    return out;
  }

  if (s.phase == InterdictionPhase::Warning) {
    s.warningRemainingSec = std::max(0.0, s.warningRemainingSec - dt);
    if (s.warningRemainingSec <= 0.0) {
      s.phase = InterdictionPhase::Active;
      s.activeRemainingSec = std::max(0.0, params.activeSec);
      out.beganActive = true;
    }
  }

  if (s.phase == InterdictionPhase::Active) {
    s.activeRemainingSec = std::max(0.0, s.activeRemainingSec - dt);

    driftEscapeDir(s, dt, shipForwardWorld, pirateStrength, params);

    const math::Vec3d fwd = safeNormalize(shipForwardWorld, {0.0, 0.0, 1.0});
    const math::Vec3d esc = safeNormalize(s.escapeDir, {0.0, 0.0, 1.0});

    const double align = std::clamp(math::dot(fwd, esc), -1.0, 1.0);
    out.alignment = align;

    double gain01 = 0.0;
    if (align > params.alignStartCos) {
      const double denom = std::max(1e-6, params.alignFullCos - params.alignStartCos);
      gain01 = std::clamp((align - params.alignStartCos) / denom, 0.0, 1.0);
    }
    out.gain01 = gain01;

    const double close = std::clamp(pirateCloseness01, 0.0, 1.0);
    const double strength = std::clamp(pirateStrength, 0.0, 3.0);
    const double pull = (params.pullBasePerSec + params.pullExtraPerSec * close) * strength;
    out.pullPerSec = pull;

    const double delta = (gain01 * params.playerGainPerSec - pull) * dt;
    s.meter = std::clamp(s.meter + delta, 0.0, 1.0);

    if (s.meter >= 1.0) {
      out.evaded = true;
      out.endedThisFrame = true;
      s = InterdictionState{};
    } else if (s.meter <= 0.0 || s.activeRemainingSec <= 0.0) {
      out.failed = true;
      out.endedThisFrame = true;
      s = InterdictionState{};
    }
  }

  out.phase = s.phase;
  out.warningRemainingSec = s.warningRemainingSec;
  out.activeRemainingSec = s.activeRemainingSec;
  out.meter = s.meter;
  out.escapeDir = s.escapeDir;
  return out;
}

} // namespace stellar::sim
