#pragma once

#include "stellar/core/Random.h"
#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Interdiction (supercruise tether minigame)
// -----------------------------------------------------------------------------
// A lightweight, deterministic minigame used during supercruise.
//
// Goals:
//  - Keep the math out of the SDL prototype so behavior is testable.
//  - Allow game code to scale difficulty by pirate strength/closeness.
//  - Provide enough state for a HUD indicator (escape vector + meter).
//
// The module does not decide *when* interdictions happen; it only provides:
//  - a helper to compute a per-frame trigger chance from simple scalars
//  - begin/step helpers to run the minigame once triggered

enum class InterdictionPhase : core::u8 {
  None = 0,
  Warning = 1,
  Active = 2,
};

struct InterdictionTriggerParams {
  // Baseline hazard rate (events/sec) when closeness=0 and cargoValue=0.
  double baseRatePerSec{0.08};

  // Additional hazard when a pirate is close.
  // Rate contribution is `closenessRatePerSec * pow(closeness01, closenessPower)`.
  double closenessRatePerSec{0.75};
  double closenessPower{2.0};

  // Additional hazard from valuable cargo. cargoRisk01 is computed as
  // clamp(cargoValueCr / cargoRiskValueCr, 0, 2).
  double cargoRatePerSec{0.22};
  double cargoRiskValueCr{9000.0};

  // Clamp per-step probability to keep spikes under control when dt is large.
  double maxChancePerStep{0.85};
};

struct InterdictionParams {
  // Phase durations (real seconds).
  double warningSec{2.2};
  double activeSec{11.0};

  // Starting tug-of-war meter (0..1). Higher means the victim is already
  // partially "caught" when the tether engages.
  double startMeter{0.42};

  // Alignment model: gain ramps from 0 at alignStartCos to 1 at alignFullCos.
  double alignStartCos{0.72};
  double alignFullCos{0.94};
  double playerGainPerSec{0.78};

  // Pirate pull model (meter/sec). pull = (pullBase + pullExtra * closeness01) * pirateStrength.
  double pullBasePerSec{0.26};
  double pullExtraPerSec{0.22};

  // Escape vector drift while active. Makes the minigame less "set and forget".
  // Drift is scaled by pirateStrength.
  double escapeDriftRadPerSec{0.16};
  double escapeJitterRadPerSec{0.10};
  double escapeMaxStepRad{0.22};

  // When the escape vector ends up behind the player, nudge it forward.
  double minForwardCos{0.20};

  // Initial escape vector: pick dot(forward, dir) in [minInitialCos, maxInitialCos].
  double minInitialCos{0.35};
  double maxInitialCos{0.85};
};

struct InterdictionState {
  InterdictionPhase phase{InterdictionPhase::None};
  double warningRemainingSec{0.0};
  double activeRemainingSec{0.0};
  double meter{0.0};

  // World-space escape direction (unit vector). The player should keep
  // their forward vector aligned to this direction.
  math::Vec3d escapeDir{0.0, 0.0, 1.0};

  core::u64 pirateId{0};
  core::SplitMix64 rng{0};
};

struct InterdictionStepOutput {
  bool evaded{false};
  bool failed{false};
  bool submitted{false};
  bool beganActive{false};
  bool endedThisFrame{false};

  double alignment{-1.0}; // dot(shipForward, escapeDir)
  double gain01{0.0};
  double pullPerSec{0.0};

  // Snapshot (for HUD)
  InterdictionPhase phase{InterdictionPhase::None};
  double warningRemainingSec{0.0};
  double activeRemainingSec{0.0};
  double meter{0.0};
  math::Vec3d escapeDir{0.0, 0.0, 1.0};
};

// Convert a pirate distance to a normalized closeness scalar.
//  - distance >= maxRangeKm => 0
//  - distance <= 0 => 1
double interdictionCloseness01(double pirateDistanceKm, double maxRangeKm = 240000.0);

// Compute a per-step probability of starting an interdiction.
// The caller can use this to roll a warning phase start.
double interdictionTriggerChance(double dtSec,
                                 double pirateCloseness01,
                                 double cargoValueCr,
                                 double pirateStrength,
                                 const InterdictionTriggerParams& params = {});

// Convenience wrapper that samples `rng`.
bool rollInterdictionStart(core::SplitMix64& rng,
                           double dtSec,
                           double pirateCloseness01,
                           double cargoValueCr,
                           double pirateStrength,
                           const InterdictionTriggerParams& params = {});

// Begin an interdiction minigame.
// The generated escape direction is deterministic from (universeSeed, pirateId, dayStamp)
// and the ship's current basis (forward/up).
InterdictionState beginInterdiction(core::u64 universeSeed,
                                    core::u64 pirateId,
                                    double timeDays,
                                    const math::Vec3d& shipForwardWorld,
                                    const math::Vec3d& shipUpWorld,
                                    double pirateCloseness01,
                                    double pirateStrength,
                                    const InterdictionParams& params = {});

// Advance the interdiction state.
//
// Inputs:
//  - shipForwardWorld: unit-ish forward direction of the victim
//  - pirateCloseness01: 0..1 (stronger pull when closer)
//  - pirateStrength: >=0 (scales pull and escape drift)
//  - submit: if true, ends immediately with `submitted=true`
InterdictionStepOutput stepInterdiction(InterdictionState& state,
                                        double dtSec,
                                        const math::Vec3d& shipForwardWorld,
                                        double pirateCloseness01,
                                        double pirateStrength,
                                        bool submit,
                                        const InterdictionParams& params = {});

inline bool interdictionInProgress(const InterdictionState& s) {
  return s.phase != InterdictionPhase::None;
}

} // namespace stellar::sim
