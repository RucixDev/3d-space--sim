#pragma once

#include "stellar/core/Types.h"

#include <algorithm>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// ElectronicWarfare — lightweight, headless radar jamming helpers
// -----------------------------------------------------------------------------
//
// The game contains a HUD radar overlay. Radar strength is now computed via
// sim::SensorModel (see SensorModel.h). Electronic warfare adds a second,
// similarly-shaped inverse-square field: **jammers**.
//
// Design goals:
//  - Deterministic + headless: usable in tests and tools.
//  - Simple scaling law: inverse-square SNR mapped to [0,1] via snr/(1+snr).
//  - No direct knowledge of game entities; callers supply distances and powers.

struct EwJammerParams {
  // Range at which jammerPower=1 yields jamming01 ~= 0.5.
  double halfRangeKm{140000.0};

  // Prevent singularities at extremely close range.
  double minDistKm{1.0};

  // How strongly jamming suppresses sensor power.
  // EffectiveSensorPower = base / (1 + suppressionGain * effectiveJamming01)
  double suppressionGain{1.15};

  // When an active sensor ping is running, reduce the effective jamming by this
  // multiplicative factor (0..1). Lower => ping "punches through" more.
  double pingJammingMult{0.45};
};

// Compute the raw jammer SNR-like value at a given distance.
// This is primarily useful when aggregating multiple jammers.
//
// snr = jammerPower * (halfRange / dist)^2
double computeJammingSnr(double distKm,
                         double jammerPower,
                         const EwJammerParams& params = {});

// Map SNR into [0,1] via snr/(1+snr).
// Safe for non-finite values.
double jamming01FromSnr(double snr);

// Convenience for a single jammer.
inline double computeJamming01(double distKm,
                               double jammerPower,
                               const EwJammerParams& params = {}) {
  return jamming01FromSnr(computeJammingSnr(distKm, jammerPower, params));
}

// Apply a jamming level to a base sensorPower.
//
// pingActive:
//   When true, the jamming is reduced by params.pingJammingMult.
//
// Returns an effective sensor power (>= 0).
double applyJammingToSensorPower(double baseSensorPower,
                                 double jamming01,
                                 bool pingActive,
                                 const EwJammerParams& params = {});

// --- Ghost blips (UI helpers) -------------------------------------------------
//
// Jamming is often experienced as noise/false returns. We provide a deterministic
// ghost blip generator so the HUD can draw plausible noise *without* keeping its
// own RNG state.

struct EwGhostBlip {
  // In the observer's local frame (x right, z forward). y is intentionally
  // omitted to keep the HUD 2D.
  double xKm{0.0};
  double zKm{0.0};

  // Suggested alpha/strength in [0,1].
  double strength01{0.0};

  // Stable id so UI can vary icons.
  core::u64 id{0};
};

struct EwGhostParams {
  // Upper bound for number of ghost blips.
  int maxBlips{22};

  // Strength range for the noise field. The output strength is additionally
  // scaled by jamming01.
  double minStrength01{0.10};
  double maxStrength01{0.55};

  // How quickly ghost blips drift (in km per second) in radar space.
  double driftKmPerSec{1200.0};

  // How frequently the ghost field re-seeds.
  // 0.15 Hz => ~6.7s per pattern.
  double reseedHz{0.15};
};

// Generate a set of ghost blips.
//
// seed:
//   Stable caller-supplied seed (e.g. hash(systemSeed, shipId)).
// timeSec:
//   Real-time seconds (used only for drift + reseed buckets).
// rangeKm:
//   Radar range (used to keep blips inside the display volume).
// jamming01:
//   Aggregate jamming level in [0,1]. If <= 0, returns empty.
std::vector<EwGhostBlip> generateGhostBlips(core::u64 seed,
                                           double timeSec,
                                           double rangeKm,
                                           double jamming01,
                                           const EwGhostParams& params = {});

} // namespace stellar::sim
