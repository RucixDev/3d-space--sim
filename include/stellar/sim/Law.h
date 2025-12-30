#pragma once

#include "stellar/core/Types.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

// A lightweight, deterministic per-faction "law profile".
//
// The core sim already has contraband legality masks (see Contraband.h). This module
// adds a *severity* layer so different factions can feel meaningfully different:
//   - Strict factions scan more often and issue steeper fines.
//   - More corrupt factions are more likely to offer bribes.
//
// Design goals:
//   - Deterministic across runs (seed + factionId)
//   - Stable, bounded values (no extreme multipliers)
//   - Usable in both headless/tests and the SDL prototype
struct LawProfile {
  // Multiplies scan start rate (1.0 = baseline).
  double scanStrictness{1.0};

  // 0..1 multiplier for whether a police scan might offer a bribe.
  // 0 => never; 1 => baseline chance.
  double corruption{0.5};

  // Contraband enforcement: fine = base + illegalValueCr * rate.
  double fineBaseCr{200.0};
  double fineRate{0.65};

  // Bribe offer (to keep cargo): bribe = base + illegalValueCr * rate.
  double bribeBaseCr{250.0};
  double bribeRate{0.55};

  // Reputation penalty model:
  //   penalty = -clamp(repBase + illegalValueCr / repDiv, repMin, repMax)
  double repBase{2.0};
  double repDiv{3000.0};
  double repMin{2.0};
  double repMax{10.0};

  // Multiplier applied to rep penalties if the player *evades* the scan window.
  double evadeRepMult{1.35};

  // ---- Helpers ----
  double contrabandFineCr(double illegalValueCr) const {
    illegalValueCr = std::max(0.0, illegalValueCr);
    return std::max(0.0, fineBaseCr + illegalValueCr * fineRate);
  }

  double contrabandBribeCr(double illegalValueCr) const {
    illegalValueCr = std::max(0.0, illegalValueCr);
    // Keep legacy behavior: bribes have a small minimum so they don't become trivial.
    return std::max(200.0, bribeBaseCr + illegalValueCr * bribeRate);
  }

  double contrabandRepPenalty(double illegalValueCr) const {
    illegalValueCr = std::max(0.0, illegalValueCr);
    const double p = repBase + illegalValueCr / std::max(1.0, repDiv);
    return -std::clamp(p, repMin, repMax);
  }

  double contrabandEvadeRepPenalty(double illegalValueCr) const {
    illegalValueCr = std::max(0.0, illegalValueCr);
    const double p = repBase * evadeRepMult + illegalValueCr / std::max(1.0, repDiv * 0.85);
    // Slightly looser max bound when evading.
    return -std::clamp(p, repMin * evadeRepMult, repMax * 1.20);
  }
};

// Deterministically generate a per-faction law profile.
//
// Faction 0 is treated as "Independent" and returns a neutral baseline profile.
LawProfile lawProfile(core::u64 universeSeed, core::u32 factionId);

} // namespace stellar::sim
