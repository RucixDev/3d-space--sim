#include "stellar/sim/Law.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>

namespace stellar::sim {

LawProfile lawProfile(core::u64 universeSeed, core::u32 factionId) {
  LawProfile p{};

  if (factionId == 0) {
    // Neutral baseline; still deterministic and bounded.
    p.scanStrictness = 1.0;
    p.corruption = 0.50;
    return p;
  }

  // Seed per (universe, faction) and keep stable across versions.
  core::u64 s = core::hashCombine(universeSeed, static_cast<core::u64>(factionId));
  s = core::hashCombine(s, core::seedFromText("law_profile_v1"));
  core::SplitMix64 rng(s);

  // Authority: 0 = lax, 1 = strict.
  const double authority = rng.nextDouble();

  // Corruption raw: 0 = honest, 1 = very corrupt.
  const double corruptionRaw = rng.nextDouble();

  // Strictness scales scan frequency (used as a multiplier on scan start probability).
  // Keep this bounded so the smuggling loop remains playable.
  p.scanStrictness = 0.70 + 0.95 * authority; // [0.70, 1.65]

  // Corruption inversely correlates with authority, but has some variation.
  // (Strict factions are less likely to offer bribes.)
  p.corruption = (0.25 + 0.75 * corruptionRaw) * (1.05 - 0.70 * authority);
  p.corruption = std::clamp(p.corruption, 0.05, 1.0);

  // Fine schedule. Keep close to legacy defaults but allow meaningful variance.
  p.fineBaseCr = 140.0 + 280.0 * authority;   // [140, 420]
  p.fineRate = 0.45 + 0.55 * authority;       // [0.45, 1.00]

  // Bribes are usually slightly cheaper than paying the fine + losing the goods.
  // Keep them in a reasonable band so they feel like a risky shortcut, not a free pass.
  p.bribeBaseCr = 180.0 + 240.0 * authority;  // [180, 420]
  p.bribeRate = 0.35 + 0.45 * authority;      // [0.35, 0.80]

  // Keep bribe slope under the fine slope so bribe doesn't become strictly worse
  // at large cargo values (since it's also a corruption check / gameplay choice).
  p.bribeRate = std::min(p.bribeRate, p.fineRate * 0.92);

  // Rep penalties: stricter factions penalize more for the same illegal value.
  p.repBase = 1.8 + 2.8 * authority;          // [1.8, 4.6]
  p.repDiv = 4600.0 - 2100.0 * authority;     // [2500, 4600]
  p.repMin = std::max(1.0, std::round(p.repBase));
  p.repMax = 9.0 + 3.0 * authority;           // [9, 12]

  // Evading a scan demand is always worse.
  p.evadeRepMult = 1.25 + 0.25 * authority;   // [1.25, 1.50]

  return p;
}

} // namespace stellar::sim
