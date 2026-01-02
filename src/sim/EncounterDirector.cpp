#include "stellar/sim/EncounterDirector.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

EncounterDirectorState makeEncounterDirector(core::u64 universeSeed, double timeDays) {
  EncounterDirectorState st{};

  // Mix the universe seed with a stable salt so encounter cadence is stable
  // even if other systems reseed their own RNGs with the same universe seed.
  const core::u64 salt = core::fnv1a64("EncounterDirector");
  st.rng.reseed(core::hashCombine(universeSeed, salt));

  // "Soon after start" defaults similar to the prototype's inline values.
  st.nextPirateSpawnDays = timeDays + 0.01;
  st.nextTraderSpawnDays = timeDays + 0.008;
  st.nextPoliceSpawnDays = timeDays + 0.006;

  return st;
}

int planPirateSpawn(EncounterDirectorState& state,
                    const EncounterDirectorContext& ctx,
                    const EncounterDirectorCounts& counts) {
  if (!ctx.combatSimEnabled) return 0;

  const double piracy = clamp01(ctx.piracy01);
  const double sec = clamp01(ctx.security01);

  // Max pirates shifts with system danger/security, but keeps the old baseline.
  const int maxPirates = std::clamp((int)std::round(4.0
                                                  + 4.0 * (piracy - 0.5)
                                                  - 3.0 * (sec - 0.5)),
                                   1, 8);

  // Pirates occasionally (baseline threat).
  if (ctx.timeDays >= state.nextPirateSpawnDays && counts.alivePirates < maxPirates && counts.aliveTotal < 14) {
    // Frequency multiplier: higher piracy => faster spawns, higher security => slower spawns.
    // Clamp so the loop stays playable and bounded.
    const double freq = std::clamp(1.0 + 1.1 * (piracy - 0.5) - 0.9 * (sec - 0.5), 0.35, 2.5);

    double minSec = 120.0 / freq;
    double maxSec = 220.0 / freq;

    // Keep bounds sane.
    minSec = std::clamp(minSec, 35.0, 420.0);
    maxSec = std::clamp(maxSec, 70.0, 520.0);
    if (maxSec < minSec) maxSec = minSec;

    state.nextPirateSpawnDays = ctx.timeDays + (state.rng.range(minSec, maxSec) / 86400.0); // every ~2-4 minutes (baseline)
    const int cap = std::max(1, std::min({14 - counts.aliveTotal, maxPirates - counts.alivePirates, 3}));
    return cap;
  }

  return 0;
}

int planTraderSpawn(EncounterDirectorState& state,
                    const EncounterDirectorContext& ctx,
                    const EncounterDirectorCounts& counts) {
  if (!ctx.combatSimEnabled) return 0;

  const double traffic = clamp01(ctx.traffic01);
  const double piracy = clamp01(ctx.piracy01);

  const int maxTraders = std::clamp((int)std::round(3.0
                                                   + 3.0 * (traffic - 0.5)
                                                   - 2.0 * (piracy - 0.5)),
                                   0, 7);

  // Traders / traffic: gives you something to pirate.
  if (ctx.timeDays >= state.nextTraderSpawnDays && counts.aliveTraders < maxTraders && ctx.stationCount > 0 && counts.aliveTotal < 14) {
    const double freq = std::clamp(1.0 + 1.0 * (traffic - 0.5) - 0.6 * (piracy - 0.5), 0.35, 2.5);

    double minSec = 70.0 / freq;
    double maxSec = 140.0 / freq;

    minSec = std::clamp(minSec, 30.0, 360.0);
    maxSec = std::clamp(maxSec, 60.0, 520.0);
    if (maxSec < minSec) maxSec = minSec;

    state.nextTraderSpawnDays = ctx.timeDays + (state.rng.range(minSec, maxSec) / 86400.0);
    return 1;
  }

  return 0;
}

int computeDesiredPolice(const EncounterDirectorContext& ctx,
                         const EncounterDirectorCounts& counts) {
  if (ctx.localFactionId == 0) return 0;

  const double sec = clamp01(ctx.security01);
  const double traffic = clamp01(ctx.traffic01);

  int desiredPolice = 1;
  if (ctx.playerWantedHere) desiredPolice += 2;
  if (ctx.localBountyCr > 1200.0) desiredPolice += 1;
  if (ctx.localBountyCr > 4500.0) desiredPolice += 1;
  if (ctx.policeHeat > 3.0) desiredPolice += 1;
  if (counts.alivePirates > 0) desiredPolice += 1;
  if (ctx.localRep < -25.0) desiredPolice += 1;

  // System-level context: stable/high-security/high-traffic systems tend to have
  // more patrol presence.
  desiredPolice += (int)std::round((sec - 0.5) * 2.0);
  desiredPolice += (int)std::round((traffic - 0.5) * 1.0);

  // If there is a governing faction, keep a minimal baseline patrol presence.
  // (Independent space is handled by the early return above.)
  desiredPolice = std::clamp(desiredPolice, 1, 7);
  return desiredPolice;
}

PoliceSpawnPlan planPoliceSpawn(EncounterDirectorState& state,
                                const EncounterDirectorContext& ctx,
                                const EncounterDirectorCounts& counts) {
  PoliceSpawnPlan out{};

  if (ctx.localFactionId == 0) return out;

  out.desiredPolice = computeDesiredPolice(ctx, counts);

  if (!ctx.combatSimEnabled) return out;

  // If recently alerted by a crime, tighten the spawn interval.
  const bool alert = ctx.policeAlert;
  const double heat = ctx.policeHeat;
  const double sec = clamp01(ctx.security01);
  const double baseMinSec = alert ? 12.0 : 55.0;
  const double baseMaxSec = alert ? 24.0 : 95.0;

  // As heat rises, response gets a bit quicker (but stays bounded).
  const double spawnMinSec = std::clamp(baseMinSec - heat * (alert ? 1.1 : 1.8), 4.0, baseMinSec);
  const double spawnMaxSec = std::clamp(baseMaxSec - heat * (alert ? 1.6 : 2.4), 10.0, baseMaxSec);

  // System security tweaks response time (bounded, default 1.0 at sec=0.5).
  const double respMul = std::clamp(1.0 - 0.60 * (sec - 0.5), 0.55, 1.45);
  const double spawnMinAdjSec = std::clamp(spawnMinSec * respMul, 3.0, 180.0);
  const double spawnMaxAdjSec = std::clamp(spawnMaxSec * respMul, 6.0, 300.0);

  if (ctx.timeDays >= state.nextPoliceSpawnDays && counts.alivePolice < out.desiredPolice && counts.aliveTotal < 16) {
    const double lo = std::min(spawnMinAdjSec, spawnMaxAdjSec);
    const double hi = std::max(spawnMinAdjSec, spawnMaxAdjSec);
    state.nextPoliceSpawnDays = ctx.timeDays + (state.rng.range(lo, hi) / 86400.0);
    const int cap = std::max(1, std::min({16 - counts.aliveTotal, out.desiredPolice - counts.alivePolice, 3}));
    out.spawnMaxCount = cap;
  }

  return out;
}

} // namespace stellar::sim
