#include "stellar/sim/EncounterDirector.h"

#include "stellar/core/Hash.h"

#include <algorithm>

namespace stellar::sim {

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

  // Pirates occasionally (baseline threat).
  if (ctx.timeDays >= state.nextPirateSpawnDays && counts.alivePirates < 4 && counts.aliveTotal < 14) {
    state.nextPirateSpawnDays = ctx.timeDays + (state.rng.range(120.0, 220.0) / 86400.0); // every ~2-4 minutes
    const int cap = std::max(1, std::min({14 - counts.aliveTotal, 4 - counts.alivePirates, 3}));
    return cap;
  }

  return 0;
}

int planTraderSpawn(EncounterDirectorState& state,
                    const EncounterDirectorContext& ctx,
                    const EncounterDirectorCounts& counts) {
  if (!ctx.combatSimEnabled) return 0;

  // Traders / traffic: gives you something to pirate.
  if (ctx.timeDays >= state.nextTraderSpawnDays && counts.aliveTraders < 3 && ctx.stationCount > 0 && counts.aliveTotal < 14) {
    state.nextTraderSpawnDays = ctx.timeDays + (state.rng.range(70.0, 140.0) / 86400.0);
    return 1;
  }

  return 0;
}

int computeDesiredPolice(const EncounterDirectorContext& ctx,
                         const EncounterDirectorCounts& counts) {
  if (ctx.localFactionId == 0) return 0;

  int desiredPolice = 1;
  if (ctx.playerWantedHere) desiredPolice += 2;
  if (ctx.localBountyCr > 1200.0) desiredPolice += 1;
  if (ctx.localBountyCr > 4500.0) desiredPolice += 1;
  if (ctx.policeHeat > 3.0) desiredPolice += 1;
  if (counts.alivePirates > 0) desiredPolice += 1;
  if (ctx.localRep < -25.0) desiredPolice += 1;
  desiredPolice = std::clamp(desiredPolice, 0, 7);
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
  const double baseMinSec = alert ? 12.0 : 55.0;
  const double baseMaxSec = alert ? 24.0 : 95.0;

  // As heat rises, response gets a bit quicker (but stays bounded).
  const double spawnMinSec = std::clamp(baseMinSec - heat * (alert ? 1.1 : 1.8), 4.0, baseMinSec);
  const double spawnMaxSec = std::clamp(baseMaxSec - heat * (alert ? 1.6 : 2.4), 10.0, baseMaxSec);

  if (ctx.timeDays >= state.nextPoliceSpawnDays && counts.alivePolice < out.desiredPolice && counts.aliveTotal < 16) {
    state.nextPoliceSpawnDays = ctx.timeDays + (state.rng.range(spawnMinSec, spawnMaxSec) / 86400.0);
    const int cap = std::max(1, std::min({16 - counts.aliveTotal, out.desiredPolice - counts.alivePolice, 3}));
    out.spawnMaxCount = cap;
  }

  return out;
}

} // namespace stellar::sim
