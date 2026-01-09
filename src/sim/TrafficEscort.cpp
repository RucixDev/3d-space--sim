#include "stellar/sim/TrafficEscort.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double roundTo(double x, double step) {
  if (step <= 0.0) return x;
  return std::round(x / step) * step;
}

TrafficEscortPlan planTrafficEscortContract(core::u64 universeSeed,
                                            SystemId systemId,
                                            core::u64 convoyId,
                                            double timeDays,
                                            double timeToArriveDays,
                                            double cargoValueCr,
                                            double piracy01,
                                            double security01,
                                            double contest01,
                                            bool piratesPresent) {
  TrafficEscortPlan p{};

  piracy01 = std::clamp(piracy01, 0.0, 1.0);
  security01 = std::clamp(security01, 0.0, 1.0);
  contest01 = std::clamp(contest01, 0.0, 1.0);
  cargoValueCr = std::max(0.0, cargoValueCr);

  // Day-stable randomness.
  const core::u64 day = (core::u64)std::max(0.0, std::floor(timeDays));
  const core::u64 seed = core::hashCombine(
    core::hashCombine(universeSeed ^ core::fnv1a64("traffic_escort_contract"), (core::u64)systemId),
    core::hashCombine(convoyId, day));
  core::SplitMix64 rng(seed);

  // Normalize cargo value into a 0..1-ish score.
  const double value01 = std::clamp(cargoValueCr / 50000.0, 0.0, 1.0);

  // Risk: piracy pressure + valuable cargo, reduced by effective security.
  double risk = 0.14 + 0.65 * piracy01 + 0.35 * value01 - 0.55 * security01 + 0.15 * contest01;
  if (piratesPresent) risk += 0.15;
  risk = std::clamp(risk, 0.0, 1.0);
  p.risk01 = risk;

  // Some extreme edge cases should not offer contracts (e.g. no time remaining).
  const double tta = std::max(0.0, timeToArriveDays);
  if (tta <= 0.0) {
    p.offer = false;
    p.durationDays = 0.0;
    p.maxRangeKm = 0.0;
    p.rewardCr = 0.0;
    p.bonusPerPirateCr = 0.0;
    p.repReward = 0.0;
    return p;
  }

  // Range requirement: higher risk => tighter formation.
  p.maxRangeKm = std::clamp(180000.0 - 65000.0 * risk, 110000.0, 190000.0);

  // Duration requirement: scale with risk, then add a small per-convoy variation.
  double durSec = 90.0 + 360.0 * risk;
  durSec *= rng.range(0.85, 1.15);
  durSec = std::clamp(durSec, 75.0, 600.0);

  // Never require longer than the remaining trip.
  const double maxDurSec = tta * 86400.0;
  durSec = std::min(durSec, maxDurSec);
  p.durationDays = durSec / 86400.0;

  // Reward: depends on cargo value, risk, and duration.
  const double durMul = std::clamp(durSec / 240.0, 0.55, 1.55);
  double reward = (180.0 + cargoValueCr * (0.012 + 0.045 * risk) + 650.0 * risk) * durMul;
  if (piratesPresent) reward *= 1.10;
  reward = std::clamp(reward, 200.0, 9000.0);
  p.rewardCr = roundTo(reward, 10.0);

  // Bonus per pirate kill (defense vouchers). Kept moderate to avoid farming.
  double bonus = 70.0 + 180.0 * risk;
  if (piratesPresent) bonus *= 1.15;
  bonus = std::clamp(bonus, 40.0, 500.0);
  p.bonusPerPirateCr = roundTo(bonus, 5.0);

  // Reputation reward (0.5 increments, capped).
  double rep = 0.8 + 2.7 * risk;
  if (piratesPresent) rep += 0.3;
  rep = std::clamp(rep, 0.5, 4.0);
  p.repReward = roundTo(rep, 0.5);

  return p;
}

} // namespace stellar::sim
