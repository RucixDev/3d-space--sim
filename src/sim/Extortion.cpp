#include "stellar/sim/Extortion.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double clamp01(double x) { return std::clamp(x, 0.0, 1.0); }

ExtortionDemandPlan planExtortionDemand(core::u64 universeSeed,
                                        core::u64 targetId,
                                        double targetCargoValueCr,
                                        double attackerStrength,
                                        double defenderStrength,
                                        double systemSecurity01,
                                        bool policeNearby) {
  ExtortionDemandPlan out{};

  // Sanitize inputs.
  if (!std::isfinite(targetCargoValueCr)) targetCargoValueCr = 0.0;
  if (!std::isfinite(attackerStrength)) attackerStrength = 0.0;
  if (!std::isfinite(defenderStrength)) defenderStrength = 0.0;
  if (!std::isfinite(systemSecurity01)) systemSecurity01 = 0.0;

  targetCargoValueCr = std::max(0.0, targetCargoValueCr);
  attackerStrength = std::max(0.0, attackerStrength);
  defenderStrength = std::max(0.0, defenderStrength);
  systemSecurity01 = clamp01(systemSecurity01);

  // If there's effectively no cargo, don't offer extortion.
  if (!(targetCargoValueCr > 90.0)) {
    return out;
  }

  // Deterministic local RNG to add a bit of variety without making tests flaky.
  // The seed is stable for a given (universeSeed, targetId).
  const core::u64 kSalt = core::fnv1a64("extortion_plan_v1");
  core::SplitMix64 rng(core::hashCombine(universeSeed ^ kSalt, targetId));

  // Relative threat: only the ratio matters.
  const double denom = std::max(1e-6, defenderStrength);
  double ratio = attackerStrength / denom;
  ratio = std::clamp(ratio, 0.15, 12.0);

  // Smoothly map ratio into 0..1.
  // ratio=1 => 0.5, ratio>>1 => ~1, ratio<<1 => ~0.
  const double threat01 = clamp01(0.5 + 0.5 * std::tanh(std::log(ratio)));

  // Demand fraction:
  //  - more threat => demand more
  //  - more security => demand a bit less (more risk)
  double demandFrac = 0.12 + 0.22 * threat01 - 0.08 * systemSecurity01;
  demandFrac *= rng.range(0.90, 1.12);
  demandFrac = std::clamp(demandFrac, 0.10, 0.42);

  double demanded = targetCargoValueCr * demandFrac;
  demanded = std::clamp(demanded, 80.0, 12000.0);
  // Round to a nice UI-friendly multiple.
  demanded = std::round(demanded / 10.0) * 10.0;

  // Compliance:
  //  - rises with threat
  //  - falls with security (confidence in police response)
  //  - falls with nearby authorities (immediate backup)
  //  - rises slightly with larger cargo value (a bigger shipper is more willing to pay to avoid delay)
  const double cargo01 = clamp01(targetCargoValueCr / 35000.0);
  double comply = 0.22 + 0.70 * threat01 - 0.35 * systemSecurity01 + 0.08 * cargo01;
  if (policeNearby) comply -= 0.22;
  comply *= rng.range(0.92, 1.08);
  comply = std::clamp(comply, 0.05, 0.95);

  // Suggested timings:
  //  - higher security => shorter cooldown (more willing to re-engage comms / bluff)
  //  - stronger threat => longer flee (they panic)
  const double cooldownSec = std::clamp(55.0 + 55.0 * systemSecurity01 + rng.range(-10.0, 12.0), 30.0, 140.0);
  const double fleeSec = std::clamp(95.0 + 115.0 * threat01 + rng.range(-20.0, 35.0), 60.0, 260.0);

  out.offer = true;
  out.demandedValueCr = demanded;
  out.complyChance01 = comply;
  out.suggestedCooldownSec = cooldownSec;
  out.suggestedFleeSec = fleeSec;
  return out;
}

} // namespace stellar::sim
