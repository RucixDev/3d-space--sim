#include "stellar/sim/Distress.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double clamp01(double x) { return std::clamp(x, 0.0, 1.0); }

DistressPlan planDistressEncounter(core::u64 universeSeed,
                                   SystemId systemId,
                                   core::u64 signalId,
                                   double timeDays,
                                   core::u32 localFactionId) {
  DistressPlan p{};

  // Keep generation stable across versions.
  core::u64 s = core::hashCombine(universeSeed, core::seedFromText("distress_plan_v1"));
  s = core::hashCombine(s, (core::u64)systemId);
  s = core::hashCombine(s, (core::u64)localFactionId);
  s = core::hashCombine(s, signalId);
  // Mix in the integer day only (so small dt jitter doesn't change content).
  const core::u64 day = (core::u64)std::max(0.0, std::floor(timeDays));
  s = core::hashCombine(s, day);

  core::SplitMix64 rng(s);

  p.payerFactionId = localFactionId;

  // Legitimacy: most distress calls are real, but "pirate bait" exists.
  const bool legit = rng.nextDouble() < 0.64;

  // Ambush chance is higher for illegitimate calls, but not zero for real calls.
  const double ambushChance = legit ? 0.22 : 0.82;
  p.ambush = rng.nextDouble() < ambushChance;
  p.pirateCount = p.ambush ? (legit ? (1 + rng.range(0, 2)) : (2 + rng.range(0, 3))) : 0;

  // Some illegitimate calls are *only* an ambush, but occasionally the trap has a "victim".
  p.hasVictim = legit || (rng.nextDouble() < 0.25);

  // Risk shaping.
  p.risk = 0.20;
  if (!legit) p.risk += 0.25;
  if (p.ambush) p.risk += 0.45;
  if (p.pirateCount >= 3) p.risk += 0.12;
  p.risk += rng.range(-0.08, 0.12);
  p.risk = clamp01(p.risk);

  if (!p.hasVictim) {
    p.scenario = DistressScenario::Ambush;
    p.needCommodity = econ::CommodityId::Food;
    p.needUnits = 0.0;
    p.rewardCr = 0.0;
    p.repReward = 0.0;
    return p;
  }

  // Pick a request category.
  // (Keep early-game friendly: avoid contraband / specialty items.)
  const double r = rng.nextDouble();
  if (r < 0.36) {
    p.scenario = DistressScenario::Supplies;
  } else if (r < 0.60) {
    p.scenario = DistressScenario::Fuel;
  } else if (r < 0.82) {
    p.scenario = DistressScenario::Medical;
  } else {
    p.scenario = DistressScenario::Mechanical;
  }

  // Choose requested commodity + units.
  switch (p.scenario) {
    case DistressScenario::Supplies: {
      p.needCommodity = (rng.nextDouble() < 0.65) ? econ::CommodityId::Food : econ::CommodityId::Water;
      p.needUnits = (double)(8 + rng.range(0, 15)); // 8..23
    } break;
    case DistressScenario::Fuel: {
      p.needCommodity = econ::CommodityId::Fuel;
      p.needUnits = (double)(5 + rng.range(0, 11)); // 5..16
    } break;
    case DistressScenario::Medical: {
      p.needCommodity = econ::CommodityId::Medicine;
      p.needUnits = (double)(3 + rng.range(0, 6)); // 3..9
      p.risk = clamp01(p.risk + 0.10);
    } break;
    case DistressScenario::Mechanical: {
      p.needCommodity = (rng.nextDouble() < 0.55) ? econ::CommodityId::Machinery : econ::CommodityId::Metals;
      p.needUnits = (double)(3 + rng.range(0, 8)); // 3..11
      p.risk = clamp01(p.risk + 0.06);
    } break;
    default: break;
  }

  // Reward curve: pay above market value; add risk premium for ambushy calls.
  const auto def = econ::commodityDef(p.needCommodity);
  const double goodsValue = std::max(0.0, p.needUnits) * std::max(0.0, def.basePrice);

  const double mult = 1.65 + 0.85 * p.risk + rng.range(-0.10, 0.20);
  const double base = 180.0 + 520.0 * p.risk + rng.range(0.0, 180.0);
  p.rewardCr = std::max(0.0, base + goodsValue * mult);

  // Keep payouts in a readable/prototype-friendly band.
  p.rewardCr = std::clamp(p.rewardCr, 160.0, 4800.0);

  // Rep: small positive bump for helping civilians.
  p.repReward = std::clamp(0.40 + 1.10 * p.risk, 0.25, 2.75);

  return p;
}

} // namespace stellar::sim
