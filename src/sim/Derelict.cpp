#include "stellar/sim/Derelict.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace stellar::sim {

static double clamp01(double x) { return std::clamp(x, 0.0, 1.0); }

static econ::CommodityId pickWeightedCommodity(core::SplitMix64& rng,
                                             const std::vector<econ::CommodityId>& items) {
  if (items.empty()) return econ::CommodityId::Machinery;

  // Weight by base price so higher-value salvage is (slightly) rarer.
  double total = 0.0;
  for (const auto cid : items) {
    const double w = std::max(1.0, econ::commodityDef(cid).basePrice);
    total += 1.0 / std::sqrt(w); // inverse-ish weight
  }

  const double r = rng.nextDouble() * total;
  double acc = 0.0;
  for (const auto cid : items) {
    const double w = std::max(1.0, econ::commodityDef(cid).basePrice);
    acc += 1.0 / std::sqrt(w);
    if (r <= acc) return cid;
  }
  return items.back();
}

DerelictPlan planDerelictEncounter(core::u64 universeSeed,
                                   SystemId systemId,
                                   core::u64 signalId,
                                   double timeDays,
                                   double piracy01,
                                   double security01,
                                   double contest01,
                                   bool missionSite,
                                   bool includeDayStamp) {
  DerelictPlan p{};

  // Keep generation stable across versions.
  core::u64 s = core::hashCombine(universeSeed, core::seedFromText("derelict_plan_v1"));
  s = core::hashCombine(s, (core::u64)systemId);
  s = core::hashCombine(s, signalId);
  s = core::hashCombine(s, (core::u64)(missionSite ? 1u : 0u));
  if (includeDayStamp) {
    const core::u64 day = (core::u64)std::max(0.0, std::floor(timeDays));
    s = core::hashCombine(s, day);
  }

  // Mix in the coarse security knobs so different systems feel distinct even if
  // the signal ids happen to align.
  auto q = [&](double x) -> core::u64 {
    // Quantize to avoid tiny floating drift affecting determinism.
    const long long v = (long long)std::llround(std::clamp(x, 0.0, 1.0) * 1000.0);
    return (core::u64)std::max(0ll, std::min(1000ll, v));
  };
  s = core::hashCombine(s, q(piracy01));
  s = core::hashCombine(s, q(security01));
  s = core::hashCombine(s, q(contest01));

  core::SplitMix64 rng(s);

  piracy01 = clamp01(piracy01);
  security01 = clamp01(security01);
  contest01 = clamp01(contest01);

  // --- Wreck class (size bucket) ---
  {
    // Mission sites skew larger so they feel "worth the trip".
    const double r = rng.nextDouble();
    if (missionSite) {
      if (r < 0.18) p.wreckClass = 0;
      else if (r < 0.70) p.wreckClass = 1;
      else p.wreckClass = 2;
    } else {
      if (r < 0.56) p.wreckClass = 0;
      else if (r < 0.88) p.wreckClass = 1;
      else p.wreckClass = 2;
    }
  }

  // --- Scenario + ambush selection ---
  // Ambush chance rises with piracy and contestedness, and is higher for mission sites.
  double ambushP = 0.10 + 0.50 * piracy01 + 0.18 * contest01 - 0.22 * security01;
  if (missionSite) ambushP += 0.15;
  if (p.wreckClass >= 2) ambushP += 0.06;
  ambushP = std::clamp(ambushP, 0.02, 0.78);

  p.ambush = (rng.nextDouble() < ambushP);

  if (p.ambush) {
    p.scenario = DerelictScenario::PirateTrap;
  } else {
    // Non-ambush derelicts still vary in flavor based on system conditions.
    const double r = rng.nextDouble();
    const double smugP = std::clamp(0.08 + 0.30 * piracy01 - 0.10 * security01, 0.0, 0.45);
    const double milP  = std::clamp(0.10 + 0.28 * security01 - 0.06 * piracy01, 0.0, 0.40);

    if (r < smugP) p.scenario = DerelictScenario::SmugglerCache;
    else if (r < smugP + milP) p.scenario = DerelictScenario::MilitaryHulk;
    else p.scenario = DerelictScenario::CivilianWreck;
  }

  // --- Loot table ---
  std::vector<econ::CommodityId> loot;
  loot.reserve(8);

  switch (p.scenario) {
    case DerelictScenario::CivilianWreck:
      loot = {econ::CommodityId::Metals,
              econ::CommodityId::Machinery,
              econ::CommodityId::Electronics,
              econ::CommodityId::Fuel,
              econ::CommodityId::Medicine};
      break;
    case DerelictScenario::SmugglerCache:
      loot = {econ::CommodityId::Luxury,
              econ::CommodityId::Stimulants,
              econ::CommodityId::Weapons,
              econ::CommodityId::Electronics};
      break;
    case DerelictScenario::MilitaryHulk:
      loot = {econ::CommodityId::Weapons,
              econ::CommodityId::Medicine,
              econ::CommodityId::Electronics,
              econ::CommodityId::Machinery};
      break;
    case DerelictScenario::PirateTrap:
      loot = {econ::CommodityId::Weapons,
              econ::CommodityId::Luxury,
              econ::CommodityId::Stimulants,
              econ::CommodityId::Electronics,
              econ::CommodityId::Machinery};
      break;
    default:
      loot = {econ::CommodityId::Machinery};
      break;
  }

  p.salvageCommodity = pickWeightedCommodity(rng, loot);

  // --- Salvage quantity + pod split ---
  // Units scale with wreckClass, and traps tend to be a bit juicier.
  {
    const int baseMin = (p.wreckClass == 0) ? 3 : (p.wreckClass == 1 ? 6 : 10);
    const int baseMax = (p.wreckClass == 0) ? 11 : (p.wreckClass == 1 ? 18 : 26);

    int units = baseMin + rng.range(0, std::max(0, baseMax - baseMin));
    if (p.scenario == DerelictScenario::PirateTrap) {
      units = (int)std::llround((double)units * (1.15 + rng.range(0.0, 0.55)));
    }
    // Keep salvage in a prototype-friendly band.
    units = std::clamp(units, 2, 36);
    p.salvageUnits = (double)units;
  }

  {
    int pods = 2 + rng.range(0, 2);
    if (p.wreckClass == 1) pods = 3 + rng.range(0, 2);
    if (p.wreckClass >= 2) pods = 4 + rng.range(0, 2);
    pods = std::clamp(pods, 2, 6);
    p.salvagePods = pods;
  }

  // --- Data core ---
  {
    double dataP = 0.18 + 0.10 * (double)p.wreckClass;
    if (p.scenario == DerelictScenario::MilitaryHulk) dataP += 0.18;
    if (p.scenario == DerelictScenario::SmugglerCache) dataP += 0.10;
    if (missionSite) dataP += 0.07;
    dataP += 0.10 * security01;
    dataP = std::clamp(dataP, 0.05, 0.80);

    p.hasDataCore = (rng.nextDouble() < dataP);
    p.dataUnits = p.hasDataCore ? 1.0 : 0.0;
    p.dataCommodity = econ::CommodityId::Electronics;
  }

  // --- Pirate pack shaping ---
  if (p.ambush) {
    int c = 2;
    if (missionSite) c += 1;
    if (piracy01 > 0.65) c += 1;
    if (p.wreckClass >= 2 && rng.nextDouble() < 0.60) c += 1;
    if (security01 < 0.35 && rng.nextDouble() < 0.45) c += 1;

    c = std::clamp(c, 2, 5);
    p.pirateCount = c;
  } else {
    p.pirateCount = 0;
  }

  // --- Risk score ---
  {
    double r = 0.12;
    r += 0.45 * piracy01;
    r += 0.18 * contest01;
    r -= 0.28 * security01;
    r += 0.06 * (double)p.wreckClass;
    if (p.scenario == DerelictScenario::SmugglerCache) r += 0.05;
    if (p.scenario == DerelictScenario::MilitaryHulk) r += 0.08;
    if (p.ambush) r += 0.38;
    if (p.pirateCount >= 4) r += 0.08;
    r += rng.range(-0.06, 0.10);
    p.risk = clamp01(r);
  }

  // Edge safety.
  if (p.salvageUnits < 1e-6) {
    p.hasSalvage = false;
    p.salvagePods = 0;
  }
  if (!p.hasDataCore) {
    p.dataUnits = 0.0;
  }

  return p;
}

} // namespace stellar::sim
