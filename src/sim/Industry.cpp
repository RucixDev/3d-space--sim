#include "stellar/sim/Industry.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cctype>
#include <cmath>

namespace stellar::sim {

static constexpr core::u32 stationBit(econ::StationType t) {
  return 1u << static_cast<core::u32>(t);
}

const IndustryRecipeTable& industryRecipeTable() {
  // NOTE: Recipes are intentionally lossy in *value* terms (see input ratios)
  // to avoid trivial "buy inputs, process, sell outputs" loops at the same
  // station. They primarily exist as a mass/value-density transform and as a
  // way to smooth local shortages.
  static const IndustryRecipeTable kTable = {
      // id, code, name, desc, stationMask, inA, inAUnits, inB, inBUnits, out, outUnits, timeDays, feeCr
      IndustryRecipeDef{IndustryRecipeId::SmeltOre,
                        "SMELT_ORE",
                        "Smelt ore",
                        "Convert raw ore into refined metals (lossy, reduces mass).",
                        stationBit(econ::StationType::Refinery),
                        econ::CommodityId::Ore,
                        4.0,
                        econ::CommodityId::Food,
                        0.0,
                        econ::CommodityId::Metals,
                        1.0,
                        0.010,
                        18.0},

      IndustryRecipeDef{IndustryRecipeId::MachineParts,
                        "MACHINE_PARTS",
                        "Machine parts",
                        "Machine metals into general-purpose machinery components.",
                        stationBit(econ::StationType::Industrial) | stationBit(econ::StationType::Shipyard),
                        econ::CommodityId::Metals,
                        3.0,
                        econ::CommodityId::Electronics,
                        0.5,
                        econ::CommodityId::Machinery,
                        1.0,
                        0.020,
                        32.0},

      IndustryRecipeDef{IndustryRecipeId::CircuitFab,
                        "CIRCUIT_FAB",
                        "Circuit fab",
                        "Fabricate electronics from machinery and refined metals.",
                        stationBit(econ::StationType::Research) | stationBit(econ::StationType::Industrial),
                        econ::CommodityId::Machinery,
                        2.0,
                        econ::CommodityId::Metals,
                        2.0,
                        econ::CommodityId::Electronics,
                        1.0,
                        0.030,
                        55.0},

      IndustryRecipeDef{IndustryRecipeId::AssembleWeapons,
                        "ASSEMBLE_WEAPONS",
                        "Assemble weapons",
                        "Assemble weapon crates from machinery, metals, and electronics.",
                        stationBit(econ::StationType::Industrial) | stationBit(econ::StationType::Shipyard),
                        econ::CommodityId::Machinery,
                        3.0,
                        econ::CommodityId::Metals,
                        5.0,
                        econ::CommodityId::Weapons,
                        1.0,
                        0.045,
                        110.0},
  };
  return kTable;
}

const IndustryRecipeDef& industryRecipeDef(IndustryRecipeId id) {
  const std::size_t idx = static_cast<std::size_t>(id);
  const auto& tbl = industryRecipeTable();
  if (idx >= tbl.size()) return tbl[0];
  return tbl[idx];
}

const IndustryRecipeDef* findIndustryRecipe(IndustryRecipeId id) {
  const std::size_t idx = static_cast<std::size_t>(id);
  const auto& tbl = industryRecipeTable();
  if (idx >= tbl.size()) return nullptr;
  return &tbl[idx];
}

static bool ieq(std::string_view a, std::string_view b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const unsigned char ca = static_cast<unsigned char>(a[i]);
    const unsigned char cb = static_cast<unsigned char>(b[i]);
    if (std::toupper(ca) != std::toupper(cb)) return false;
  }
  return true;
}

const IndustryRecipeDef* findIndustryRecipeByCode(std::string_view code) {
  for (const auto& r : industryRecipeTable()) {
    if (ieq(code, r.code)) return &r;
  }
  return nullptr;
}

bool stationSupportsIndustry(econ::StationType stationType, const IndustryRecipeDef& recipe) {
  const core::u32 bit = stationBit(stationType);
  return (recipe.stationTypeMask & bit) != 0u;
}

bool stationSupportsIndustry(econ::StationType stationType, IndustryRecipeId recipeId) {
  return stationSupportsIndustry(stationType, industryRecipeDef(recipeId));
}

std::vector<const IndustryRecipeDef*> availableIndustryRecipes(econ::StationType stationType) {
  std::vector<const IndustryRecipeDef*> out;
  for (const auto& r : industryRecipeTable()) {
    if (stationSupportsIndustry(stationType, r)) out.push_back(&r);
  }
  return out;
}

IndustryStationModifiers stationIndustryModifiers(StationId stationId) {
  // Deterministic per-station "quality" modifiers.
  // Keep these ranges fairly narrow so economy balance doesn't swing wildly.
  const core::u64 seed = core::hashCombine((core::u64)stationId, core::fnv1a64("industry_mods_v1"));
  core::SplitMix64 rng(seed);

  IndustryStationModifiers m;
  m.yieldMul = rng.range(0.92, 1.08);
  m.speedMul = rng.range(0.90, 1.12);
  m.feeMul = rng.range(0.92, 1.15);
  return m;
}

static double clampRep(double rep) { return std::clamp(rep, -100.0, 100.0); }
static double repNorm(double rep) { return clampRep(rep) / 100.0; }

IndustryQuote quoteIndustryOrder(const IndustryRecipeDef& recipe,
                                 StationId stationId,
                                 econ::StationType stationType,
                                 double batches,
                                 double effectiveStationFeeRate,
                                 double rep) {
  IndustryQuote q{};
  q.recipe = recipe.id;
  q.stationId = stationId;
  q.stationType = stationType;

  if (!std::isfinite(batches)) batches = 0.0;
  q.batches = std::max(0.0, batches);

  // If the station can't run it, return a zero quote.
  if (!stationSupportsIndustry(stationType, recipe)) {
    q.batches = 0.0;
    return q;
  }

  q.inputA = recipe.inputA;
  q.inputAUnits = recipe.inputAUnits * q.batches;
  q.inputB = recipe.inputB;
  q.inputBUnits = recipe.inputBUnits * q.batches;
  q.output = recipe.output;

  q.mods = stationIndustryModifiers(stationId);

  // Small rep perks (kept intentionally subtle).
  // - Friendly factions run the line slightly faster and with a bit less waste.
  const double rn = repNorm(rep);
  q.mods.yieldMul *= (1.0 + 0.06 * rn);
  q.mods.speedMul *= (1.0 + 0.04 * rn);

  // Output + timing.
  q.outputUnits = recipe.outputUnits * q.batches * q.mods.yieldMul;

  // Higher speedMul => shorter time.
  const double speed = std::max(0.25, q.mods.speedMul);
  q.timeDays = (recipe.baseTimeDays * q.batches) / speed;

  // Service fee: baseline * batches, with station fee multiplier.
  if (!std::isfinite(effectiveStationFeeRate)) effectiveStationFeeRate = 0.0;
  effectiveStationFeeRate = std::clamp(effectiveStationFeeRate, 0.0, 1.0);
  q.serviceFeeCr = recipe.baseServiceFeeCr * q.batches * q.mods.feeMul * (1.0 + effectiveStationFeeRate);

  // Defensive clamps.
  if (!std::isfinite(q.outputUnits)) q.outputUnits = 0.0;
  if (!std::isfinite(q.timeDays)) q.timeDays = 0.0;
  if (!std::isfinite(q.serviceFeeCr)) q.serviceFeeCr = 0.0;
  q.outputUnits = std::max(0.0, q.outputUnits);
  q.timeDays = std::max(0.0, q.timeDays);
  q.serviceFeeCr = std::max(0.0, q.serviceFeeCr);

  return q;
}

} // namespace stellar::sim
