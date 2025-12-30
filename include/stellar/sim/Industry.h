#pragma once

#include "stellar/core/Random.h"
#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/econ/Economy.h"
#include "stellar/sim/Celestial.h"

#include <array>
#include <cstddef>
#include <string_view>
#include <vector>

namespace stellar::sim {

// Lightweight "industry" / processing system.
//
// This is intentionally small and deterministic:
//  - recipes are stable & persisted by explicit numeric ids
//  - stations provide mild deterministic modifiers (speed/yield/fees)
//  - gameplay can convert commodities over time via "orders"
//
// The game client can use this to expose station services like:
//   Ore -> Metals (Refinery)
//   Metals+Machinery -> Weapons (Industrial)
// without hard-coding recipe math into UI code.

// IMPORTANT: IndustryRecipeId is persisted to save files.
// Always use explicit, stable numeric values here.
// Do not reorder existing members.
enum class IndustryRecipeId : core::u8 {
  SmeltOre        = 0,
  MachineParts    = 1,
  CircuitFab      = 2,
  AssembleWeapons = 3,
  Count           = 4,
};

struct IndustryRecipeDef {
  IndustryRecipeId id{IndustryRecipeId::SmeltOre};

  // Stable text code (useful for tooling/CLI flags).
  const char* code{""};
  const char* name{""};
  const char* desc{""};

  // Bitmask of econ::StationType values that can run this recipe.
  // Bit i corresponds to static_cast<core::u8>(StationType).
  core::u32 stationTypeMask{0};

  // Up to two inputs (B can be omitted by setting inputBUnits=0).
  econ::CommodityId inputA{econ::CommodityId::Food};
  double inputAUnits{0.0};
  econ::CommodityId inputB{econ::CommodityId::Food};
  double inputBUnits{0.0};

  // Single output commodity.
  econ::CommodityId output{econ::CommodityId::Food};
  double outputUnits{0.0};

  // Baseline processing time & service fee per batch.
  double baseTimeDays{0.0};
  double baseServiceFeeCr{0.0};
};

struct IndustryStationModifiers {
  // Multiply output by yieldMul.
  double yieldMul{1.0};
  // Multiply speed by speedMul (higher => faster).
  double speedMul{1.0};
  // Multiply service fees by feeMul.
  double feeMul{1.0};
};

struct IndustryQuote {
  IndustryRecipeId recipe{IndustryRecipeId::SmeltOre};
  StationId stationId{0};
  econ::StationType stationType{econ::StationType::Outpost};

  double batches{0.0};

  econ::CommodityId inputA{econ::CommodityId::Food};
  double inputAUnits{0.0};
  econ::CommodityId inputB{econ::CommodityId::Food};
  double inputBUnits{0.0};
  econ::CommodityId output{econ::CommodityId::Food};
  double outputUnits{0.0};

  double serviceFeeCr{0.0};
  double timeDays{0.0};

  IndustryStationModifiers mods{};
};

// Persistable, lightweight "industry order".
// The order stores a fully specified contract (inputs/outputs), so it remains valid
// even if future patches rebalance recipe defaults.
struct IndustryOrder {
  core::u64 id{0};
  IndustryRecipeId recipe{IndustryRecipeId::SmeltOre};
  StationId stationId{0};

  econ::CommodityId inputA{econ::CommodityId::Food};
  double inputAUnits{0.0};
  econ::CommodityId inputB{econ::CommodityId::Food};
  double inputBUnits{0.0};
  econ::CommodityId output{econ::CommodityId::Food};
  double outputUnits{0.0};

  double submittedDay{0.0};
  double readyDay{0.0};
  bool claimed{false};
};

constexpr std::size_t kIndustryRecipeCount = static_cast<std::size_t>(IndustryRecipeId::Count);
using IndustryRecipeTable = std::array<IndustryRecipeDef, kIndustryRecipeCount>;

const IndustryRecipeTable& industryRecipeTable();
const IndustryRecipeDef& industryRecipeDef(IndustryRecipeId id);
const IndustryRecipeDef* findIndustryRecipe(IndustryRecipeId id);

const IndustryRecipeDef* findIndustryRecipeByCode(std::string_view code);

bool stationSupportsIndustry(econ::StationType stationType, const IndustryRecipeDef& recipe);
bool stationSupportsIndustry(econ::StationType stationType, IndustryRecipeId recipeId);

std::vector<const IndustryRecipeDef*> availableIndustryRecipes(econ::StationType stationType);

// Deterministic station-specific modifiers (derived from station id).
IndustryStationModifiers stationIndustryModifiers(StationId stationId);

// Quote a processing order for a station.
//
// Parameters:
//  - batches: number of recipe batches (>=0)
//  - effectiveStationFeeRate: already adjusted for reputation, if desired (0..1)
//  - rep: player rep [-100,+100]; used to give small yield/speed perks
IndustryQuote quoteIndustryOrder(const IndustryRecipeDef& recipe,
                                 StationId stationId,
                                 econ::StationType stationType,
                                 double batches,
                                 double effectiveStationFeeRate,
                                 double rep);

} // namespace stellar::sim
