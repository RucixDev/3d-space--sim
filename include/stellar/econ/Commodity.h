#pragma once

#include "stellar/core/Types.h"

#include <array>
#include <string_view>

namespace stellar::econ {

enum class CommodityId : core::u16 {
  Food = 0,
  Water,
  Ore,
  Metals,
  Fuel,
  Machinery,
  Medicine,
  Electronics,
  Luxury,

  // Higher-value / specialty goods.
  Weapons,
  Stimulants,

  Count
};

constexpr std::size_t kCommodityCount = static_cast<std::size_t>(CommodityId::Count);

struct CommodityDef {
  CommodityId id{};
  const char* code{};     // short symbol for save files / UI
  const char* name{};     // display name
  double basePrice{};     // "credits" per unit (mid price baseline)
  double massKg{};        // kg per unit (for cargo capacity later)
};

const std::array<CommodityDef, kCommodityCount>& commodityTable();
const CommodityDef& commodityDef(CommodityId id);
std::string_view commodityName(CommodityId id);
std::string_view commodityCode(CommodityId id);

// Parse helpers (case-insensitive).
//
// These are primarily intended for headless tooling / CLI inputs.
//
// Examples:
//  - tryParseCommodityCode("FOOD") => Food
//  - tryParseCommodity("water")    => Water
//
// Returns true on success and writes the parsed id to `out`.
bool tryParseCommodityCode(std::string_view code, CommodityId& out);
bool tryParseCommodity(std::string_view codeOrName, CommodityId& out);

} // namespace stellar::econ
