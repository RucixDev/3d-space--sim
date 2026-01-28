#pragma once

#include "stellar/core/Types.h"

#include <array>
#include <string_view>

namespace stellar::econ {

// Small helper wrapper for const char* literals used throughout economy defs.
//
// Motivation:
//  - some UI code historically used `.name.c_str()`
//  - other call sites treat names as raw C strings (ImGui, printf-style)
//  - other call sites want a `std::string_view`
//
// This keeps the underlying storage as a plain const char* while providing a
// safe `.c_str()` and implicit `std::string_view` / `const char*` views.
struct CStrView {
  const char* ptr{nullptr};

  constexpr CStrView() = default;
  constexpr CStrView(const char* s) : ptr(s) {}

  constexpr const char* c_str() const { return ptr ? ptr : ""; }

  constexpr operator const char*() const { return c_str(); }
  constexpr operator std::string_view() const { return std::string_view{c_str()}; }
};

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
  CStrView name{};        // display name
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
