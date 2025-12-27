#include "stellar/econ/Commodity.h"

#include <cctype>
#include <string>

namespace stellar::econ {

static const std::array<CommodityDef, kCommodityCount> kTable = {{
  {CommodityId::Food,        "FOOD", "Food",        12.0,  1.0},
  {CommodityId::Water,       "H2O",  "Water",        6.0,  1.0},
  {CommodityId::Ore,         "ORE",  "Ore",          9.0,  4.0},
  {CommodityId::Metals,      "MET",  "Metals",      22.0,  3.0},
  {CommodityId::Fuel,        "FUEL", "Fuel",        35.0,  1.0},
  {CommodityId::Machinery,   "MACH", "Machinery",   55.0,  2.0},
  {CommodityId::Medicine,    "MED",  "Medicine",    75.0,  0.5},
  {CommodityId::Electronics, "ELEC", "Electronics", 85.0,  0.5},
  {CommodityId::Luxury,      "LUX",  "Luxury",     120.0,  0.2},

  {CommodityId::Weapons,     "ARMS", "Weapons",    150.0,  1.6},
  {CommodityId::Stimulants,  "STIM", "Stimulants", 165.0,  0.2},
}};

const std::array<CommodityDef, kCommodityCount>& commodityTable() { return kTable; }

const CommodityDef& commodityDef(CommodityId id) {
  return kTable[static_cast<std::size_t>(id)];
}

std::string_view commodityName(CommodityId id) {
  return commodityDef(id).name;
}

std::string_view commodityCode(CommodityId id) {
  return commodityDef(id).code;
}

static std::string normalizeToken(std::string_view s, bool keepDigits = true) {
  // Uppercase, strip whitespace. Optionally strip digits for name matching.
  std::string out;
  out.reserve(s.size());
  for (const unsigned char uc : s) {
    if (std::isspace(uc)) continue;
    if (!keepDigits && std::isdigit(uc)) continue;
    out.push_back((char)std::toupper(uc));
  }
  return out;
}

static std::string normalizeName(std::string_view s) {
  // Uppercase, keep alphanumerics only.
  std::string out;
  out.reserve(s.size());
  for (const unsigned char uc : s) {
    if (std::isalnum(uc)) out.push_back((char)std::toupper(uc));
  }
  return out;
}

bool tryParseCommodityCode(std::string_view code, CommodityId& out) {
  const std::string tok = normalizeToken(code);
  if (tok.empty()) return false;

  for (const auto& def : commodityTable()) {
    if (tok == def.code) {
      out = def.id;
      return true;
    }
  }
  return false;
}

bool tryParseCommodity(std::string_view codeOrName, CommodityId& out) {
  if (tryParseCommodityCode(codeOrName, out)) return true;

  const std::string tok = normalizeName(codeOrName);
  if (tok.empty()) return false;

  for (const auto& def : commodityTable()) {
    if (tok == normalizeName(def.name)) {
      out = def.id;
      return true;
    }
  }

  return false;
}

} // namespace stellar::econ
