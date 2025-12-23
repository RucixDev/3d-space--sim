#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/Celestial.h"
#include "stellar/sim/Faction.h"

#include <array>
#include <cstddef>
#include <optional>
#include <string_view>
#include <vector>

namespace stellar::sim {

// A small set of trade goods to support early gameplay.
enum class Commodity : stellar::core::u8 {
  Food,
  Water,
  Ore,
  Metals,
  Electronics,
  Medicine,
  Machinery,
  Luxury,
  Fuel,
  Count
};

constexpr std::size_t kCommodityCount = static_cast<std::size_t>(Commodity::Count);

std::string_view commodityKey(Commodity c);
std::string_view commodityName(Commodity c);

struct MarketOffer {
  Commodity commodity = Commodity::Food;
  int supply = 0;
  double price = 0.0; // credits per unit
};

struct Market {
  std::vector<MarketOffer> offers;

  const MarketOffer* find(Commodity c) const {
    for (const auto& o : offers) {
      if (o.commodity == c) return &o;
    }
    return nullptr;
  }
};

// Simple cargo hold for the player.
struct CargoHold {
  std::array<int, kCommodityCount> units{};

  int& operator[](Commodity c) { return units[static_cast<std::size_t>(c)]; }
  int operator[](Commodity c) const { return units[static_cast<std::size_t>(c)]; }
};

// Deterministic market generator.
// Markets are derived from universeSeed + system id + (integer) day.
class MarketGenerator {
public:
  explicit MarketGenerator(stellar::core::u64 universeSeed) : m_seed(universeSeed) {}

  Market generate(const StarSystem& sys, const Faction& faction, std::int64_t day) const;

private:
  stellar::core::u64 m_seed = 0;
};

} // namespace stellar::sim
