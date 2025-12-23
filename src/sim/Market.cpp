#include "stellar/sim/Market.h"

#include "stellar/core/Random.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace stellar::sim {
namespace {

constexpr std::array<std::string_view, kCommodityCount> kKeys = {
  "food",
  "water",
  "ore",
  "metals",
  "electronics",
  "medicine",
  "machinery",
  "luxury",
  "fuel",
};

constexpr std::array<std::string_view, kCommodityCount> kNames = {
  "Food",
  "Water",
  "Ore",
  "Metals",
  "Electronics",
  "Medicine",
  "Machinery",
  "Luxury Goods",
  "Fuel",
};

constexpr std::array<double, kCommodityCount> kBasePrices = {
  12.0,  // Food
  8.0,   // Water
  25.0,  // Ore
  55.0,  // Metals
  130.0, // Electronics
  160.0, // Medicine
  110.0, // Machinery
  220.0, // Luxury
  45.0,  // Fuel
};

double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

struct SystemFactors {
  double agriculture = 1.0;
  double mining = 1.0;
  double refining = 1.0;
  double highTech = 1.0;
  double medicine = 1.0;
  double fuel = 1.0;
  double luxury = 1.0;
};

SystemFactors factorsFromSystem(const StarSystem& sys) {
  int rocky = 0;
  int ice = 0;
  int gas = 0;
  int belts = 0;
  int temperateRocky = 0;

  for (const auto& p : sys.planets) {
    switch (p.type) {
      case PlanetType::Rocky:
        ++rocky;
        if (p.equilibriumTempK > 240.0 && p.equilibriumTempK < 330.0) {
          ++temperateRocky;
        }
        break;
      case PlanetType::Ice:
        ++ice;
        break;
      case PlanetType::GasGiant:
      case PlanetType::IceGiant:
        ++gas;
        break;
      case PlanetType::AsteroidBelt:
        ++belts;
        break;
    }
  }

  SystemFactors f;

  // Very simple heuristics: more bodies -> more activity.
  const double bodyCount = static_cast<double>(sys.planets.size());
  const double activity = 0.85 + 0.03 * bodyCount;

  f.agriculture = activity * (0.80 + 0.25 * temperateRocky);
  f.mining = activity * (0.90 + 0.15 * rocky + 0.30 * belts);
  f.refining = activity * (0.90 + 0.10 * rocky);
  f.fuel = activity * (0.95 + 0.08 * gas);

  // High-tech is not really a physical factor; keep it mild.
  f.highTech = activity;
  f.medicine = activity;
  f.luxury = activity;

  // Avoid extreme factors.
  f.agriculture = std::clamp(f.agriculture, 0.75, 2.25);
  f.mining = std::clamp(f.mining, 0.75, 2.75);
  f.refining = std::clamp(f.refining, 0.75, 2.25);
  f.fuel = std::clamp(f.fuel, 0.75, 2.25);
  return f;
}

double factionPriceMultiplier(const Faction& fac, Commodity c) {
  const double tech = clamp01(fac.techLevel);
  const double wealth = clamp01(fac.wealth);
  const double law = clamp01(fac.lawfulness);

  // Wealthier factions have higher prices (higher purchasing power / higher costs).
  const double wealthMult = 0.85 + 0.55 * wealth;
  // Higher tech reduces the price of advanced goods (efficiency), but may increase demand.
  const double techMult = 1.05 - 0.25 * tech;
  // Higher lawfulness adds a small "overhead".
  const double lawMult = 0.95 + 0.20 * law;

  switch (c) {
    case Commodity::Electronics:
    case Commodity::Machinery:
    case Commodity::Medicine:
      return wealthMult * techMult * lawMult;
    case Commodity::Luxury:
      return (0.95 + 0.70 * wealth) * lawMult;
    default:
      return wealthMult * lawMult;
  }
}

double systemSupplyMultiplier(const SystemFactors& f, Commodity c) {
  switch (c) {
    case Commodity::Food: return f.agriculture;
    case Commodity::Water: return f.agriculture;
    case Commodity::Ore: return f.mining;
    case Commodity::Metals: return f.refining;
    case Commodity::Fuel: return f.fuel;
    case Commodity::Electronics: return f.highTech;
    case Commodity::Machinery: return f.highTech;
    case Commodity::Medicine: return f.medicine;
    case Commodity::Luxury: return f.luxury;
    case Commodity::Count: break;
  }
  return 1.0;
}

double systemPriceMultiplier(const SystemFactors& f, Commodity c) {
  // If supply is high, price tends to be lower.
  const double s = systemSupplyMultiplier(f, c);
  return std::clamp(1.30 - 0.35 * (s - 1.0), 0.65, 1.65);
}

}

std::string_view commodityKey(Commodity c) {
  return kKeys[static_cast<std::size_t>(c)];
}

std::string_view commodityName(Commodity c) {
  return kNames[static_cast<std::size_t>(c)];
}

Market MarketGenerator::generate(const StarSystem& sys, const Faction& faction, std::int64_t day) const {
  Market m;
  m.offers.reserve(kCommodityCount);

  const auto fSys = factorsFromSystem(sys);

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    const auto c = static_cast<Commodity>(i);

    // Stable seed per (system, day, commodity)
    const auto s0 = stellar::core::deriveSeed(m_seed, "market");
    const auto s1 = stellar::core::deriveSeed(s0, sys.id);
    const auto s2 = stellar::core::deriveSeed(s1, static_cast<stellar::core::u64>(day));
    const auto s3 = stellar::core::deriveSeed(s2, static_cast<stellar::core::u64>(i));
    stellar::core::SplitMix64 rng(s3);

    const double base = kBasePrices[i];

    const double sysPriceMult = systemPriceMultiplier(fSys, c);
    const double facPriceMult = factionPriceMultiplier(faction, c);

    // Time noise: small daily volatility.
    const double noise = rng.uniform(-0.15, 0.15);
    const double price = std::max(0.01, base * sysPriceMult * facPriceMult * (1.0 + noise));

    // Supply scales with system productivity.
    const double supplyMult = systemSupplyMultiplier(fSys, c);
    const int supply = std::max(0, static_cast<int>(std::round(rng.uniform(60.0, 2200.0) * supplyMult)));

    m.offers.push_back(MarketOffer{
      .commodity = c,
      .supply = supply,
      .price = price,
    });
  }

  return m;
}

} // namespace stellar::sim
