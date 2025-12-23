#pragma once

#include "stellar/core/Random.h"
#include "stellar/econ/Commodity.h"

#include <array>
#include <vector>

namespace stellar::econ {

enum class StationType : core::u8 {
  Outpost = 0,
  Agricultural,
  Mining,
  Refinery,
  Industrial,
  Research,
  TradeHub,
  Shipyard,
  Count
};

struct PricePoint {
  double day{0.0};
  double price{0.0};
};

struct StationEconomyModel {
  StationType type{StationType::Outpost};

  // Units: units/day
  std::array<double, kCommodityCount> productionPerDay{};
  std::array<double, kCommodityCount> consumptionPerDay{};

  // Target inventory level
  std::array<double, kCommodityCount> desiredStock{};

  // Hard cap
  std::array<double, kCommodityCount> capacity{};

  // Price tuning
  double priceVolatility{0.9};   // how strongly inventory affects price
  double shockVolatility{0.02};  // random drift per update
};

struct StationEconomyState {
  double lastUpdateDay{0.0};
  double lastSampleDay{0.0};

  std::array<double, kCommodityCount> inventory{};
  std::array<std::vector<PricePoint>, kCommodityCount> history{};

  void clampToCapacity(const StationEconomyModel& model);
};

StationEconomyModel makeEconomyModel(StationType type, double industryBias /*-1..+1*/);

// Initialize a station's economy state deterministically.
StationEconomyState makeInitialState(const StationEconomyModel& model, core::SplitMix64& rng);

// Advance station economy to `timeDays` (in-place).
void updateEconomyTo(StationEconomyState& state,
                     const StationEconomyModel& model,
                     double timeDays,
                     core::SplitMix64& rng,
                     double sampleIntervalDays = 0.25);

// Compute a "mid" price for a commodity given current inventory.
double midPrice(const StationEconomyState& state, const StationEconomyModel& model, CommodityId id);

} // namespace stellar::econ
