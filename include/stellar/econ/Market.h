#pragma once

#include "stellar/econ/Economy.h"

namespace stellar::econ {

struct MarketQuote {
  double mid{0.0};
  double ask{0.0}; // station sells to player
  double bid{0.0}; // station buys from player
  double inventory{0.0};
};

struct TradeResult {
  bool ok{false};
  double creditsDelta{0.0}; // negative when spending
  double unitsDelta{0.0};   // positive when acquiring cargo
  const char* reason{nullptr};
};

MarketQuote quote(const StationEconomyState& state,
                  const StationEconomyModel& model,
                  CommodityId id,
                  double bidAskSpread = 0.10);

TradeResult buy(StationEconomyState& state,
                const StationEconomyModel& model,
                CommodityId id,
                double units,
                double& playerCredits,
                double bidAskSpread = 0.10,
                double stationFeeRate = 0.0);

TradeResult sell(StationEconomyState& state,
                 const StationEconomyModel& model,
                 CommodityId id,
                 double units,
                 double& playerCredits,
                 double bidAskSpread = 0.10,
                 double stationFeeRate = 0.0);

} // namespace stellar::econ
