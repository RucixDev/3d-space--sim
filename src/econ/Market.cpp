#include "stellar/econ/Market.h"

#include "stellar/econ/Commodity.h"
#include "stellar/econ/Economy.h"

#include <algorithm>
#include <cmath>

namespace stellar::econ {

static constexpr std::size_t idx(CommodityId id) { return static_cast<std::size_t>(id); }

double takeInventory(StationEconomyState& state,
                     const StationEconomyModel& model,
                     CommodityId id,
                     double units) {
  if (units <= 0.0 || !std::isfinite(units)) return 0.0;

  const double cap = std::max(0.0, model.capacity[idx(id)]);
  double& inv = state.inventory[idx(id)];
  if (!std::isfinite(inv)) inv = 0.0;
  inv = std::clamp(inv, 0.0, cap);

  const double taken = std::min(inv, units);
  inv = std::clamp(inv - taken, 0.0, cap);
  return taken;
}

double addInventory(StationEconomyState& state,
                    const StationEconomyModel& model,
                    CommodityId id,
                    double units) {
  if (units <= 0.0 || !std::isfinite(units)) return 0.0;

  const double cap = std::max(0.0, model.capacity[idx(id)]);
  double& inv = state.inventory[idx(id)];
  if (!std::isfinite(inv)) inv = 0.0;
  inv = std::max(0.0, inv);

  const double space = std::max(0.0, cap - inv);
  const double added = std::min(space, units);
  inv = std::min(cap, inv + added);
  return added;
}

MarketQuote quote(const StationEconomyState& state,
                  const StationEconomyModel& model,
                  CommodityId id,
                  double bidAskSpread) {
  MarketQuote q{};
  q.mid = midPrice(state, model, id);

  // Defensive clamps: tooling/gameplay callers sometimes pass through user input.
  if (!std::isfinite(bidAskSpread)) bidAskSpread = 0.0;
  bidAskSpread = std::clamp(bidAskSpread, 0.0, 1.0);

  const double half = bidAskSpread * 0.5;
  q.ask = q.mid * (1.0 + half);
  q.bid = q.mid * (1.0 - half);

  const double cap = std::max(0.0, model.capacity[idx(id)]);
  double inv = state.inventory[idx(id)];
  if (!std::isfinite(inv)) inv = 0.0;
  q.inventory = std::clamp(inv, 0.0, cap);
  return q;
}

TradeResult buy(StationEconomyState& state,
                const StationEconomyModel& model,
                CommodityId id,
                double units,
                double& playerCredits,
                double bidAskSpread,
                double stationFeeRate) {
  if (units <= 0.0 || !std::isfinite(units)) return {false, 0.0, 0.0, "units<=0"};
  if (!std::isfinite(playerCredits)) return {false, 0.0, 0.0, "credits not finite"};

  // Clamp fee/spread for safety.
  if (!std::isfinite(stationFeeRate)) stationFeeRate = 0.0;
  stationFeeRate = std::clamp(stationFeeRate, 0.0, 1.0);

  const double cap = std::max(0.0, model.capacity[idx(id)]);
  double& inv = state.inventory[idx(id)];
  if (!std::isfinite(inv)) inv = 0.0;
  inv = std::clamp(inv, 0.0, cap);
  const auto q = quote(state, model, id, bidAskSpread);

  if (q.inventory + 1e-9 < units) return {false, 0.0, 0.0, "station out of stock"};

  const double total = q.ask * units * (1.0 + stationFeeRate);
  if (playerCredits + 1e-9 < total) return {false, 0.0, 0.0, "insufficient credits"};

  playerCredits -= total;
  inv = std::clamp(inv - units, 0.0, cap);

  return {true, -total, units, nullptr};
}

TradeResult sell(StationEconomyState& state,
                 const StationEconomyModel& model,
                 CommodityId id,
                 double units,
                 double& playerCredits,
                 double bidAskSpread,
                 double stationFeeRate) {
  if (units <= 0.0 || !std::isfinite(units)) return {false, 0.0, 0.0, "units<=0"};
  if (!std::isfinite(playerCredits)) return {false, 0.0, 0.0, "credits not finite"};

  if (!std::isfinite(stationFeeRate)) stationFeeRate = 0.0;
  stationFeeRate = std::clamp(stationFeeRate, 0.0, 1.0);

  const auto q = quote(state, model, id, bidAskSpread);

  const double cap = std::max(0.0, model.capacity[idx(id)]);
  double cur = state.inventory[idx(id)];
  if (!std::isfinite(cur)) cur = 0.0;
  cur = std::clamp(cur, 0.0, cap);
  if (cur + units > cap + 1e-6) return {false, 0.0, 0.0, "station storage full"};

  const double payout = q.bid * units * (1.0 - stationFeeRate);
  playerCredits += payout;
  state.inventory[idx(id)] = std::clamp(cur + units, 0.0, cap);

  return {true, payout, -units, nullptr};
}

} // namespace stellar::econ
