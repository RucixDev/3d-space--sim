#include "stellar/econ/Market.h"

#include "stellar/econ/Commodity.h"
#include "stellar/econ/Economy.h"

#include <algorithm>

namespace stellar::econ {

static constexpr std::size_t idx(CommodityId id) { return static_cast<std::size_t>(id); }

double takeInventory(StationEconomyState& state,
                     const StationEconomyModel& model,
                     CommodityId id,
                     double units) {
  (void)model;
  if (units <= 0.0) return 0.0;

  double& inv = state.inventory[idx(id)];
  inv = std::max(0.0, inv);

  const double taken = std::min(inv, units);
  inv = std::max(0.0, inv - taken);
  return taken;
}

double addInventory(StationEconomyState& state,
                    const StationEconomyModel& model,
                    CommodityId id,
                    double units) {
  if (units <= 0.0) return 0.0;

  const double cap = std::max(0.0, model.capacity[idx(id)]);
  double& inv = state.inventory[idx(id)];
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
  const double half = bidAskSpread * 0.5;
  q.ask = q.mid * (1.0 + half);
  q.bid = q.mid * (1.0 - half);
  q.inventory = state.inventory[idx(id)];
  return q;
}

TradeResult buy(StationEconomyState& state,
                const StationEconomyModel& model,
                CommodityId id,
                double units,
                double& playerCredits,
                double bidAskSpread,
                double stationFeeRate) {
  if (units <= 0.0) return {false, 0.0, 0.0, "units<=0"};
  const auto q = quote(state, model, id, bidAskSpread);

  if (q.inventory + 1e-9 < units) return {false, 0.0, 0.0, "station out of stock"};

  const double total = q.ask * units * (1.0 + std::max(0.0, stationFeeRate));
  if (playerCredits + 1e-9 < total) return {false, 0.0, 0.0, "insufficient credits"};

  playerCredits -= total;
  state.inventory[idx(id)] = std::max(0.0, state.inventory[idx(id)] - units);

  return {true, -total, units, nullptr};
}

TradeResult sell(StationEconomyState& state,
                 const StationEconomyModel& model,
                 CommodityId id,
                 double units,
                 double& playerCredits,
                 double bidAskSpread,
                 double stationFeeRate) {
  if (units <= 0.0) return {false, 0.0, 0.0, "units<=0"};

  const auto q = quote(state, model, id, bidAskSpread);

  const double cap = model.capacity[idx(id)];
  const double cur = state.inventory[idx(id)];
  if (cur + units > cap + 1e-6) return {false, 0.0, 0.0, "station storage full"};

  const double payout = q.bid * units * (1.0 - std::max(0.0, stationFeeRate));
  playerCredits += payout;
  state.inventory[idx(id)] = std::min(cap, cur + units);

  return {true, payout, -units, nullptr};
}

} // namespace stellar::econ
