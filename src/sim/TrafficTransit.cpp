#include "stellar/sim/TrafficTransit.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/econ/Market.h"
#include "stellar/sim/Traffic.h"
#include "stellar/sim/TrafficLedger.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <vector>

namespace stellar::sim {
namespace {

static constexpr std::size_t idx(econ::CommodityId id) { return static_cast<std::size_t>(id); }

struct BestCommodity {
  econ::CommodityId id{econ::CommodityId::Food};
  double score{0.0};
};

BestCommodity pickCommodity(const Station& src,
                           const econ::StationEconomyState& srcState,
                           const Station& dst,
                           const econ::StationEconomyState& dstState,
                           core::SplitMix64& rng) {
  std::array<BestCommodity, 3> top{};
  for (auto& t : top) { t.score = 0.0; }

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const econ::CommodityId cid = static_cast<econ::CommodityId>(i);

    const double prodS = src.economyModel.productionPerDay[i];
    const double consS = src.economyModel.consumptionPerDay[i];
    const double prodD = dst.economyModel.productionPerDay[i];
    const double consD = dst.economyModel.consumptionPerDay[i];

    const double netSrc = prodS - consS;
    const double netNeed = consD - prodD;
    if (netSrc <= 0.0 || netNeed <= 0.0) continue;

    const double invS = srcState.inventory[i];
    const double invD = dstState.inventory[i];
    const double desiredS = std::max(1e-9, src.economyModel.desiredStock[i]);
    const double desiredD = std::max(1e-9, dst.economyModel.desiredStock[i]);

    // If source is extremely dry, don't ship it.
    if (invS < desiredS * 0.22) continue;

    const double surplus = std::max(0.0, (invS - desiredS) / desiredS);
    const double shortage = std::max(0.0, (desiredD - invD) / desiredD);

    // Score favors natural producer->consumer flows, but is modulated by current
    // surplus/shortage so inventories can influence routes.
    double score = (netSrc * netNeed);
    score *= (0.55 + 1.25 * std::min(2.0, surplus));
    score *= (0.55 + 1.25 * std::min(2.0, shortage));

    // Small randomness so it doesn't collapse to a single lane.
    score *= (0.92 + 0.16 * rng.nextDouble());

    // Insert into a tiny "top N" list.
    if (score > top[0].score) {
      top[2] = top[1];
      top[1] = top[0];
      top[0] = {cid, score};
    } else if (score > top[1].score) {
      top[2] = top[1];
      top[1] = {cid, score};
    } else if (score > top[2].score) {
      top[2] = {cid, score};
    }
  }

  // If we found nothing producer->consumer, fall back to "anything with inventory".
  if (top[0].score <= 1e-9) {
    // Try a handful of random commodities.
    for (int tries = 0; tries < 6; ++tries) {
      const econ::CommodityId cid =
          static_cast<econ::CommodityId>(rng.range<int>(0, (int)econ::kCommodityCount - 1));
      const std::size_t i = idx(cid);
      const double invS = srcState.inventory[i];
      if (invS <= 1e-6) continue;
      top[0] = {cid, 1.0};
      break;
    }
  }

  // Choose among the top 3 with weighted probability (keeps some variety).
  const double s0 = top[0].score;
  const double s1 = top[1].score;
  const double s2 = top[2].score;
  const double sum = s0 + s1 + s2;
  if (sum > 1e-9) {
    const double r = rng.nextDouble() * sum;
    if (r < s0) return top[0];
    if (r < s0 + s1) return top[1];
    return top[2];
  }

  return top[0];
}

static const Station* findStationById(const StarSystem& system, StationId id) {
  for (const auto& st : system.stations) {
    if (st.id == id) return &st;
  }
  return nullptr;
}

// Deliver all shipments whose arriveDay <= upToDay into their destination station inventories,
// removing them from the ledger afterwards.
//
// This is intentionally deterministic and safe to call frequently.
static void deliverArrivedShipments(Universe& universe,
                                    const StarSystem& system,
                                    double upToDay,
                                    TrafficLedger& ledger) {
  if (ledger.shipments.empty()) return;
  if (upToDay < 0.0) return;

  const SystemId sysId = system.stub.id;

  std::vector<TrafficShipment*> due;
  due.reserve(ledger.shipments.size());

  for (auto& sh : ledger.shipments) {
    if (sh.systemId != sysId) continue;

    // Backward compat / defensive: rebuild missing schedule fields deterministically.
    (void)hydrateShipmentScheduleFromId(sh, system, ledger.params);

    if (sh.units <= 1e-9) continue;
    if (sh.arriveDay <= upToDay) due.push_back(&sh);
  }

  if (due.empty()) return;

  std::sort(due.begin(), due.end(), [](const TrafficShipment* a, const TrafficShipment* b) {
    if (a->arriveDay != b->arriveDay) return a->arriveDay < b->arriveDay;
    return a->id < b->id;
  });

  std::unordered_set<core::u64> delivered;
  delivered.reserve(due.size() * 2 + 1);

  for (TrafficShipment* sh : due) {
    const Station* dst = findStationById(system, sh->toStation);
    if (!dst) continue;

    auto& dstState = universe.stationEconomy(*dst, sh->arriveDay);
    (void)econ::addInventory(dstState, dst->economyModel, sh->commodity, sh->units);
    delivered.insert(sh->id);
  }

  if (delivered.empty()) return;

  ledger.shipments.erase(
      std::remove_if(ledger.shipments.begin(), ledger.shipments.end(),
                     [&](const TrafficShipment& s) {
                       if (s.systemId != sysId) return false;
                       return delivered.find(s.id) != delivered.end();
                     }),
      ledger.shipments.end());
}

static void simulateNpcTradeTrafficTransitImpl(Universe& universe,
                                               const StarSystem& system,
                                               double timeDays,
                                               int& lastDay,
                                               int kMaxBackfillDays,
                                               TrafficLedger& ledger) {
  if (system.stations.size() < 2) return;
  if (timeDays < 0.0) return;

  const SystemId sysId = system.stub.id;
  const int currentDay = (int)std::floor(timeDays);

  if (currentDay <= lastDay) {
    // Still allow deliveries to resolve within the day, and keep station economies up to date.
    deliverArrivedShipments(universe, system, timeDays, ledger);
    for (const auto& st : system.stations) (void)universe.stationEconomy(st, timeDays);
    ledger.prune(timeDays);
    return;
  }

  kMaxBackfillDays = std::max(1, kMaxBackfillDays);
  int startDay = lastDay + 1;
  if (currentDay - startDay + 1 > kMaxBackfillDays) {
    startDay = currentDay - kMaxBackfillDays + 1;
  }

  // Build stable station pointers once.
  std::vector<const Station*> stations;
  stations.reserve(system.stations.size());
  for (const auto& st : system.stations) stations.push_back(&st);

  for (int day = startDay; day <= currentDay; ++day) {
    // Seed per (universe, system, day)
    core::u64 s = core::hashCombine(universe.seed(), static_cast<core::u64>(sysId));
    s = core::hashCombine(s, core::seedFromText("traffic_transit"));
    s = core::hashCombine(s, static_cast<core::u64>(day));
    core::SplitMix64 rng(s);

    // Update station economies to a stable time within the day.
    const double t = (double)day + 0.5;

    // Deliver any shipments that have arrived before this simulation slice.
    deliverArrivedShipments(universe, system, t, ledger);

    std::vector<econ::StationEconomyState*> states;
    states.reserve(stations.size());
    for (const Station* st : stations) {
      states.push_back(&universe.stationEconomy(*st, t));
    }

    const int n = (int)stations.size();
    const int baseRuns = std::clamp(1 + n, 2, 12);
    const int runs = rng.range<int>(std::max(1, baseRuns - 1), std::min(16, baseRuns + 2));

    for (int r = 0; r < runs; ++r) {
      int srcIdx = rng.range<int>(0, n - 1);
      int dstIdx = rng.range<int>(0, n - 1);
      if (n > 1) {
        int guard = 0;
        while (dstIdx == srcIdx && guard++ < 6) dstIdx = rng.range<int>(0, n - 1);
      }
      if (dstIdx == srcIdx) continue;

      const Station& src = *stations[srcIdx];
      const Station& dst = *stations[dstIdx];
      econ::StationEconomyState& srcState = *states[srcIdx];
      econ::StationEconomyState& dstState = *states[dstIdx];

      const BestCommodity pick = pickCommodity(src, srcState, dst, dstState, rng);
      const econ::CommodityId cid = pick.id;
      const std::size_t i = idx(cid);

      const double invS = srcState.inventory[i];
      const double invD = dstState.inventory[i];
      const double desiredS = std::max(1e-9, src.economyModel.desiredStock[i]);
      const double desiredD = std::max(1e-9, dst.economyModel.desiredStock[i]);
      const double capD = std::max(1e-9, dst.economyModel.capacity[i]);

      // Compute a "reasonable" shipment size.
      const double surplus = std::max(0.0, invS - desiredS * 0.25);
      const double need = std::max(0.0, (desiredD * 1.05) - invD);
      if (surplus <= 1e-6 || need <= 1e-6) continue;

      double units = std::min(surplus, need);
      units = std::min(units, desiredS * 0.45);
      units = std::min(units, capD * 0.30);
      units *= rng.range<double>(0.25, 0.85);

      if (units <= 1e-4) continue;

      const double taken = econ::takeInventory(srcState, src.economyModel, cid, units);
      if (taken <= 1e-6) continue;

      // Record in-flight shipment (do NOT credit destination inventory yet).
      ledger.record(makeNpcTradeShipment(universe.seed(), system, day, r, src, dst, cid, taken, ledger.params));
    }
  }

  // Deliver anything that has arrived by the actual current time (e.g., late-in-day arrivals).
  deliverArrivedShipments(universe, system, timeDays, ledger);

  // Advance all affected stations to the actual current time.
  for (const auto& st : system.stations) (void)universe.stationEconomy(st, timeDays);

  lastDay = currentDay;

  ledger.prune(timeDays);
}

} // namespace

void simulateNpcTradeTrafficTransit(Universe& universe,
                                    const StarSystem& system,
                                    double timeDays,
                                    std::unordered_map<SystemId, int>& lastTrafficDayBySystem,
                                    int kMaxBackfillDays,
                                    TrafficLedger* ledger) {
  if (system.stations.size() < 2) return;
  if (timeDays < 0.0) return;

  if (!ledger) {
    simulateNpcTradeTraffic(universe, system, timeDays, lastTrafficDayBySystem, kMaxBackfillDays, nullptr);
    return;
  }

  const SystemId sysId = system.stub.id;
  const int currentDay = (int)std::floor(timeDays);

  auto it = lastTrafficDayBySystem.find(sysId);
  if (it == lastTrafficDayBySystem.end()) {
    // First time we see this system: don't backfill the whole universe.
    // Match the legacy behavior by doing a tiny warm start (simulate one day).
    lastTrafficDayBySystem.emplace(sysId, std::max(-1, currentDay - 1));
    it = lastTrafficDayBySystem.find(sysId);
  }

  simulateNpcTradeTrafficTransitImpl(universe, system, timeDays, it->second, kMaxBackfillDays, *ledger);
}

void simulateNpcTradeTrafficTransit(Universe& universe,
                                    const StarSystem& system,
                                    double timeDays,
                                    std::vector<SystemTrafficStamp>& trafficStamps,
                                    int kMaxBackfillDays,
                                    TrafficLedger* ledger) {
  if (system.stations.size() < 2) return;
  if (timeDays < 0.0) return;

  if (!ledger) {
    simulateNpcTradeTraffic(universe, system, timeDays, trafficStamps, kMaxBackfillDays, nullptr);
    return;
  }

  const SystemId sysId = system.stub.id;
  const int currentDay = (int)std::floor(timeDays);

  // Find the first stamp entry for this system (keeping order stable) and also
  // compute the max day across duplicates (if any).
  SystemTrafficStamp* first = nullptr;
  int maxDay = std::numeric_limits<int>::min();
  for (auto& t : trafficStamps) {
    if (t.systemId != sysId) continue;
    if (!first) first = &t;
    maxDay = std::max(maxDay, t.dayStamp);
  }

  if (!first) {
    // Warm start (matches unordered_map overload behavior).
    trafficStamps.push_back(SystemTrafficStamp{sysId, std::max(-1, currentDay - 1)});
    first = &trafficStamps.back();
  } else {
    first->dayStamp = maxDay;
  }

  simulateNpcTradeTrafficTransitImpl(universe, system, timeDays, first->dayStamp, kMaxBackfillDays, *ledger);

  // Remove any duplicate entries for this system (keep the first occurrence).
  bool kept = false;
  for (auto it = trafficStamps.begin(); it != trafficStamps.end();) {
    if (it->systemId == sysId) {
      if (!kept) {
        kept = true;
        ++it;
      } else {
        it = trafficStamps.erase(it);
      }
    } else {
      ++it;
    }
  }
}

} // namespace stellar::sim
