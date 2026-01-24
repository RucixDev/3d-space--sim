#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Universe.h"

#include <unordered_map>
#include <vector>

namespace stellar::sim {

struct TrafficLedger; // fwd (required for transit-mode)

// ----------------------------------------------------------------------------
// Traffic model selection
// ----------------------------------------------------------------------------
//
// The original NPC traffic model (SupplyDemand) nudges markets by shipping
// goods from deterministic producers toward deterministic consumers.
//
// For more "space-trader"-like behavior, MarketArbitrage routes cargo using
// current market quotes (bid/ask + station fees), so traffic tends to follow
// profit gradients (surplus -> shortage) and reacts quickly to player-driven
// supply shocks.
//
// Hybrid blends both: it follows profit when present, but still maintains
// producer->consumer baseline flows so low-liquidity commodities continue to
// move.
enum class TrafficTradeModel : core::u8 {
  SupplyDemand = 0,
  MarketArbitrage = 1,
  Hybrid = 2,
};

// Human-friendly label (stable string literal).
const char* trafficTradeModelName(TrafficTradeModel m);

// Clamp/convert from an int (e.g., from config / UI).
TrafficTradeModel trafficTradeModelFromInt(core::i64 v);

// Tuning parameters for the transit-mode traffic sim.
//
// Notes:
//  - Parameters are intentionally "soft" (weights / clamps) so saves remain
//    compatible and so the model can be tuned live via CVars.
//  - The defaults are conservative and aim to keep the economy stable.
struct TrafficTransitTradeParams {
  TrafficTradeModel model{TrafficTradeModel::SupplyDemand};

  // Used by MarketArbitrage/Hybrid when converting mid prices to bid/ask quotes.
  // Range: [0, 1].
  double bidAskSpread{0.12};

  // Minimum per-unit profit required for a cargo to be considered in
  // MarketArbitrage/Hybrid. Set <= 0 to disable the threshold.
  double minProfitPerUnitCr{0.0};

  // Distance penalty (km^-1). Arbitrage scores are divided by:
  //   (1 + distancePenaltyPerKm * distKm)
  // This prevents far-apart stations from dominating purely due to price
  // gradients when the system has many stations spread across large orbits.
  double distancePenaltyPerKm{1e-8};

  // Hybrid weights: score = profitWeight * profitScore + flowWeight * flowScore
  double profitWeight{1.0};
  double flowWeight{0.35};

  // Shipment sizing clamps (fractions).
  double maxUnitsFracOfSrcDesired{0.45};
  double maxUnitsFracOfDstCapacity{0.30};

  // Random shipment size scale (uniform range). Must be within [0, 1].
  double randomUnitsMinFrac{0.25};
  double randomUnitsMaxFrac{0.85};
};

// Transit-mode variant of NPC trade traffic.
//
// Unlike simulateNpcTradeTraffic(), which applies inventory deltas instantly,
// this function records shipments into a TrafficLedger and only credits cargo
// to the destination once the shipment's arrival time has been reached.
//
// The ledger therefore represents "in-flight" cargo that can be visualized as
// convoys and interdicted by gameplay systems.
void simulateNpcTradeTrafficTransit(Universe& universe,
                                    const StarSystem& system,
                                    double timeDays,
                                    std::unordered_map<SystemId, int>& lastTrafficDayBySystem,
                                    int kMaxBackfillDays = 14,
                                    TrafficLedger* ledger = nullptr);

// SaveGame-friendly overload.
//
// Convenience for callers that already store traffic stamps in the SaveGame format.
// The vector will be updated in-place (the matching system entry is added/updated).
void simulateNpcTradeTrafficTransit(Universe& universe,
                                    const StarSystem& system,
                                    double timeDays,
                                    std::vector<SystemTrafficStamp>& trafficStamps,
                                    int kMaxBackfillDays = 14,
                                    TrafficLedger* ledger = nullptr);

// Extended API that allows selecting the traffic model and tuning weights.
//
// Existing callers can keep using simulateNpcTradeTrafficTransit(); this is
// an opt-in path for gameplay/tools.
void simulateNpcTradeTrafficTransitEx(Universe& universe,
                                      const StarSystem& system,
                                      double timeDays,
                                      std::unordered_map<SystemId, int>& lastTrafficDayBySystem,
                                      const TrafficTransitTradeParams& params,
                                      int kMaxBackfillDays = 14,
                                      TrafficLedger* ledger = nullptr);

void simulateNpcTradeTrafficTransitEx(Universe& universe,
                                      const StarSystem& system,
                                      double timeDays,
                                      std::vector<SystemTrafficStamp>& trafficStamps,
                                      const TrafficTransitTradeParams& params,
                                      int kMaxBackfillDays = 14,
                                      TrafficLedger* ledger = nullptr);

} // namespace stellar::sim
