#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Universe.h"

#include <unordered_map>
#include <vector>

namespace stellar::sim {

struct TrafficLedger; // fwd (required for transit-mode)

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

} // namespace stellar::sim
