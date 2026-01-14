#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Economy.h"
#include "stellar/proc/TradeProfile.h"

namespace stellar::proc {

// Apply a system TradeProfile to a station's economy model.
//
// This is deliberately gameplay-oriented (not a realistic macro-economy):
//  - high hub/population => deeper inventories + lower volatility
//  - high exportScore => stronger production for already-produced goods
//  - high importScore => stronger consumption for already-consumed goods
//  - high lawlessness => higher shock volatility
//
// The tuning is deterministic per station via `stationSeed`.
econ::StationEconomyModel tuneEconomyModelForTradeProfile(core::u64 stationSeed,
                                                         const econ::StationEconomyModel& base,
                                                         const TradeProfile& profile);

// Suggest a local station fee rate (tax / tariff) from a base faction fee.
//
// Deterministic per station via `stationSeed`. Used in proc-gen so trade hubs can
// feel more competitive and lawless systems can feel "riskier".
double tuneStationFeeRateForTradeProfile(core::u64 stationSeed,
                                        double baseFeeRate,
                                        econ::StationType stationType,
                                        const TradeProfile& profile);

} // namespace stellar::proc
