#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/proc/TradeProfile.h"

#include <array>

namespace stellar::proc {

// -----------------------------------------------------------------------------
// Procedural macro price model (TradeProfile -> price multipliers)
// -----------------------------------------------------------------------------
//
// TradeProfile captures a galaxy-coherent supply/demand fingerprint for each
// commodity via exportScore/importScore.
//
// This module turns that signal into a *deterministic* per-system price
// multiplier field that can be used for:
//  - cheap "macro" trade planning (without simulating station inventories)
//  - UI overlays / intel probes
//  - later: seeding NPC logistics volumes / missions
//
// It is intentionally stylized and gameplay-oriented, not a macro-econ sim.

struct TradePriceModelParams {
  // Scarcity is measured as (importScore - exportScore) in [-1, +1].
  // We map this to a multiplicative factor using exp(logScarcity * scarcity).
  //
  // logScarcity=0.55 yields a raw range of ~[0.58, 1.73].
  double logScarcity{0.55};

  // Damp price extremes for high-pop / high-hub systems.
  // Effective amplitude: amp = clamp(ampBase - ampHub*hub - ampPop*population).
  double ampBase{0.95};
  double ampHub{0.55};
  double ampPop{0.35};
  double ampMin{0.20};
  double ampMax{0.95};

  // Small per-commodity micro-variation (multiplicative), deterministic from
  // systemSeed + commodity index. Value is the maximum absolute deviation.
  // Example: jitter=0.04 -> ±4%.
  double jitter{0.04};

  // Clamp multipliers to keep the macro model sane.
  double minMul{0.55};
  double maxMul{1.80};
};

// Deterministic price multipliers for all commodities in a system.
//
// `systemSeed` should be stable for the system (SystemStub::seed is ideal).
std::array<double, econ::kCommodityCount>
estimateTradePriceMultipliers(core::u64 systemSeed,
                              const TradeProfile& profile,
                              const TradePriceModelParams& params = {});

// Convenience: compute a single commodity multiplier.
double estimateTradePriceMultiplier(core::u64 systemSeed,
                                    const TradeProfile& profile,
                                    econ::CommodityId commodity,
                                    const TradePriceModelParams& params = {});

// A macro "market quote" derived from the TradeProfile field.
//
// Prices are in credits per unit.
struct TradePriceQuote {
  double mid{0.0};
  double ask{0.0};
  double bid{0.0};
  double multiplier{1.0};
};

TradePriceQuote estimateTradePriceQuote(core::u64 systemSeed,
                                        const TradeProfile& profile,
                                        econ::CommodityId commodity,
                                        double bidAskSpread = 0.10,
                                        const TradePriceModelParams& params = {});

} // namespace stellar::proc
