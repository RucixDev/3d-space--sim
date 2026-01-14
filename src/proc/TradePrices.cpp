#include "stellar/proc/TradePrices.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {

static inline double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

static inline double clampSigned1(double v) {
  if (v < -1.0) return -1.0;
  if (v > 1.0) return 1.0;
  return v;
}

std::array<double, econ::kCommodityCount>
estimateTradePriceMultipliers(core::u64 systemSeed,
                              const TradeProfile& profile,
                              const TradePriceModelParams& params) {
  std::array<double, econ::kCommodityCount> out{};

  // Defensive clamps.
  TradePriceModelParams p = params;
  if (!std::isfinite(p.logScarcity)) p.logScarcity = 0.55;
  if (!std::isfinite(p.ampBase)) p.ampBase = 0.95;
  if (!std::isfinite(p.ampHub)) p.ampHub = 0.55;
  if (!std::isfinite(p.ampPop)) p.ampPop = 0.35;
  if (!std::isfinite(p.ampMin)) p.ampMin = 0.20;
  if (!std::isfinite(p.ampMax)) p.ampMax = 0.95;
  if (!std::isfinite(p.jitter)) p.jitter = 0.0;
  if (!std::isfinite(p.minMul)) p.minMul = 0.55;
  if (!std::isfinite(p.maxMul)) p.maxMul = 1.80;

  p.jitter = std::clamp(p.jitter, 0.0, 0.25);
  if (p.minMul > p.maxMul) std::swap(p.minMul, p.maxMul);
  p.ampMin = std::max(0.0, p.ampMin);
  p.ampMax = std::max(p.ampMin, p.ampMax);

  const double hub = clamp01(profile.hub);
  const double pop = clamp01(profile.population);

  const double amp = std::clamp(p.ampBase - p.ampHub * hub - p.ampPop * pop, p.ampMin, p.ampMax);

  // Stable micro-variation base seed.
  const core::u64 baseSeed = core::hashCombine(systemSeed, core::fnv1a64("trade_prices_v1"));

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double expS = clamp01(profile.exportScore[i]);
    const double impS = clamp01(profile.importScore[i]);
    const double scarcity = clampSigned1(impS - expS);

    // Raw multiplicative factor.
    const double raw = std::exp(p.logScarcity * scarcity);

    // Compress extremes for large hub/pop systems.
    double mul = 1.0 + (raw - 1.0) * amp;

    // Per-commodity micro jitter.
    if (p.jitter > 0.0) {
      core::SplitMix64 rng(core::hashCombine(baseSeed, static_cast<core::u64>(i) * 0x9E3779B97F4A7C15ull));
      const double j = rng.nextDouble() * 2.0 - 1.0; // [-1, 1]
      mul *= (1.0 + p.jitter * j);
    }

    mul = stellar::math::clamp(mul, p.minMul, p.maxMul);
    out[i] = mul;
  }

  return out;
}

double estimateTradePriceMultiplier(core::u64 systemSeed,
                                    const TradeProfile& profile,
                                    econ::CommodityId commodity,
                                    const TradePriceModelParams& params) {
  const auto m = estimateTradePriceMultipliers(systemSeed, profile, params);
  return m[static_cast<std::size_t>(commodity)];
}

TradePriceQuote estimateTradePriceQuote(core::u64 systemSeed,
                                        const TradeProfile& profile,
                                        econ::CommodityId commodity,
                                        double bidAskSpread,
                                        const TradePriceModelParams& params) {
  TradePriceQuote q{};
  q.multiplier = estimateTradePriceMultiplier(systemSeed, profile, commodity, params);

  const auto& def = econ::commodityDef(commodity);
  q.mid = def.basePrice * q.multiplier;

  if (!std::isfinite(bidAskSpread)) bidAskSpread = 0.0;
  bidAskSpread = std::clamp(bidAskSpread, 0.0, 1.0);
  const double half = bidAskSpread * 0.5;
  q.ask = q.mid * (1.0 + half);
  q.bid = q.mid * (1.0 - half);
  return q;
}

} // namespace stellar::proc
