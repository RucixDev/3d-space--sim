#include "stellar/proc/TradeEconomy.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/econ/Commodity.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {

static constexpr std::size_t cidx(econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

static double clamp01(double x) {
  return std::clamp(x, 0.0, 1.0);
}

static double stableJitter01(core::u64 seed) {
  core::SplitMix64 rng(seed);
  return rng.nextDouble();
}

static double stationTypeMarketMul(econ::StationType t, const TradeProfile& p) {
  using econ::StationType;
  // Base market depth is driven by population + hub-ness.
  double m = 0.55 + 1.40 * clamp01(p.population) + 0.80 * clamp01(p.hub);

  // Station-type specific tuning. Keep within sane bounds.
  switch (t) {
    case StationType::TradeHub:
      m *= 1.10 + 0.70 * clamp01(p.hub);
      break;
    case StationType::Shipyard:
      m *= 0.95 + 0.35 * clamp01(p.industry);
      break;
    case StationType::Industrial:
      m *= 0.95 + 0.30 * clamp01(p.industry);
      break;
    case StationType::Refinery:
      m *= 0.92 + 0.35 * clamp01(p.industry);
      break;
    case StationType::Research:
      m *= 0.88 + 0.40 * clamp01(p.technology);
      break;
    case StationType::Agricultural:
      m *= 0.90 + 0.20 * clamp01(p.population);
      break;
    case StationType::Mining:
      m *= 0.88 + 0.20 * clamp01(p.population);
      break;
    case StationType::Outpost:
    default:
      m *= 0.72 + 0.28 * clamp01(p.population);
      break;
  }

  return std::clamp(m, 0.35, 3.25);
}

econ::StationEconomyModel tuneEconomyModelForTradeProfile(core::u64 stationSeed,
                                                         const econ::StationEconomyModel& base,
                                                         const TradeProfile& profile) {
  econ::StationEconomyModel m = base;

  const core::u64 seed0 = core::hashCombine(stationSeed, core::seedFromText("trade_econ"));
  core::SplitMix64 rng(seed0);

  const double marketMul = stationTypeMarketMul(m.type, profile);

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double ex = clamp01(profile.exportScore[i]);
    const double im = clamp01(profile.importScore[i]);

    // Mild commodity-specific stock preference, then a small station micro-variance.
    const double desireBias = 0.70 + 0.45 * std::max(ex, im); // [0.70, 1.15]
    const double microStock = 0.92 + 0.16 * rng.nextDouble(); // [0.92, 1.08]
    const double stockMul = marketMul * desireBias * microStock;

    m.desiredStock[i] *= stockMul;
    m.capacity[i] *= stockMul;

    // Preserve a reasonable cap ratio even if upstream models change.
    m.capacity[i] = std::max(m.capacity[i], m.desiredStock[i] * 1.05);

    // Production/consumption modulation based on export/import desirability.
    const double microFlow = 0.94 + 0.12 * rng.nextDouble(); // [0.94, 1.06]

    if (m.productionPerDay[i] > 0.0) {
      const double prodMul = (0.75 + 0.90 * ex) * microFlow;
      m.productionPerDay[i] *= prodMul;
    }
    if (m.consumptionPerDay[i] > 0.0) {
      const double consMul = (0.75 + 0.90 * im) * microFlow;
      m.consumptionPerDay[i] *= consMul;
    }
  }

  // Volatility: bigger, more connected markets tend to be more stable.
  // Also allow lawlessness to increase drift (shocks).
  {
    const double hub = clamp01(profile.hub);
    const double pop = clamp01(profile.population);
    const double law = clamp01(profile.lawlessness);

    double volMul = 1.20 - 0.55 * hub - 0.45 * pop;
    if (m.type == econ::StationType::TradeHub) volMul *= 0.85;
    if (m.type == econ::StationType::Outpost) volMul *= 1.10;
    volMul *= (0.96 + 0.08 * stableJitter01(core::hashCombine(seed0, 0xA11CEull)));

    m.priceVolatility = std::clamp(m.priceVolatility * volMul, 0.30, 1.55);

    double shockMul = 0.65 + 1.70 * law;
    shockMul *= (1.05 - 0.40 * hub - 0.25 * pop);
    shockMul = std::clamp(shockMul, 0.25, 2.75);
    m.shockVolatility = std::max(0.0, m.shockVolatility * shockMul);
  }

  return m;
}

double tuneStationFeeRateForTradeProfile(core::u64 stationSeed,
                                        double baseFeeRate,
                                        econ::StationType stationType,
                                        const TradeProfile& profile) {
  double fee = std::max(0.0, baseFeeRate);
  const double hub = clamp01(profile.hub);
  const double pop = clamp01(profile.population);
  const double law = clamp01(profile.lawlessness);

  // Competitive hubs tend to compress fees; lawless systems inflate them.
  double mul = 1.0;
  mul *= (1.05 - 0.35 * hub - 0.20 * pop);
  mul *= (1.00 + 0.35 * law);

  // Station type adjustments.
  if (stationType == econ::StationType::TradeHub) mul *= 0.85;
  if (stationType == econ::StationType::Outpost) mul *= 1.10;

  // Micro variation to avoid perfectly uniform taxes.
  const double jitter = stableJitter01(core::hashCombine(stationSeed, core::seedFromText("fee")));
  mul *= (0.95 + 0.10 * jitter);

  fee *= mul;
  return std::clamp(fee, 0.0, 0.25);
}

} // namespace stellar::proc
