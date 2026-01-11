#include "stellar/sim/MarketAnalysis.h"

#include "stellar/math/Math.h"
#include "stellar/econ/Commodity.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace stellar::sim {

static bool finite(double x) {
  return std::isfinite(x);
}

PriceTrendStats analyzePriceHistory(const std::vector<stellar::econ::PricePoint>& hist,
                                   double nowDay,
                                   double windowDays) {
  PriceTrendStats out{};

  if (hist.empty()) return out;
  if (!finite(nowDay)) nowDay = hist.back().day;
  if (!finite(windowDays)) windowDays = 0.0;

  const double minDay = (windowDays > 0.0) ? (nowDay - windowDays) : -std::numeric_limits<double>::infinity();

  // Locate the first sample in-range.
  std::size_t start = 0;
  while (start < hist.size() && hist[start].day < minDay) ++start;
  if (start >= hist.size()) return out;

  // Locate the last sample at/before nowDay.
  std::size_t end = hist.size();
  while (end > start && hist[end - 1].day > nowDay) --end;
  if (end <= start) return out;

  const std::size_t n = end - start;
  out.samples = n;

  // Basic stats.
  double minP = std::numeric_limits<double>::infinity();
  double maxP = -std::numeric_limits<double>::infinity();
  double sum = 0.0;
  std::size_t validCount = 0;

  for (std::size_t i = start; i < end; ++i) {
    const double p = hist[i].price;
    if (!finite(p)) continue;
    minP = std::min(minP, p);
    maxP = std::max(maxP, p);
    sum += p;
    ++validCount;
  }

  if (validCount < 2) {
    // Still return "valid" if we have at least one finite sample.
    if (validCount == 1) {
      out.valid = true;
      out.firstDay = hist[start].day;
      out.lastDay = hist[end - 1].day;
      out.firstPrice = hist[start].price;
      out.lastPrice = hist[end - 1].price;
      out.minPrice = hist[start].price;
      out.maxPrice = hist[start].price;
      out.meanPrice = hist[start].price;
      out.stdDevPrice = 0.0;
      out.slopePerDay = 0.0;
      out.r2 = 0.0;
      out.pctChange = 0.0;
      out.volatilityPerDay = 0.0;
    }
    return out;
  }

  const double mean = sum / (double)validCount;

  // Variance of price (for z-score / window summary).
  double var = 0.0;
  for (std::size_t i = start; i < end; ++i) {
    const double p = hist[i].price;
    if (!finite(p)) continue;
    const double d = p - mean;
    var += d * d;
  }
  var /= std::max(1.0, (double)validCount - 1.0);

  // Linear regression on (day, price).
  double meanX = 0.0;
  double meanY = 0.0;
  {
    double sx = 0.0;
    double sy = 0.0;
    std::size_t c = 0;
    for (std::size_t i = start; i < end; ++i) {
      const double x = hist[i].day;
      const double y = hist[i].price;
      if (!finite(x) || !finite(y)) continue;
      sx += x;
      sy += y;
      ++c;
    }
    if (c >= 2) {
      meanX = sx / (double)c;
      meanY = sy / (double)c;
    }
  }

  double sxx = 0.0;
  double sxy = 0.0;
  double syy = 0.0;
  std::size_t regCount = 0;
  for (std::size_t i = start; i < end; ++i) {
    const double x = hist[i].day;
    const double y = hist[i].price;
    if (!finite(x) || !finite(y)) continue;

    const double dx = x - meanX;
    const double dy = y - meanY;
    sxx += dx * dx;
    sxy += dx * dy;
    syy += dy * dy;
    ++regCount;
  }

  double slope = 0.0;
  double r2 = 0.0;
  if (regCount >= 2 && sxx > 1e-12) {
    slope = sxy / sxx;

    // r^2 = (cov^2) / (varX*varY)
    const double denom = sxx * syy;
    if (denom > 1e-18) {
      const double r = sxy / std::sqrt(denom);
      r2 = r * r;
    }
  }

  // Volatility estimate from log returns per day.
  double vol = 0.0;
  {
    std::vector<double> rs;
    rs.reserve(n);

    for (std::size_t i = start + 1; i < end; ++i) {
      const double p0 = hist[i - 1].price;
      const double p1 = hist[i].price;
      const double t0 = hist[i - 1].day;
      const double t1 = hist[i].day;

      const double dt = t1 - t0;
      if (!finite(p0) || !finite(p1) || !finite(dt)) continue;
      if (dt <= 1e-9) continue;
      if (p0 <= 1e-12 || p1 <= 1e-12) continue;

      const double r = std::log(p1 / p0) / dt;
      if (!finite(r)) continue;
      rs.push_back(r);
    }

    if (rs.size() >= 2) {
      double m = 0.0;
      for (double r : rs) m += r;
      m /= (double)rs.size();

      double vv = 0.0;
      for (double r : rs) {
        const double d = r - m;
        vv += d * d;
      }
      vv /= std::max(1.0, (double)rs.size() - 1.0);
      vol = std::sqrt(std::max(0.0, vv));
    } else {
      vol = 0.0;
    }
  }

  // Assemble.
  out.valid = true;
  out.firstDay = hist[start].day;
  out.lastDay = hist[end - 1].day;
  out.firstPrice = hist[start].price;
  out.lastPrice = hist[end - 1].price;

  out.minPrice = minP;
  out.maxPrice = maxP;
  out.meanPrice = mean;
  out.stdDevPrice = std::sqrt(std::max(0.0, var));

  out.slopePerDay = slope;
  out.r2 = r2;

  if (finite(out.firstPrice) && std::abs(out.firstPrice) > 1e-12) {
    out.pctChange = (out.lastPrice - out.firstPrice) / out.firstPrice * 100.0;
  } else {
    out.pctChange = 0.0;
  }

  out.volatilityPerDay = vol;

  return out;
}

static constexpr std::size_t cidx(stellar::econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

double forecastMidPrice(const stellar::econ::StationEconomyState& state,
                        const stellar::econ::StationEconomyModel& model,
                        stellar::econ::CommodityId commodity,
                        double horizonDays) {
  if (!finite(horizonDays)) horizonDays = 0.0;
  horizonDays = std::max(0.0, horizonDays);

  const std::size_t i = cidx(commodity);
  const double cap = std::max(0.0, model.capacity[i]);

  double inv = state.inventory[i];
  if (!finite(inv)) inv = 0.0;
  inv = std::clamp(inv, 0.0, cap);

  const double netPerDay = model.productionPerDay[i] - model.consumptionPerDay[i];
  double inv2 = inv + netPerDay * horizonDays;
  inv2 = std::clamp(inv2, 0.0, cap);

  // Mirror econ::midPrice() logic for a single commodity.
  const auto& def = stellar::econ::commodityDef(commodity);
  const double desired = model.desiredStock[i];
  if (desired <= 1e-9) return def.basePrice;

  const double ratio = inv2 / desired;
  double factor = 1.0 + model.priceVolatility * (1.0 - ratio);

  if (inv2 < desired * 0.05) factor *= 1.4;

  factor = stellar::math::clamp(factor, 0.2, 5.0);
  return def.basePrice * factor;
}

} // namespace stellar::sim
