#include "stellar/sim/MarketAnalysis.h"
#include "stellar/econ/Market.h"

#include <cmath>
#include <iostream>
#include <vector>

static bool approx(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

int test_market_analysis() {
  int fails = 0;

  using stellar::econ::PricePoint;

  // --- Constant price: zero slope, zero change, zero volatility ---
  {
    std::vector<PricePoint> hist;
    for (int i = 0; i < 20; ++i) {
      hist.push_back(PricePoint{(double)i, 100.0});
    }

    const auto st = stellar::sim::analyzePriceHistory(hist, /*nowDay=*/19.0, /*windowDays=*/10.0);
    if (!st.valid) {
      std::cerr << "[test_market_analysis] expected valid stats for constant series\n";
      ++fails;
    } else {
      if (!approx(st.slopePerDay, 0.0, 1e-9)) {
        std::cerr << "[test_market_analysis] constant series slope not ~0: " << st.slopePerDay << "\n";
        ++fails;
      }
      if (!approx(st.pctChange, 0.0, 1e-9)) {
        std::cerr << "[test_market_analysis] constant series pctChange not ~0: " << st.pctChange << "\n";
        ++fails;
      }
      if (!approx(st.volatilityPerDay, 0.0, 1e-9)) {
        std::cerr << "[test_market_analysis] constant series volatility not ~0: " << st.volatilityPerDay << "\n";
        ++fails;
      }
      if (st.stdDevPrice > 1e-9) {
        std::cerr << "[test_market_analysis] constant series stdDevPrice not ~0: " << st.stdDevPrice << "\n";
        ++fails;
      }
    }
  }

  // --- Linear trend: slope ~= 10, r2 ~= 1, pctChange ~= 100% ---
  {
    std::vector<PricePoint> hist;
    for (int i = 0; i <= 10; ++i) {
      hist.push_back(PricePoint{(double)i, 100.0 + 10.0 * (double)i});
    }

    const auto st = stellar::sim::analyzePriceHistory(hist, /*nowDay=*/10.0, /*windowDays=*/0.0);
    if (!st.valid) {
      std::cerr << "[test_market_analysis] expected valid stats for linear series\n";
      ++fails;
    } else {
      if (!approx(st.slopePerDay, 10.0, 1e-9)) {
        std::cerr << "[test_market_analysis] linear series slope mismatch: " << st.slopePerDay << "\n";
        ++fails;
      }
      if (st.r2 < 0.999999) {
        std::cerr << "[test_market_analysis] linear series r2 too low: " << st.r2 << "\n";
        ++fails;
      }
      if (!approx(st.pctChange, 100.0, 1e-9)) {
        std::cerr << "[test_market_analysis] linear series pctChange mismatch: " << st.pctChange << "\n";
        ++fails;
      }
    }
  }

  // --- Forecast: positive net production should reduce price ---
  {
    using namespace stellar::econ;
    StationEconomyModel m{};
    StationEconomyState s{};

    for (std::size_t i = 0; i < kCommodityCount; ++i) {
      m.capacity[i] = 200.0;
      m.desiredStock[i] = 100.0;
      m.productionPerDay[i] = 0.0;
      m.consumptionPerDay[i] = 0.0;
      s.inventory[i] = 0.0;
    }

    const auto cid = CommodityId::Food;
    const std::size_t iFood = static_cast<std::size_t>(cid);

    m.priceVolatility = 1.0;
    m.productionPerDay[iFood] = 50.0;
    m.consumptionPerDay[iFood] = 0.0;

    // Start below desired -> high price.
    s.inventory[iFood] = 50.0;

    const double midNow = quote(s, m, cid, 0.10).mid;
    const double midForecast = stellar::sim::forecastMidPrice(s, m, cid, 1.0);

    if (!(midForecast + 1e-9 < midNow)) {
      std::cerr << "[test_market_analysis] expected forecast mid < current mid when net production positive (" << midForecast << " vs " << midNow << ")\n";
      ++fails;
    }
  }

  // --- Forecast: net consumption should increase price ---
  {
    using namespace stellar::econ;
    StationEconomyModel m{};
    StationEconomyState s{};

    for (std::size_t i = 0; i < kCommodityCount; ++i) {
      m.capacity[i] = 200.0;
      m.desiredStock[i] = 100.0;
      m.productionPerDay[i] = 0.0;
      m.consumptionPerDay[i] = 0.0;
      s.inventory[i] = 0.0;
    }

    const auto cid = CommodityId::Food;
    const std::size_t iFood = static_cast<std::size_t>(cid);

    m.priceVolatility = 1.0;
    m.productionPerDay[iFood] = 0.0;
    m.consumptionPerDay[iFood] = 40.0;

    // Start at desired -> baseline price.
    s.inventory[iFood] = 100.0;

    const double midNow = quote(s, m, cid, 0.10).mid;
    const double midForecast = stellar::sim::forecastMidPrice(s, m, cid, 1.0);

    if (!(midForecast > midNow + 1e-9)) {
      std::cerr << "[test_market_analysis] expected forecast mid > current mid when net consumption positive (" << midForecast << " vs " << midNow << ")\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_market_analysis] PASS\n";
  }

  return fails;
}
