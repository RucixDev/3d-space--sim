#pragma once

#include "stellar/econ/Economy.h"

#include <vector>

namespace stellar::sim {

// Summary statistics for a slice of a commodity price history.
//
// This is intended for UI (dashboards, tooltips) and lightweight gameplay logic.
struct PriceTrendStats {
  bool valid{false};

  // Window bounds actually used (clamped to available history).
  double firstDay{0.0};
  double lastDay{0.0};

  double firstPrice{0.0};
  double lastPrice{0.0};

  double minPrice{0.0};
  double maxPrice{0.0};
  double meanPrice{0.0};
  double stdDevPrice{0.0};

  // Linear regression slope (credits/unit/day) and goodness-of-fit.
  double slopePerDay{0.0};
  double r2{0.0};

  // Percent change over the window.
  double pctChange{0.0};

  // Volatility estimate: standard deviation of log-returns per day.
  // Example: 0.02 ~= 2%/day.
  double volatilityPerDay{0.0};

  std::size_t samples{0};
};

// Analyze a price history series, considering only points in the last `windowDays`.
//
// - `hist` is typically StationEconomyState::history[commodity]
// - `nowDay` is the current simulation day (timeDays)
// - If `windowDays <= 0`, the full history is used.
PriceTrendStats analyzePriceHistory(const std::vector<stellar::econ::PricePoint>& hist,
                                   double nowDay,
                                   double windowDays);

// Deterministic "expected" mid-price forecast at horizonDays.
//
// Uses only the station model's net production/consumption (mean shock = 0) to
// advance inventory, then evaluates the mid price curve. This is not intended
// to be perfectly accurate, but gives the player a meaningful directional hint.
//
// Returns the forecast mid price in credits/unit.
double forecastMidPrice(const stellar::econ::StationEconomyState& state,
                        const stellar::econ::StationEconomyModel& model,
                        stellar::econ::CommodityId commodity,
                        double horizonDays);

} // namespace stellar::sim
