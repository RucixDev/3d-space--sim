#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/Celestial.h"

#include <array>

namespace stellar::proc {

// A lightweight, deterministic “macro‑economy fingerprint” for a star system.
//
// Design goals:
//  - Galaxy‑coherent: nearby systems tend to have similar profiles.
//  - Seed‑stable: changes only when the universe seed changes.
//  - Cheap: suitable for bulk queries (UI tables, route analysis).
//
// The profile is intentionally generic. Future rounds can use it to drive:
//  - station archetype selection (mining/agri/industrial hubs)
//  - per‑commodity supply/demand baselines
//  - piracy/law enforcement ambient behavior
struct TradeProfile {
  // Coarse macro factors, each normalized to [0,1].
  double resources{0.0};    // ore/metals/fuel potential
  double agriculture{0.0};  // food/water surplus potential
  double industry{0.0};     // manufacturing depth
  double technology{0.0};   // research/advanced goods bias
  double hub{0.0};          // connectivity / “trade gravity”
  double lawlessness{0.0};  // piracy / volatility
  double population{0.0};   // market depth
  double wealth{0.0};       // luxury/tech purchasing power

  // Commodity‑level export/import desirability in [0,1].
  std::array<double, econ::kCommodityCount> exportScore{};
  std::array<double, econ::kCommodityCount> importScore{};

  // A deterministic signature (useful for caching/UI diffing).
  core::u64 signature{0};
};

// Generate a galaxy‑coherent trade profile for a given system stub.
//
// Note: This does NOT mutate the universe/economy; it is purely procedural.
TradeProfile generateTradeProfile(core::u64 universeSeed, const sim::SystemStub& stub);

} // namespace stellar::proc
