#include "stellar/sim/SystemEventEconomy.h"

#include "stellar/core/Random.h"
#include "stellar/econ/Economy.h"

#include "test_harness.h"

#include <algorithm>

namespace {

constexpr std::size_t idx(stellar::econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

} // namespace

// Validates that SystemEventEconomy produces intuitively-signed directional deltas for key
// commodities, and that start-of-event inventory shocks stay within station capacity.
int test_system_event_economy() {
  int failures = 0;

  using stellar::econ::CommodityId;
  using stellar::econ::StationType;
  using stellar::sim::SystemEvent;
  using stellar::sim::SystemEventKind;

  const auto baseModel = stellar::econ::makeEconomyModel(StationType::TradeHub, 0.0);

  // --- Extra net-per-day sign checks ---
  {
    // CivilUnrest drains essentials and depresses luxury demand (inventory builds up).
    SystemEvent ev{};
    ev.active = true;
    ev.kind = SystemEventKind::CivilUnrest;
    ev.severity01 = 1.0;

    const auto extra = stellar::sim::systemEventExtraNetPerDay(baseModel, ev);

    CHECK(extra[idx(CommodityId::Food)] < 0.0);
    CHECK(extra[idx(CommodityId::Water)] < 0.0);
    CHECK(extra[idx(CommodityId::Medicine)] < 0.0);
    CHECK(extra[idx(CommodityId::Luxury)] > 0.0);
  }

  {
    // ResearchBreakthrough increases high-tech surplus and consumes raw inputs.
    SystemEvent ev{};
    ev.active = true;
    ev.kind = SystemEventKind::ResearchBreakthrough;
    ev.severity01 = 1.0;

    const auto extra = stellar::sim::systemEventExtraNetPerDay(baseModel, ev);

    CHECK(extra[idx(CommodityId::Electronics)] > 0.0);
    CHECK(extra[idx(CommodityId::Machinery)] > 0.0);
    CHECK(extra[idx(CommodityId::Ore)] < 0.0);
    CHECK(extra[idx(CommodityId::Metals)] < 0.0);
  }

  {
    // TradeBoom tends to build bulk inventory while pulling luxury.
    SystemEvent ev{};
    ev.active = true;
    ev.kind = SystemEventKind::TradeBoom;
    ev.severity01 = 1.0;

    const auto extra = stellar::sim::systemEventExtraNetPerDay(baseModel, ev);

    CHECK(extra[idx(CommodityId::Food)] > 0.0);
    CHECK(extra[idx(CommodityId::Ore)] > 0.0);
    CHECK(extra[idx(CommodityId::Luxury)] < 0.0);
  }

  // --- Start-of-event shock sanity ---
  {
    const auto model = stellar::econ::makeEconomyModel(StationType::TradeHub, 0.0);
    stellar::core::SplitMix64 rng{12345};
    auto state = stellar::econ::makeInitialState(model, rng);

    // PirateRaid steals inventory.
    {
      SystemEvent ev{};
      ev.active = true;
      ev.kind = SystemEventKind::PirateRaid;
      ev.severity01 = 1.0;

      const double before = state.inventory[idx(CommodityId::Luxury)];
      stellar::sim::applySystemEventInventoryShock(state, model, ev);
      const double after = state.inventory[idx(CommodityId::Luxury)];

      CHECK(after < before);
      CHECK(after >= 0.0);
      CHECK(after <= model.capacity[idx(CommodityId::Luxury)] + 1e-6);
    }

    // SecurityCrackdown can inject seized contraband without exceeding capacity.
    {
      SystemEvent ev{};
      ev.active = true;
      ev.kind = SystemEventKind::SecurityCrackdown;
      ev.severity01 = 1.0;

      state.inventory[idx(CommodityId::Stimulants)] = 0.0;
      const double before = state.inventory[idx(CommodityId::Stimulants)];
      stellar::sim::applySystemEventInventoryShock(state, model, ev);
      const double after = state.inventory[idx(CommodityId::Stimulants)];

      CHECK(after > before);
      CHECK(after <= model.capacity[idx(CommodityId::Stimulants)] + 1e-6);
    }
  }

  return failures;
}
