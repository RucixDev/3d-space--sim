#include "stellar/sim/SystemEventEconomy.h"

#include "stellar/core/Random.h"
#include "stellar/econ/Economy.h"

#include <catch2/catch_test_macros.hpp>

namespace {

constexpr std::size_t idx(stellar::econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

} // namespace

TEST_CASE("SystemEventEconomy: extra net per day matches expected sign for key commodities") {
  using stellar::econ::CommodityId;
  using stellar::econ::StationType;
  using stellar::sim::SystemEvent;
  using stellar::sim::SystemEventKind;

  const auto baseModel = stellar::econ::makeEconomyModel(StationType::TradeHub, 0.0);

  SECTION("CivilUnrest drains essentials and depresses luxury") {
    SystemEvent ev{};
    ev.active = true;
    ev.kind = SystemEventKind::CivilUnrest;
    ev.severity01 = 1.0;

    const auto extra = stellar::sim::systemEventExtraNetPerDay(baseModel, ev);

    REQUIRE(extra[idx(CommodityId::Food)] < 0.0);
    REQUIRE(extra[idx(CommodityId::Water)] < 0.0);
    REQUIRE(extra[idx(CommodityId::Medicine)] < 0.0);
    // Unrest: luxury tends to accumulate (demand drops) => positive net.
    REQUIRE(extra[idx(CommodityId::Luxury)] > 0.0);
  }

  SECTION("ResearchBreakthrough increases high-tech surplus and consumes raw inputs") {
    SystemEvent ev{};
    ev.active = true;
    ev.kind = SystemEventKind::ResearchBreakthrough;
    ev.severity01 = 1.0;

    const auto extra = stellar::sim::systemEventExtraNetPerDay(baseModel, ev);

    REQUIRE(extra[idx(CommodityId::Electronics)] > 0.0);
    REQUIRE(extra[idx(CommodityId::Machinery)] > 0.0);
    REQUIRE(extra[idx(CommodityId::Ore)] < 0.0);
    REQUIRE(extra[idx(CommodityId::Metals)] < 0.0);
  }

  SECTION("TradeBoom tends to build bulk inventory while pulling luxury") {
    SystemEvent ev{};
    ev.active = true;
    ev.kind = SystemEventKind::TradeBoom;
    ev.severity01 = 1.0;

    const auto extra = stellar::sim::systemEventExtraNetPerDay(baseModel, ev);

    REQUIRE(extra[idx(CommodityId::Food)] > 0.0);
    REQUIRE(extra[idx(CommodityId::Ore)] > 0.0);
    REQUIRE(extra[idx(CommodityId::Luxury)] < 0.0);
  }
}

TEST_CASE("SystemEventEconomy: start-of-event shocks stay within capacity") {
  using stellar::econ::CommodityId;
  using stellar::econ::StationType;
  using stellar::sim::SystemEvent;
  using stellar::sim::SystemEventKind;

  const auto model = stellar::econ::makeEconomyModel(StationType::TradeHub, 0.0);
  stellar::core::SplitMix64 rng{12345};
  auto state = stellar::econ::makeInitialState(model, rng);

  SECTION("PirateRaid steals inventory") {
    SystemEvent ev{};
    ev.active = true;
    ev.kind = SystemEventKind::PirateRaid;
    ev.severity01 = 1.0;

    const double before = state.inventory[idx(CommodityId::Luxury)];
    stellar::sim::applySystemEventInventoryShock(state, model, ev);
    const double after = state.inventory[idx(CommodityId::Luxury)];

    REQUIRE(after < before);
    REQUIRE(after >= 0.0);
    REQUIRE(after <= model.capacity[idx(CommodityId::Luxury)] + 1e-6);
  }

  SECTION("SecurityCrackdown can inject seized contraband without exceeding capacity") {
    SystemEvent ev{};
    ev.active = true;
    ev.kind = SystemEventKind::SecurityCrackdown;
    ev.severity01 = 1.0;

    state.inventory[idx(CommodityId::Stimulants)] = 0.0;
    const double before = state.inventory[idx(CommodityId::Stimulants)];
    stellar::sim::applySystemEventInventoryShock(state, model, ev);
    const double after = state.inventory[idx(CommodityId::Stimulants)];

    REQUIRE(after > before);
    REQUIRE(after <= model.capacity[idx(CommodityId::Stimulants)] + 1e-6);
  }
}

int test_system_event_economy() {
  return 0;
}
