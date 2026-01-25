#include "stellar/sim/SystemEventEconomy.h"

#include "stellar/econ/Commodity.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

namespace stellar::sim {

static constexpr std::size_t idx(stellar::econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

// Weights encode *net inventory pressure* for each commodity:
//   negative => inventory tends to drop => prices tend to rise
//   positive => inventory tends to rise => prices tend to fall
static std::array<double, econ::kCommodityCount> weightsFor(SystemEventKind k) {
  std::array<double, econ::kCommodityCount> w{};
  w.fill(0.0);

  auto W = [&](econ::CommodityId id, double v) { w[idx(id)] = v; };

  switch (k) {
    case SystemEventKind::None:
      break;

    case SystemEventKind::TradeBoom:
      // Healthy logistics: essentials become plentiful, luxury/high-tech get snapped up.
      W(econ::CommodityId::Food,       +0.12);
      W(econ::CommodityId::Water,      +0.12);
      W(econ::CommodityId::Fuel,       +0.08);
      W(econ::CommodityId::Medicine,   +0.06);
      W(econ::CommodityId::Metals,     +0.10);
      W(econ::CommodityId::Ore,        +0.10);
      W(econ::CommodityId::Electronics,-0.18);
      W(econ::CommodityId::Machinery,  -0.18);
      W(econ::CommodityId::Luxury,     -0.35);
      W(econ::CommodityId::Stimulants, -0.10);
      W(econ::CommodityId::Weapons,    -0.05);
      break;

    case SystemEventKind::TradeBust:
      // Import slump: essentials get tighter, luxury stacks up unsold.
      W(econ::CommodityId::Food,       -0.12);
      W(econ::CommodityId::Water,      -0.10);
      W(econ::CommodityId::Fuel,       -0.18);
      W(econ::CommodityId::Medicine,   -0.15);
      W(econ::CommodityId::Metals,     -0.05);
      W(econ::CommodityId::Ore,        -0.05);
      W(econ::CommodityId::Electronics,+0.10);
      W(econ::CommodityId::Machinery,  +0.10);
      W(econ::CommodityId::Luxury,     +0.25);
      W(econ::CommodityId::Stimulants, +0.15);
      W(econ::CommodityId::Weapons,    +0.05);
      break;

    case SystemEventKind::PirateRaid:
      // Disrupted lanes + panic buying.
      W(econ::CommodityId::Food,       -0.25);
      W(econ::CommodityId::Water,      -0.25);
      W(econ::CommodityId::Fuel,       -0.40);
      W(econ::CommodityId::Medicine,   -0.35);
      W(econ::CommodityId::Weapons,    -0.45);
      W(econ::CommodityId::Luxury,     -0.30);
      W(econ::CommodityId::Electronics,-0.20);
      W(econ::CommodityId::Machinery,  -0.15);
      W(econ::CommodityId::Stimulants, -0.10);
      W(econ::CommodityId::Metals,     -0.10);
      W(econ::CommodityId::Ore,        -0.05);
      break;

    case SystemEventKind::SecurityCrackdown:
      // Contraband gets dumped/confiscated; security spending rises.
      W(econ::CommodityId::Electronics,-0.12);
      W(econ::CommodityId::Machinery,  -0.08);
      W(econ::CommodityId::Weapons,    +0.18);
      W(econ::CommodityId::Stimulants, +0.30);
      W(econ::CommodityId::Fuel,       -0.05);
      W(econ::CommodityId::Medicine,   -0.05);
      W(econ::CommodityId::Metals,     -0.05);
      W(econ::CommodityId::Luxury,     -0.05);
      break;

    case SystemEventKind::CivilUnrest:
      // Essentials run short; luxuries slump.
      W(econ::CommodityId::Food,       -0.55);
      W(econ::CommodityId::Water,      -0.55);
      W(econ::CommodityId::Medicine,   -0.45);
      W(econ::CommodityId::Fuel,       -0.35);
      W(econ::CommodityId::Weapons,    -0.25);
      W(econ::CommodityId::Stimulants, -0.10);
      W(econ::CommodityId::Luxury,     +0.20);
      W(econ::CommodityId::Electronics,+0.05);
      W(econ::CommodityId::Machinery,  +0.05);
      break;

    case SystemEventKind::ResearchBreakthrough:
      // Tech output surges; raw inputs get consumed.
      W(econ::CommodityId::Electronics,+0.35);
      W(econ::CommodityId::Machinery,  +0.25);
      W(econ::CommodityId::Medicine,   +0.20);
      W(econ::CommodityId::Metals,     -0.20);
      W(econ::CommodityId::Ore,        -0.15);
      W(econ::CommodityId::Luxury,     -0.05);
      W(econ::CommodityId::Stimulants, -0.10);
      W(econ::CommodityId::Fuel,       -0.05);
      break;

  }

  return w;
}

std::array<double, econ::kCommodityCount> systemEventExtraNetPerDay(
    const econ::StationEconomyModel& model,
    const SystemEvent& ev,
    const SystemEventEconomyParams& params) {

  std::array<double, econ::kCommodityCount> out{};
  out.fill(0.0);

  if (!ev.active) return out;

  const double sev = stellar::math::clamp(ev.severity01, 0.0, std::max(0.0, params.maxSeverity01));
  if (sev <= 1e-9) return out;

  const auto w = weightsFor(ev.kind);

  // Station-type modulation so events feel different between hubs/outposts.
  double stationMod = 1.0;
  switch (model.type) {
    case econ::StationType::Outpost:    stationMod = 0.85; break;
    case econ::StationType::TradeHub:   stationMod = 1.15; break;
    case econ::StationType::Shipyard:   stationMod = 1.10; break;
    case econ::StationType::Agricultural:
    case econ::StationType::Mining:
    case econ::StationType::Refinery:
    case econ::StationType::Industrial:
    case econ::StationType::Research:
    case econ::StationType::Count:
      stationMod = 1.0;
      break;
  }

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double wi = w[i];
    if (std::abs(wi) <= 1e-9) continue;

    const double desired = std::max(1e-9, model.desiredStock[i]);
    double flow = params.flowScaleFracOfDesiredPerDay * sev * desired * wi;

    // Make local producers/consumers amplify the relevant direction.
    double localMod = 1.0;
    if (flow > 0.0 && model.productionPerDay[i] > 0.0) localMod *= 1.25;
    if (flow < 0.0 && model.consumptionPerDay[i] > 0.0) localMod *= 1.25;

    flow *= stationMod * localMod;

    const double cap = std::max(0.0, params.maxAbsFlowFracOfDesiredPerDay) * desired;
    flow = stellar::math::clamp(flow, -cap, cap);

    out[i] = flow;
  }

  return out;
}

void applySystemEventInventoryShock(econ::StationEconomyState& state,
                                   const econ::StationEconomyModel& model,
                                   const SystemEvent& ev,
                                   const SystemEventEconomyParams& params) {
  if (!ev.active) return;

  const double sev = stellar::math::clamp(ev.severity01, 0.0, std::max(0.0, params.maxSeverity01));
  if (sev <= 1e-9) return;

  auto mulInv = [&](econ::CommodityId id, double mul) {
    const std::size_t i = idx(id);
    state.inventory[i] = std::max(0.0, state.inventory[i] * mul);
  };

  auto addInv = [&](econ::CommodityId id, double addUnits) {
    const std::size_t i = idx(id);
    state.inventory[i] = std::min(model.capacity[i], std::max(0.0, state.inventory[i] + addUnits));
  };

  switch (ev.kind) {
    case SystemEventKind::PirateRaid: {
      const double frac = stellar::math::clamp(params.pirateRaidTheftFrac, 0.0, 0.95) * sev;
      // High-value / strategic goods suffer most.
      mulInv(econ::CommodityId::Luxury,      1.0 - frac * 1.00);
      mulInv(econ::CommodityId::Weapons,     1.0 - frac * 0.85);
      mulInv(econ::CommodityId::Fuel,        1.0 - frac * 0.75);
      mulInv(econ::CommodityId::Medicine,    1.0 - frac * 0.70);
      mulInv(econ::CommodityId::Electronics, 1.0 - frac * 0.65);
      mulInv(econ::CommodityId::Food,        1.0 - frac * 0.55);
      mulInv(econ::CommodityId::Water,       1.0 - frac * 0.55);
      mulInv(econ::CommodityId::Stimulants,  1.0 - frac * 0.40);
      break;
    }

    case SystemEventKind::CivilUnrest: {
      const double frac = stellar::math::clamp(params.civilUnrestLossFrac, 0.0, 0.95) * sev;
      mulInv(econ::CommodityId::Food,      1.0 - frac * 1.00);
      mulInv(econ::CommodityId::Water,     1.0 - frac * 1.00);
      mulInv(econ::CommodityId::Medicine,  1.0 - frac * 0.80);
      mulInv(econ::CommodityId::Fuel,      1.0 - frac * 0.60);
      mulInv(econ::CommodityId::Weapons,   1.0 - frac * 0.35);
      mulInv(econ::CommodityId::Luxury,    1.0 - frac * 0.25);
      break;
    }

    case SystemEventKind::SecurityCrackdown: {
      const double frac = stellar::math::clamp(params.securityCrackdownDumpFrac, 0.0, 1.0) * sev;
      addInv(econ::CommodityId::Stimulants, frac * model.desiredStock[idx(econ::CommodityId::Stimulants)]);
      addInv(econ::CommodityId::Weapons,    frac * model.desiredStock[idx(econ::CommodityId::Weapons)]);
      break;
    }

    case SystemEventKind::ResearchBreakthrough: {
      const double frac = stellar::math::clamp(params.researchBreakthroughSurplusFrac, 0.0, 1.0) * sev;
      addInv(econ::CommodityId::Electronics, frac * model.desiredStock[idx(econ::CommodityId::Electronics)]);
      addInv(econ::CommodityId::Machinery,   frac * model.desiredStock[idx(econ::CommodityId::Machinery)]);
      addInv(econ::CommodityId::Medicine,    frac * model.desiredStock[idx(econ::CommodityId::Medicine)]);
      break;
    }

    case SystemEventKind::TradeBoom:
    case SystemEventKind::TradeBust:
    case SystemEventKind::None:
      break;
    default:
      break;
  }

  state.clampToCapacity(model);
}

std::string systemEventEconomySummary(const SystemEvent& ev, int maxUp, int maxDown) {
  if (!ev.active) return {};

  const double sev = stellar::math::clamp(ev.severity01, 0.0, 1.0);
  if (sev <= 1e-9) return {};

  const auto w = weightsFor(ev.kind);

  struct Item {
    double mag;
    econ::CommodityId id;
  };

  std::vector<Item> up;
  std::vector<Item> down;
  up.reserve(econ::kCommodityCount);
  down.reserve(econ::kCommodityCount);

  // Threshold scales gently with severity so mild events don't spam.
  const double thresh = 0.08 * (0.35 + 0.65 * sev);

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double wi = w[i];
    if (std::abs(wi) < thresh) continue;

    const auto id = static_cast<econ::CommodityId>(i);
    const double mag = std::abs(wi);

    // Negative weight => price pressure up (inventory down)
    if (wi < 0.0) up.push_back(Item{mag, id});
    else down.push_back(Item{mag, id});
  }

  auto byMagDesc = [](const Item& a, const Item& b) { return a.mag > b.mag; };
  std::sort(up.begin(), up.end(), byMagDesc);
  std::sort(down.begin(), down.end(), byMagDesc);

  maxUp = std::max(0, maxUp);
  maxDown = std::max(0, maxDown);

  if ((int)up.size() > maxUp) up.resize((std::size_t)maxUp);
  if ((int)down.size() > maxDown) down.resize((std::size_t)maxDown);

  if (up.empty() && down.empty()) {
    return "Market impact: (neutral)";
  }

  auto join = [](const std::vector<Item>& items, const char* suffix) {
    std::ostringstream oss;
    for (std::size_t i = 0; i < items.size(); ++i) {
      const auto& it = items[i];
      if (i) oss << ", ";
      oss << stellar::econ::commodityDef(it.id).name << suffix;
    }
    return oss.str();
  };

  std::ostringstream out;
  out << "Market impact: ";

  if (!up.empty()) {
    out << join(up, "+");
  } else {
    out << "(no major spikes)";
  }

  if (!down.empty()) {
    out << " | " << join(down, "-");
  }

  return out.str();
}

} // namespace stellar::sim
