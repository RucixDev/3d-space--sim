#include "stellar/sim/Warehouse.h"

#include "stellar/sim/Reputation.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double cargoMassKg(const std::array<double, econ::kCommodityCount>& cargo) {
  double kg = 0.0;
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double u = cargo[i];
    if (!(u > 0.0)) continue;
    const auto id = static_cast<econ::CommodityId>(i);
    kg += u * econ::commodityDef(id).massKg;
  }
  return kg;
}

static double stationCapacityByTypeKg(econ::StationType t) {
  using econ::StationType;
  switch (t) {
    case StationType::Outpost:       return 5000.0;
    case StationType::Agricultural:  return 9000.0;
    case StationType::Mining:        return 9000.0;
    case StationType::Refinery:      return 11000.0;
    case StationType::Industrial:    return 14000.0;
    case StationType::Research:      return 8000.0;
    case StationType::TradeHub:      return 26000.0;
    case StationType::Shipyard:      return 20000.0;
    default:                         return 6000.0;
  }
}

static double stationFeeMultiplier(econ::StationType t) {
  using econ::StationType;
  // Higher multiplier => higher storage price.
  // This is a small gameplay flavor effect; the station tariff and reputation do the
  // heavy lifting.
  switch (t) {
    case StationType::Outpost:       return 1.15;
    case StationType::Agricultural:  return 0.95;
    case StationType::Mining:        return 1.05;
    case StationType::Refinery:      return 1.10;
    case StationType::Industrial:    return 1.08;
    case StationType::Research:      return 0.90;
    case StationType::TradeHub:      return 0.80;
    case StationType::Shipyard:      return 0.92;
    default:                         return 1.00;
  }
}

StationStorage* findStorage(std::vector<StationStorage>& storage, StationId stationId) {
  for (auto& e : storage) {
    if (e.stationId == stationId) return &e;
  }
  return nullptr;
}

const StationStorage* findStorage(const std::vector<StationStorage>& storage, StationId stationId) {
  for (const auto& e : storage) {
    if (e.stationId == stationId) return &e;
  }
  return nullptr;
}

StationStorage& getOrCreateStorage(std::vector<StationStorage>& storage,
                                  const Station& station,
                                  double timeDays) {
  if (auto* e = findStorage(storage, station.id)) return *e;

  StationStorage e{};
  e.stationId = station.id;
  e.stationType = station.type;
  e.factionId = station.factionId;
  e.feeRate = std::clamp(station.feeRate, 0.0, 1.0);
  e.lastFeeDay = timeDays;
  e.feesDueCr = 0.0;
  e.cargo.fill(0.0);

  storage.push_back(e);
  return storage.back();
}

double storageMassKg(const StationStorage& entry) {
  return cargoMassKg(entry.cargo);
}

double storageCapacityKg(const StationStorage& entry) {
  // Station capacity is purely gameplay; it is not part of station economy caps.
  return stationCapacityByTypeKg(entry.stationType);
}

double storageFeeRateCrPerKgPerDay(const StationStorage& entry, double rep) {
  // Base storage rate (credits per kg per day) before multipliers.
  // Tuned so that early-game storage is affordable but not free.
  constexpr double kBaseCrPerKgPerDay = 0.05;

  // Reputation modifies the station's tariff, which then influences storage costs.
  const double tariffEff = applyReputationToFeeRate(std::clamp(entry.feeRate, 0.0, 0.25), rep);

  // Convert tariff (0..0.25) into a modest multiplier (~0.8..1.675)
  const double tariffMul = 0.80 + 3.50 * tariffEff;

  const double typeMul = stationFeeMultiplier(entry.stationType);

  const double rate = kBaseCrPerKgPerDay * typeMul * tariffMul;
  return std::clamp(rate, 0.0, 2.0);
}

double estimateStorageDailyFeeCr(const StationStorage& entry, double rep) {
  const double massKg = storageMassKg(entry);
  if (!(massKg > 0.0)) return 0.0;

  const double rate = storageFeeRateCrPerKgPerDay(entry, rep);
  return massKg * rate;
}

void accrueStorageFees(StationStorage& entry, double timeDays, double rep) {
  // Handle time reset (load game) gracefully.
  if (!(timeDays > entry.lastFeeDay)) {
    entry.lastFeeDay = timeDays;
    return;
  }

  const double dtDays = timeDays - entry.lastFeeDay;
  if (!(dtDays > 0.0)) {
    entry.lastFeeDay = timeDays;
    return;
  }

  const double daily = estimateStorageDailyFeeCr(entry, rep);
  const double add = daily * dtDays;
  if (add > 0.0 && std::isfinite(add)) {
    entry.feesDueCr += add;
  }
  entry.lastFeeDay = timeDays;
}

WarehouseResult payStorageFees(std::vector<StationStorage>& storage,
                              const Station& station,
                              double timeDays,
                              double rep,
                              double& credits,
                              double amountCr) {
  WarehouseResult r{};

  if (!(amountCr > 0.0) || !std::isfinite(amountCr)) {
    r.ok = false;
    r.reason = "invalid";
    return r;
  }

  StationStorage* e = findStorage(storage, station.id);
  if (!e) {
    r.ok = false;
    r.reason = "no_storage";
    return r;
  }

  accrueStorageFees(*e, timeDays, rep);

  const double due = std::max(0.0, e->feesDueCr);
  if (!(due > 1e-6)) {
    r.ok = true;
    r.unitsMoved = 0.0;
    r.creditsPaid = 0.0;
    return r;
  }

  const double pay = std::max(0.0, std::min({amountCr, due, credits}));
  if (!(pay > 1e-6)) {
    r.ok = false;
    r.reason = "invalid";
    return r;
  }

  credits -= pay;
  e->feesDueCr = std::max(0.0, due - pay);
  if (e->feesDueCr <= 1e-3) e->feesDueCr = 0.0;

  r.ok = true;
  r.creditsPaid = pay;
  return r;
}

WarehouseResult depositToStorage(std::vector<StationStorage>& storage,
                                const Station& station,
                                double timeDays,
                                double rep,
                                std::array<double, econ::kCommodityCount>& shipCargo,
                                econ::CommodityId commodity,
                                double units) {
  WarehouseResult r{};

  const std::size_t ci = static_cast<std::size_t>(commodity);
  if (ci >= econ::kCommodityCount || !(units > 0.0) || !std::isfinite(units)) {
    r.ok = false;
    r.reason = "invalid";
    return r;
  }

  const double have = std::max(0.0, shipCargo[ci]);
  if (!(have > 1e-9)) {
    r.ok = false;
    r.reason = "empty";
    return r;
  }

  StationStorage& e = getOrCreateStorage(storage, station, timeDays);
  accrueStorageFees(e, timeDays, rep);

  const double want = std::min(units, have);

  const double massPerUnit = econ::commodityDef(commodity).massKg;
  const double capKg = storageCapacityKg(e);
  const double massNowKg = storageMassKg(e);
  const double freeKg = std::max(0.0, capKg - massNowKg);

  if (!(freeKg > 1e-9)) {
    r.ok = false;
    r.reason = "capacity";
    return r;
  }

  double maxUnitsByMass = want;
  if (massPerUnit > 1e-9) {
    maxUnitsByMass = std::min(want, freeKg / massPerUnit);
  }

  const double moved = std::max(0.0, maxUnitsByMass);
  if (!(moved > 1e-9)) {
    r.ok = false;
    r.reason = "capacity";
    return r;
  }

  shipCargo[ci] = std::max(0.0, shipCargo[ci] - moved);
  e.cargo[ci] += moved;

  r.ok = true;
  r.unitsMoved = moved;
  return r;
}

WarehouseResult withdrawFromStorage(std::vector<StationStorage>& storage,
                                  const Station& station,
                                  double timeDays,
                                  double rep,
                                  std::array<double, econ::kCommodityCount>& shipCargo,
                                  double& credits,
                                  double cargoCapacityKg,
                                  econ::CommodityId commodity,
                                  double units) {
  WarehouseResult r{};

  const std::size_t ci = static_cast<std::size_t>(commodity);
  if (ci >= econ::kCommodityCount || !(units > 0.0) || !std::isfinite(units)) {
    r.ok = false;
    r.reason = "invalid";
    return r;
  }

  StationStorage* e = findStorage(storage, station.id);
  if (!e) {
    r.ok = false;
    r.reason = "no_storage";
    return r;
  }

  accrueStorageFees(*e, timeDays, rep);

  // Withdrawal is blocked by unpaid fees. For convenience, we auto-pay if possible.
  if (e->feesDueCr > 1e-3) {
    if (credits + 1e-6 < e->feesDueCr) {
      r.ok = false;
      r.reason = "fees_due";
      return r;
    }

    credits -= e->feesDueCr;
    r.creditsPaid = e->feesDueCr;
    e->feesDueCr = 0.0;
  }

  const double stored = std::max(0.0, e->cargo[ci]);
  if (!(stored > 1e-9)) {
    r.ok = false;
    r.reason = "empty";
    return r;
  }

  const double massPerUnit = econ::commodityDef(commodity).massKg;
  const double shipMassKg = cargoMassKg(shipCargo);
  const double freeKg = std::max(0.0, cargoCapacityKg - shipMassKg);

  if (!(freeKg > 1e-9)) {
    r.ok = false;
    r.reason = "cargo_full";
    return r;
  }

  double maxUnitsByMass = units;
  if (massPerUnit > 1e-9) {
    maxUnitsByMass = std::min(units, freeKg / massPerUnit);
  }

  const double moved = std::max(0.0, std::min(stored, maxUnitsByMass));
  if (!(moved > 1e-9)) {
    r.ok = false;
    r.reason = "cargo_full";
    return r;
  }

  e->cargo[ci] = std::max(0.0, e->cargo[ci] - moved);
  shipCargo[ci] += moved;

  r.ok = true;
  r.unitsMoved = moved;
  return r;
}

void pruneEmptyStorage(std::vector<StationStorage>& storage, double epsUnits, double epsFees) {
  storage.erase(
    std::remove_if(storage.begin(), storage.end(), [&](const StationStorage& e) {
      if (e.feesDueCr > epsFees) return false;
      for (double u : e.cargo) {
        if (u > epsUnits) return false;
      }
      return true;
    }),
    storage.end());
}

} // namespace stellar::sim
