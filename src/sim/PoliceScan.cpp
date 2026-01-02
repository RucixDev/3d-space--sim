#include "stellar/sim/PoliceScan.h"

#include "stellar/econ/Commodity.h"
#include "stellar/sim/Contraband.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static int clampHoldMk(int mk) {
  if (mk < 0) return 0;
  if (mk > 3) return 3;
  return mk;
}

double smuggleHoldScanRateMult(int smuggleHoldMk) {
  switch (clampHoldMk(smuggleHoldMk)) {
    case 0: return 1.0;
    case 1: return 0.72;
    case 2: return 0.50;
    default: return 0.35;
  }
}

double smuggleHoldScanDurationMult(int smuggleHoldMk) {
  switch (clampHoldMk(smuggleHoldMk)) {
    case 0: return 1.0;
    case 1: return 1.15;
    case 2: return 1.35;
    default: return 1.60;
  }
}

double cargoScanStartRatePerSec(bool hasContraband, const LawProfile& law, int smuggleHoldMk) {
  const double strict = std::clamp(law.scanStrictness, 0.5, 2.0);
  const double baseRatePerSec = (hasContraband ? 0.06 : 0.012) * strict;
  return baseRatePerSec * smuggleHoldScanRateMult(smuggleHoldMk);
}

double cargoScanDurationSecStation(bool hasContraband, int smuggleHoldMk) {
  const double baseDur = hasContraband ? 7.0 : 4.0;
  return baseDur * smuggleHoldScanDurationMult(smuggleHoldMk);
}

double cargoScanDurationSecPolice(int smuggleHoldMk) {
  const double baseDur = 6.0;
  return baseDur * smuggleHoldScanDurationMult(smuggleHoldMk);
}

IllegalCargoScanResult scanIllegalCargoMask(core::u32 illegalMask,
                                           const std::array<double, econ::kCommodityCount>& cargoUnits,
                                           const std::array<double, econ::kCommodityCount>* midPriceOverrideCr) {
  IllegalCargoScanResult r{};
  r.scannedIllegalUnits.fill(0.0);

  const core::u32 mask = illegalMask;
  if (mask == 0u) return r;

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    if ((mask & ((core::u32)1u << (core::u32)i)) == 0u) continue;

    const double units = std::max(0.0, cargoUnits[i]);
    if (units <= 1e-9) continue;

    const auto cid = static_cast<econ::CommodityId>(i);
    double mid = econ::commodityDef(cid).basePrice;
    if (midPriceOverrideCr) {
      mid = std::max(0.0, (*midPriceOverrideCr)[i]);
      if (mid <= 1e-9) mid = econ::commodityDef(cid).basePrice;
    }

    r.illegalValueCr += units * mid;
    r.scannedIllegalUnits[i] = units;
  }

  return r;
}

IllegalCargoScanResult scanIllegalCargo(core::u64 universeSeed,
                                       core::u32 factionId,
                                       const std::array<double, econ::kCommodityCount>& cargoUnits,
                                       const std::array<double, econ::kCommodityCount>* midPriceOverrideCr) {
  if (factionId == 0) return IllegalCargoScanResult{};

  const core::u32 mask = illegalCommodityMask(universeSeed, factionId);
  return scanIllegalCargoMask(mask, cargoUnits, midPriceOverrideCr);
}

double bribeOfferChance(const LawProfile& law, double playerRep, double playerHeat, double illegalValueCr) {
  double chance = 0.33;
  chance *= std::clamp(law.corruption, 0.0, 1.0);

  // Lower willingness to take bribes when the ship is already running hot.
  chance *= std::clamp(1.10 - (playerHeat / 110.0), 0.25, 1.10);

  // Reputation influences willingness.
  if (playerRep < -30.0) chance *= 0.35;
  if (playerRep > 55.0) chance *= 0.75;

  // Huge hauls are riskier to ignore.
  chance *= std::clamp(1.05 - (std::max(0.0, illegalValueCr) / 20000.0), 0.35, 1.05);

  return std::clamp(chance, 0.0, 1.0);
}

double bribeAmountCr(const LawProfile& law, double illegalValueCr) {
  return std::max(0.0, law.contrabandBribeCr(std::max(0.0, illegalValueCr)));
}

double bribeAmountCrRounded(const LawProfile& law, double illegalValueCr, double roundToCr) {
  double a = bribeAmountCr(law, illegalValueCr);
  if (roundToCr > 1e-9) {
    a = std::round(a / roundToCr) * roundToCr;
  }
  return a;
}

ContrabandEnforcementResult enforceContraband(const LawProfile& law,
                                             double credits,
                                             const std::array<double, econ::kCommodityCount>& cargoUnits,
                                             const std::array<double, econ::kCommodityCount>& scannedIllegalUnits,
                                             double illegalValueCr) {
  ContrabandEnforcementResult r{};
  r.cargoAfter = cargoUnits;
  r.confiscatedUnits.fill(0.0);

  // Confiscate the scanned illegal goods (or whatever remains of them).
  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    const double seen = scannedIllegalUnits[i];
    if (seen <= 1e-9) continue;

    const double have = std::max(0.0, r.cargoAfter[i]);
    const double take = std::min(have, std::max(0.0, seen));

    if (take > 1e-9) {
      r.cargoAfter[i] = std::max(0.0, have - take);
      r.confiscatedUnits[i] = take;
    }
  }

  const double v = std::max(0.0, illegalValueCr);

  r.fineCr = law.contrabandFineCr(v);
  r.paidCr = std::min(std::max(0.0, credits), r.fineCr);
  r.unpaidCr = std::max(0.0, r.fineCr - r.paidCr);

  r.creditsAfter = std::max(0.0, credits) - r.paidCr;

  r.repPenalty = law.contrabandRepPenalty(v);
  r.bountyAddedCr = r.unpaidCr;

  // Prototype escalation tuning.
  r.policeAlertSeconds = 90.0;
  r.nextPoliceSpawnDelaySeconds = 6.0;

  const double strictMult = std::clamp(0.85 + 0.15 * law.scanStrictness, 0.75, 1.15);
  r.policeHeatDelta = std::min(1.8, (0.5 + v / 5000.0) * strictMult);
  if (r.policeHeatDelta < 0.0) r.policeHeatDelta = 0.0;

  return r;
}

ContrabandEvadeResult evadeContraband(const LawProfile& law,
                                     double fineCr,
                                     double illegalValueCr,
                                     bool escalateLocal) {
  ContrabandEvadeResult r{};

  const double v = std::max(0.0, illegalValueCr);
  const double fine = std::max(0.0, fineCr);

  r.repPenalty = law.contrabandEvadeRepPenalty(v);
  r.bountyAddedCr = fine;
  r.shipHeatDelta = 5.0;

  if (escalateLocal) {
    const double strictMult = std::clamp(0.90 + 0.20 * law.scanStrictness, 0.80, 1.25);
    r.policeHeatDelta = (1.2 + std::min(2.0, fine / 2500.0)) * strictMult;
    r.policeAlertSeconds = 140.0;
    r.nextPoliceSpawnDelaySeconds = 4.0;
  }

  return r;
}

} // namespace stellar::sim
