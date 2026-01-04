#include "stellar/sim/Insurance.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

namespace {

static int clampInt(int v, int lo, int hi) {
  return std::max(lo, std::min(v, hi));
}

static int upgradeCount(double v, double base, double step) {
  if (step <= 1e-9) return 0;
  const double d = (v - base) / step;
  // Use rounding rather than floor to tolerate FP noise from hull scaling.
  const long long n = (long long)std::llround(d);
  return (int)std::max<long long>(0, std::min<long long>(n, 100000));
}

static double smuggleHoldValueCr(int mk) {
  // Base prices used by ShipyardService (see quoteSmuggleHoldUpgrade).
  static constexpr double kBaseCost[4] = {0.0, 4500.0, 8500.0, 14000.0};
  const int m = clampInt(mk, 0, 3);
  double total = 0.0;
  for (int i = 1; i <= m; ++i) total += kBaseCost[i];
  return total;
}

} // namespace

ShipValueBreakdown computeShipValue(const PlayerShipEconomyState& s) {
  ShipValueBreakdown out{};

  const HullDef& hd = hullDef(s.hull);
  out.hullCr = std::max(0.0, hd.priceCr);

  const int tMk = clampInt(s.thrusterMk, 1, 3);
  const int shMk = clampInt(s.shieldMk, 1, 3);
  const int dMk = clampInt(s.distributorMk, 1, 3);
  out.coreModulesCr = std::max(0.0, kThrusters[tMk].priceCr)
                    + std::max(0.0, kShields[shMk].priceCr)
                    + std::max(0.0, kDistributors[dMk].priceCr);

  out.weaponsCr = std::max(0.0, weaponDef(s.weaponPrimary).priceCr)
                + std::max(0.0, weaponDef(s.weaponSecondary).priceCr);

  // Shipyard upgrades. These are stored in SaveGame as the final tuned values.
  // We infer how many upgrades were applied by converting back into the
  // "Scout frame" using the current hull's multipliers.
  const double cargoBaseFrame = (hd.cargoMult > 1e-9) ? (s.cargoCapacityKg / hd.cargoMult) : s.cargoCapacityKg;
  const double fuelBaseFrame  = (hd.fuelMult  > 1e-9) ? (s.fuelMax / hd.fuelMult) : s.fuelMax;

  const int cargoUp = upgradeCount(cargoBaseFrame, 420.0, 200.0);
  const int fuelUp  = upgradeCount(fuelBaseFrame, 45.0, 10.0);
  const int paxUp   = upgradeCount((double)s.passengerSeats, 2.0, 2.0);
  const int fsdUp   = upgradeCount(s.fsdRangeLy, 18.0, 2.0);

  const double cargoValue = (double)cargoUp * 2000.0;
  const double fuelValue  = (double)fuelUp * 2500.0;
  const double paxValue   = (double)paxUp * 3000.0;
  const double fsdValue   = (double)fsdUp * 9000.0;
  const double smuggleVal = smuggleHoldValueCr(s.smuggleHoldMk);

  out.upgradesCr = cargoValue + fuelValue + paxValue + fsdValue + smuggleVal;

  out.totalCr = out.hullCr + out.coreModulesCr + out.weaponsCr + out.upgradesCr;
  return out;
}

double computeRebuyCost(const InsurancePolicy& policy, double shipValueCr) {
  const double v = std::max(0.0, shipValueCr);
  const double rate = std::max(0.0, policy.rebuyRate);
  const double minC = std::max(0.0, policy.minRebuyCr);
  return std::max(minC, v * rate);
}

RebuyOutcome applyRebuyOnDeath(const InsurancePolicy& policy,
                              PlayerShipEconomyState& ship,
                              double& credits,
                              double& debtCr) {
  RebuyOutcome out{};

  const ShipValueBreakdown val = computeShipValue(ship);
  out.shipValueCr = val.totalCr;
  out.rebuyCostCr = computeRebuyCost(policy, val.totalCr);

  credits = std::max(0.0, credits);
  debtCr = std::max(0.0, debtCr);

  const double rebuy = out.rebuyCostCr;
  const double loanHeadroom = std::max(0.0, policy.loanMaxCr - debtCr);

  if (credits + 1e-6 >= rebuy) {
    out.result = RebuyResult::Paid;
    out.paidFromCreditsCr = rebuy;
    credits -= rebuy;
    out.debtAfterCr = debtCr;
    return out;
  }

  if (credits + loanHeadroom + 1e-6 >= rebuy) {
    out.result = RebuyResult::Loan;
    out.paidFromCreditsCr = credits;
    const double shortfall = std::max(0.0, rebuy - credits);
    credits = 0.0;
    debtCr += shortfall;
    out.loanTakenCr = shortfall;
    out.debtAfterCr = debtCr;
    return out;
  }

  // Bankruptcy: can't cover rebuy even with loan headroom.
  out.result = RebuyResult::BankruptReset;
  out.shipReset = true;

  credits = std::min(std::max(0.0, credits), std::max(0.0, policy.bankruptcyKeepCreditsCr));
  debtCr = 0.0;

  ship = PlayerShipEconomyState{}; // loaner Scout
  out.debtAfterCr = debtCr;
  return out;
}

} // namespace stellar::sim
