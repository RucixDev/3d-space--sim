#include "stellar/sim/Insurance.h"

#include <cmath>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_insurance() {
  using namespace stellar::sim;

  int fails = 0;

  {
    // Starter Scout with no upgrades should still have a non-zero rebuy
    // due to the minimum.
    PlayerShipEconomyState st{};
    const auto val = computeShipValue(st);
    if (!nearly(val.totalCr, 0.0)) {
      std::cerr << "[test_insurance] expected starter ship value ~0, got " << val.totalCr << "\n";
      ++fails;
    }

    InsurancePolicy pol{};
    const double rebuy = computeRebuyCost(pol, val.totalCr);
    if (!nearly(rebuy, pol.minRebuyCr)) {
      std::cerr << "[test_insurance] expected min rebuy " << pol.minRebuyCr << ", got " << rebuy << "\n";
      ++fails;
    }
  }

  {
    // Cargo/fuel upgrade inference should be stable across hull scaling.
    PlayerShipEconomyState st{};
    st.hull = ShipHullClass::Hauler;

    // One cargo rack upgrade in the "Scout frame": 420 + 200 = 620, then scaled by hauler mult.
    const HullDef& hd = hullDef(st.hull);
    st.cargoCapacityKg = (420.0 + 200.0) * hd.cargoMult;

    // Two fuel tank upgrades: 45 + 20 = 65, scaled.
    st.fuelMax = (45.0 + 20.0) * hd.fuelMult;

    // One passenger upgrade, one FSD tuning.
    st.passengerSeats = 4;
    st.fsdRangeLy = 20.0;

    // Smuggle mk2 => 4500 + 8500
    st.smuggleHoldMk = 2;

    const auto val = computeShipValue(st);
    const double expectedUpgrades = 2000.0 + 2.0 * 2500.0 + 3000.0 + 9000.0 + (4500.0 + 8500.0);
    if (!nearly(val.upgradesCr, expectedUpgrades)) {
      std::cerr << "[test_insurance] upgrades mismatch: expected " << expectedUpgrades
                << " got " << val.upgradesCr << "\n";
      ++fails;
    }
  }

  {
    // Rebuy path: paid from credits.
    InsurancePolicy pol{};
    pol.rebuyRate = 0.10;
    pol.minRebuyCr = 0.0;

    PlayerShipEconomyState st{};
    st.hull = ShipHullClass::Fighter;
    st.thrusterMk = 3;
    st.shieldMk = 3;
    st.distributorMk = 3;
    st.weaponPrimary = WeaponType::Railgun;
    st.weaponSecondary = WeaponType::HomingMissile;

    const double shipVal = computeShipValue(st).totalCr;
    const double rebuy = computeRebuyCost(pol, shipVal);

    double credits = rebuy + 123.0;
    double debt = 0.0;
    auto out = applyRebuyOnDeath(pol, st, credits, debt);
    if (out.result != RebuyResult::Paid) {
      std::cerr << "[test_insurance] expected Paid\n";
      ++fails;
    }
    if (!nearly(credits, 123.0, 1e-6)) {
      std::cerr << "[test_insurance] credits mismatch after Paid rebuy\n";
      ++fails;
    }
    if (!nearly(debt, 0.0)) {
      std::cerr << "[test_insurance] debt should remain 0 after Paid rebuy\n";
      ++fails;
    }
  }

  {
    // Rebuy path: loan.
    InsurancePolicy pol{};
    pol.rebuyRate = 0.10;
    pol.minRebuyCr = 0.0;
    pol.loanMaxCr = 20000.0;

    PlayerShipEconomyState st{};
    st.hull = ShipHullClass::Fighter;
    st.thrusterMk = 3;
    st.shieldMk = 3;
    st.distributorMk = 3;
    st.weaponPrimary = WeaponType::Railgun;
    st.weaponSecondary = WeaponType::HomingMissile;

    const double rebuy = computeRebuyCost(pol, computeShipValue(st).totalCr);
    double credits = rebuy * 0.25;
    double debt = 0.0;

    auto out = applyRebuyOnDeath(pol, st, credits, debt);
    if (out.result != RebuyResult::Loan) {
      std::cerr << "[test_insurance] expected Loan\n";
      ++fails;
    }
    if (!nearly(credits, 0.0)) {
      std::cerr << "[test_insurance] credits should go to 0 when loaning\n";
      ++fails;
    }
    if (!(debt > 0.0)) {
      std::cerr << "[test_insurance] debt should increase when loaning\n";
      ++fails;
    }
  }

  {
    // Rebuy path: bankrupt reset (insufficient remaining loan headroom).
    InsurancePolicy pol{};
    pol.rebuyRate = 0.10;
    pol.minRebuyCr = 0.0;
    pol.loanMaxCr = 1000.0;

    PlayerShipEconomyState st{};
    st.hull = ShipHullClass::Fighter;
    st.thrusterMk = 3;
    st.shieldMk = 3;
    st.distributorMk = 3;
    st.weaponPrimary = WeaponType::Railgun;
    st.weaponSecondary = WeaponType::HomingMissile;

    const double rebuy = computeRebuyCost(pol, computeShipValue(st).totalCr);
    double credits = 0.0;
    double debt = 900.0; // leaving only 100 loan headroom

    auto out = applyRebuyOnDeath(pol, st, credits, debt);
    if (out.result != RebuyResult::BankruptReset || !out.shipReset) {
      std::cerr << "[test_insurance] expected BankruptReset\n";
      ++fails;
    }
    if (st.hull != ShipHullClass::Scout || st.thrusterMk != 1) {
      std::cerr << "[test_insurance] expected loaner Scout loadout\n";
      ++fails;
    }
    if (!nearly(debt, 0.0)) {
      std::cerr << "[test_insurance] expected debt cleared on bankruptcy\n";
      ++fails;
    }
    (void)rebuy;
  }

  return fails;
}
