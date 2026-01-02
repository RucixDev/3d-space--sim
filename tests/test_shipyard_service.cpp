#include "stellar/sim/ShipyardService.h"

#include "stellar/econ/Commodity.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {

static bool nearly(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

} // namespace

int test_shipyard_service() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;
  using namespace stellar::econ;

  const core::u64 seed = 0x12345678ull;

  Station st{};
  st.id = 42;
  st.name = "Arcadia Shipyard";
  st.type = StationType::Shipyard;
  st.factionId = 7;
  st.feeRate = 0.10;

  // --- Rep + trade-in multiplier ------------------------------------------------
  if (!nearly(tradeInMultiplier(-100.0), 0.50, 1e-6)) {
    std::cerr << "[shipyard_service] tradeInMultiplier(-100) unexpected\n";
    ++fails;
  }
  if (!nearly(tradeInMultiplier(100.0), 0.80, 1e-6)) {
    std::cerr << "[shipyard_service] tradeInMultiplier(100) unexpected\n";
    ++fails;
  }

  // --- Deterministic pricing ----------------------------------------------------
  {
    ShipyardPriceModel m{0.07, 12.0};
    const double p1 = shipyardPrice(seed, st.factionId, st.type, ShipyardItemTier::Basic, 2000.0, m);
    const double p2 = shipyardPrice(seed, st.factionId, st.type, ShipyardItemTier::Basic, 2000.0, m);
    if (!nearly(p1, p2, 1e-9)) {
      std::cerr << "[shipyard_service] shipyardPrice not deterministic\n";
      ++fails;
    }
    if (!(p1 > 0.0)) {
      std::cerr << "[shipyard_service] shipyardPrice returned non-positive\n";
      ++fails;
    }
  }

  // --- Cargo rack upgrade -------------------------------------------------------
  {
    ShipyardPriceModel m{0.05, 5.0};
    const auto q = quoteCargoRackUpgrade(seed, st, m);
    if (!q.ok) {
      std::cerr << "[shipyard_service] cargo rack quote not ok\n";
      ++fails;
    }
    if (!nearly(q.delta, 200.0, 1e-6)) {
      std::cerr << "[shipyard_service] cargo rack delta not 200kg\n";
      ++fails;
    }

    double credits = q.costCr + 10.0;
    double cap = 420.0;
    if (!applyCargoRackUpgrade(credits, cap, q)) {
      std::cerr << "[shipyard_service] applyCargoRackUpgrade failed unexpectedly\n";
      ++fails;
    }
    if (!nearly(cap, 620.0, 1e-6)) {
      std::cerr << "[shipyard_service] cargo capacity not increased correctly\n";
      ++fails;
    }
  }

  // --- Smuggle hold upgrade progression ----------------------------------------
  {
    ShipyardPriceModel m{0.00, 0.0};
    auto q1 = quoteSmuggleHoldUpgrade(seed, st, m, 0);
    auto q2 = quoteSmuggleHoldUpgrade(seed, st, m, 1);
    if (!(q1.ok && q2.ok)) {
      std::cerr << "[shipyard_service] smuggle hold quotes not ok\n";
      ++fails;
    }
    if (!(q2.costCr >= q1.costCr - 1e-6)) {
      std::cerr << "[shipyard_service] expected Mk2 smuggle upgrade to cost >= Mk1\n";
      ++fails;
    }
  }

  // --- Hull purchase: capacity scaling + cargo fit gating ----------------------
  {
    ShipyardPriceModel m{0.00, 0.0};
    const double cargoCap = 420.0;
    const double fuelMax = 45.0;

    auto q = quoteHullPurchase(seed,
                               st,
                               m,
                               ShipHullClass::Scout,
                               ShipHullClass::Hauler,
                               cargoCap,
                               fuelMax,
                               /*cargoMassKg=*/0.0);
    if (!q.ok) {
      std::cerr << "[shipyard_service] hull purchase quote not ok\n";
      ++fails;
    }
    if (!nearly(q.newCargoCapacityKg, cargoCap * (kHullDefs[(int)ShipHullClass::Hauler].cargoMult / kHullDefs[(int)ShipHullClass::Scout].cargoMult), 1e-6)) {
      std::cerr << "[shipyard_service] hull cargo capacity rescale mismatch\n";
      ++fails;
    }
    if (!nearly(q.newFuelMax, fuelMax * (kHullDefs[(int)ShipHullClass::Hauler].fuelMult / kHullDefs[(int)ShipHullClass::Scout].fuelMult), 1e-6)) {
      std::cerr << "[shipyard_service] hull fuel max rescale mismatch\n";
      ++fails;
    }

    auto qHeavy = quoteHullPurchase(seed,
                                    st,
                                    m,
                                    ShipHullClass::Scout,
                                    ShipHullClass::Hauler,
                                    cargoCap,
                                    fuelMax,
                                    /*cargoMassKg=*/9999.0);
    if (qHeavy.cargoFits) {
      std::cerr << "[shipyard_service] expected cargoFits=false for over-cap cargo\n";
      ++fails;
    }
  }

  // --- Trade-in purchase: rep should improve trade-in and lower cost ----------
  {
    ShipyardPriceModel low{0.05, -80.0};
    ShipyardPriceModel high{0.05, 80.0};
    const double oldBase = 9000.0;
    const double newBase = 18500.0;

    const auto qLow = quoteTradeInPurchase(seed,
                                          st,
                                          low,
                                          shipyardTierForPrice(oldBase),
                                          oldBase,
                                          shipyardTierForPrice(newBase),
                                          newBase);
    const auto qHigh = quoteTradeInPurchase(seed,
                                           st,
                                           high,
                                           shipyardTierForPrice(oldBase),
                                           oldBase,
                                           shipyardTierForPrice(newBase),
                                           newBase);
    if (!(qHigh.tradeInCr > qLow.tradeInCr)) {
      std::cerr << "[shipyard_service] expected higher rep to increase trade-in\n";
      ++fails;
    }
    if (!(qHigh.finalCostCr <= qLow.finalCostCr + 1e-6)) {
      std::cerr << "[shipyard_service] expected higher rep to not increase final cost\n";
      ++fails;
    }
  }

  // --- Weapon purchase: allow downgrading to free weapons without refund ------
  {
    ShipyardPriceModel m{0.10, 50.0};
    const auto q = quoteWeaponPurchase(seed, st, m, WeaponType::Railgun, WeaponType::Cannon);
    if (!q.ok) {
      std::cerr << "[shipyard_service] expected weapon purchase quote ok\n";
      ++fails;
    }
    if (!(q.finalCostCr <= 1e-6)) {
      std::cerr << "[shipyard_service] expected free downgrade to have zero final cost\n";
      ++fails;
    }
    if (!(q.tradeInCr > 0.0)) {
      std::cerr << "[shipyard_service] expected trade-in to be computed even if no refund\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_shipyard_service] PASS\n";
  }
  return fails;
}
