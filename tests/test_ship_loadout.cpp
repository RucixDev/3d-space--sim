#include "stellar/sim/ShipLoadout.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_ship_loadout() {
  int fails = 0;

  // --- Sanity: table sizes and name pointers should look valid. ---
  {
    if (sim::hullDefCount() != 3) {
      std::cerr << "[test_ship_loadout] expected 3 hull defs, got=" << sim::hullDefCount() << "\n";
      ++fails;
    }
    if (sim::weaponDefCount() != 5) {
      std::cerr << "[test_ship_loadout] expected 5 weapon defs, got=" << sim::weaponDefCount() << "\n";
      ++fails;
    }
    for (std::size_t i = 0; i < sim::hullDefCount(); ++i) {
      if (!sim::kHullDefs[i].name || sim::kHullDefs[i].name[0] == '\0') {
        std::cerr << "[test_ship_loadout] hull def " << i << " has empty name\n";
        ++fails;
      }
    }
    for (std::size_t i = 0; i < sim::weaponDefCount(); ++i) {
      if (!sim::kWeaponDefs[i].name || sim::kWeaponDefs[i].name[0] == '\0') {
        std::cerr << "[test_ship_loadout] weapon def " << i << " has empty name\n";
        ++fails;
      }
    }
  }

  // --- Derived stats match the base hull tables at Mk1. ---
  {
    const auto ds = sim::computeShipDerivedStats(sim::ShipHullClass::Scout, 1, 1, 1);
    const auto& h = sim::hullDef(sim::ShipHullClass::Scout);

    if (!approx(ds.hullMax, h.hullMax, 1e-12)) {
      std::cerr << "[test_ship_loadout] Scout hullMax mismatch. got=" << ds.hullMax
                << " expected=" << h.hullMax << "\n";
      ++fails;
    }

    const double expectedShield = h.shieldBase * sim::kShields[1].mult;
    if (!approx(ds.shieldMax, expectedShield, 1e-12)) {
      std::cerr << "[test_ship_loadout] Scout shieldMax mismatch. got=" << ds.shieldMax
                << " expected=" << expectedShield << "\n";
      ++fails;
    }
  }

  // --- Monotonicity: higher Mk improves the intended axis. ---
  {
    const auto ds1 = sim::computeShipDerivedStats(sim::ShipHullClass::Scout, 1, 1, 1);
    const auto dsT3 = sim::computeShipDerivedStats(sim::ShipHullClass::Scout, 3, 1, 1);
    const auto dsS3 = sim::computeShipDerivedStats(sim::ShipHullClass::Scout, 1, 3, 1);
    const auto dsD3 = sim::computeShipDerivedStats(sim::ShipHullClass::Scout, 1, 1, 3);

    if (!(dsT3.baseLinAccelKmS2 > ds1.baseLinAccelKmS2 + 1e-12)) {
      std::cerr << "[test_ship_loadout] thruster Mk3 should increase baseLinAccel. got=" << dsT3.baseLinAccelKmS2
                << " base=" << ds1.baseLinAccelKmS2 << "\n";
      ++fails;
    }
    if (!(dsT3.baseAngAccelRadS2 > ds1.baseAngAccelRadS2 + 1e-12)) {
      std::cerr << "[test_ship_loadout] thruster Mk3 should increase baseAngAccel. got=" << dsT3.baseAngAccelRadS2
                << " base=" << ds1.baseAngAccelRadS2 << "\n";
      ++fails;
    }

    if (!(dsS3.shieldMax > ds1.shieldMax + 1e-12)) {
      std::cerr << "[test_ship_loadout] shield Mk3 should increase shieldMax. got=" << dsS3.shieldMax
                << " base=" << ds1.shieldMax << "\n";
      ++fails;
    }

    if (!(dsD3.heatCoolRate > ds1.heatCoolRate + 1e-12)) {
      std::cerr << "[test_ship_loadout] distributor Mk3 should increase heatCoolRate. got=" << dsD3.heatCoolRate
                << " base=" << ds1.heatCoolRate << "\n";
      ++fails;
    }

    if (!(dsD3.shieldRegenPerSimMin > ds1.shieldRegenPerSimMin + 1e-12)) {
      std::cerr << "[test_ship_loadout] distributor Mk3 should increase shield regen (indirectly). got="
                << dsD3.shieldRegenPerSimMin << " base=" << ds1.shieldRegenPerSimMin << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_ship_loadout] PASS\n";
  }
  return fails;
}
