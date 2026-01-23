#include "stellar/sim/StationServices.h"

#include "stellar/econ/Economy.h"
#include "stellar/econ/Commodity.h"

#include <cmath>
#include <iostream>

namespace {

static bool nearly(double a, double b, double eps = 1e-5) {
  return std::fabs(a - b) <= eps;
}

static std::size_t idx(stellar::econ::CommodityId id) {
  return static_cast<std::size_t>(id);
}

} // namespace

int test_station_services() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;
  using namespace stellar::econ;

  StationEconomyModel model = makeEconomyModel(StationType::Outpost, 0.0);
  StationEconomyState st{};

  // --------------------------- Hull repair (stock limited) ---------------------------
  {
    st = StationEconomyState{};
    st.inventory[idx(CommodityId::Metals)] = 10.0;
    st.inventory[idx(CommodityId::Machinery)] = 2.0; // should be the limiting reagent

    double credits = 1.0e9;
    double hull = 50.0;
    const double hullMax = 100.0;

    StationServicePriceModel pm{0.10, 0.10};
    const auto q = quoteHullRepairToFull(st, model, hull, hullMax, pm, credits);
    if (!q.ok) {
      std::cerr << "[station_services] expected repair quote ok\n";
      ++fails;
    } else {
      const double expected = std::min({hullMax - hull, 10.0 / kRepairMetalsPerHull, 2.0 / kRepairMachineryPerHull});
      if (!nearly(q.hullToRepair, expected, 1e-4)) {
        std::cerr << "[station_services] repair hullToRepair mismatch\n";
        ++fails;
      }

      const double creditsBefore = credits;
      const auto r = applyHullRepairToFull(st, model, credits, hull, hullMax, pm);
      if (!r.ok) {
        std::cerr << "[station_services] expected repair apply ok\n";
        ++fails;
      }
      if (!nearly(hull, 50.0 + expected, 1e-4)) {
        std::cerr << "[station_services] repaired hull did not increase as expected\n";
        ++fails;
      }
      if (!(credits < creditsBefore)) {
        std::cerr << "[station_services] credits not deducted for repair\n";
        ++fails;
      }
      if (!(st.inventory[idx(CommodityId::Machinery)] <= 1e-5)) {
        std::cerr << "[station_services] expected machinery to be depleted\n";
        ++fails;
      }
    }
  }

  // --------------------------- Hull repair (credit limited) --------------------------
  {
    st = StationEconomyState{};
    st.inventory[idx(CommodityId::Metals)] = 1000.0;
    st.inventory[idx(CommodityId::Machinery)] = 1000.0;

    double credits = 100.0;
    double hull = 25.0;
    const double hullMax = 100.0;

    StationServicePriceModel pm{0.10, 0.10};
    const auto q = quoteHullRepairToFull(st, model, hull, hullMax, pm, credits);
    if (!q.ok) {
      std::cerr << "[station_services] expected credit-limited repair quote ok\n";
      ++fails;
    }
    if (!q.limitedByCredits) {
      std::cerr << "[station_services] expected limitedByCredits for repair\n";
      ++fails;
    }

    const auto r = applyHullRepairToFull(st, model, credits, hull, hullMax, pm);
    if (!r.ok) {
      std::cerr << "[station_services] expected credit-limited repair apply ok\n";
      ++fails;
    }
    if (!(credits <= 100.0 + 1e-6)) {
      std::cerr << "[station_services] credits went up?\n";
      ++fails;
    }
    if (!(hull > 25.0)) {
      std::cerr << "[station_services] expected some hull repaired under credit limit\n";
      ++fails;
    }
  }

  // --------------------------- Refuel (stock limited) --------------------------------
  {
    st = StationEconomyState{};
    st.inventory[idx(CommodityId::Fuel)] = 10.0;

    double credits = 1.0e9;
    double fuel = 5.0;
    const double fuelMax = 25.0;

    StationServicePriceModel pm{0.10, 0.10};
    const auto q = quoteRefuelToFull(st, model, fuel, fuelMax, pm, credits);
    if (!q.ok) {
      std::cerr << "[station_services] expected refuel quote ok\n";
      ++fails;
    }
    if (!nearly(q.fuelToBuy, 10.0, 1e-6)) {
      std::cerr << "[station_services] expected to buy exactly station fuel inventory\n";
      ++fails;
    }

    const auto r = applyRefuelToFull(st, model, credits, fuel, fuelMax, pm);
    if (!r.ok) {
      std::cerr << "[station_services] expected refuel apply ok\n";
      ++fails;
    }
    if (!nearly(fuel, 15.0, 1e-6)) {
      std::cerr << "[station_services] refuel amount mismatch\n";
      ++fails;
    }
    if (!(st.inventory[idx(CommodityId::Fuel)] <= 1e-6)) {
      std::cerr << "[station_services] expected station fuel inventory depleted\n";
      ++fails;
    }
  }

  // --------------------------- Countermeasure restock (bounded) ----------------------
  {
    st = StationEconomyState{};
    st.inventory[idx(CommodityId::Weapons)] = 10.0;
    st.inventory[idx(CommodityId::Electronics)] = 10.0;
    st.inventory[idx(CommodityId::Metals)] = 10.0;
    st.inventory[idx(CommodityId::Machinery)] = 10.0;

    double credits = 5000.0;
    const ShipHullClass hullClass = ShipHullClass::Scout;
    int flares = 0;
    int chaff = 0;
    int heat = 0;

    StationServicePriceModel pm{0.10, 0.10};
    const auto q = quoteCountermeasureRestockAll(st, model, hullClass, flares, chaff, heat, pm, credits);
    if (!q.ok) {
      std::cerr << "[station_services] expected countermeasure quote ok\n";
      ++fails;
    }
    if (!(q.buyFlares >= 0 && q.buyChaff >= 0 && q.buyHeatSinks >= 0)) {
      std::cerr << "[station_services] negative countermeasure buy plan\n";
      ++fails;
    }

    const auto r = applyCountermeasureRestockAll(st, model, credits, hullClass, flares, chaff, heat, pm);
    if (!r.ok) {
      std::cerr << "[station_services] expected countermeasure apply ok\n";
      ++fails;
    }

    const auto cap = countermeasureCapsForHull(hullClass);
    if (!(flares <= cap.flares && chaff <= cap.chaff && heat <= cap.heatSinks)) {
      std::cerr << "[station_services] countermeasure counts exceeded caps\n";
      ++fails;
    }
  }

  // --------------------------- Ordnance rearm (stock + int ammo) ---------------------
  {
    st = StationEconomyState{};
    st.inventory[idx(CommodityId::Weapons)] = 5.2;
    st.inventory[idx(CommodityId::Electronics)] = 0.0;

    double credits = 1.0e9;
    const ShipHullClass hullClass = ShipHullClass::Scout;
    core::u8 ammo = 0;

    StationServicePriceModel pm{0.10, 0.10};
    const auto q = quoteOrdnanceRearmToFull(st, model, hullClass, WeaponType::HomingMissile, ammo, pm, credits);
    if (!q.ok) {
      std::cerr << "[station_services] expected ordnance quote ok\n";
      ++fails;
    }
    if (!(q.buyAmmo == 5)) {
      std::cerr << "[station_services] expected ordnance buyAmmo to floor station weapons inventory\n";
      ++fails;
    }

    const auto r = applyOrdnanceRearmToFull(st, model, credits, hullClass, WeaponType::HomingMissile, ammo, pm);
    if (!r.ok) {
      std::cerr << "[station_services] expected ordnance apply ok\n";
      ++fails;
    }
    if (!((int)ammo == 5)) {
      std::cerr << "[station_services] ammo not updated as expected\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_station_services] PASS\n";
  }
  return fails;
}
