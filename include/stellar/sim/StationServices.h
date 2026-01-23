#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Market.h"
#include "stellar/sim/CountermeasureLoadout.h"
#include "stellar/sim/ShipLoadout.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// StationServices (headless)
// -----------------------------------------------------------------------------
//
// The prototype historically implemented "repair/refuel/restock" directly in the
// giant renderer app (apps/stellar_game/main.cpp). This module extracts the
// deterministic *business logic* into the core library so:
//   - the same rules are shared by UI, tools and tests
//   - services can be backed by the station economy inventory
//   - scarcity becomes gameplay (remote outposts can run out of parts/ammo)


struct StationServicePriceModel {
  // Effective station fee rate (0..1). The game typically computes this as
  // applyRepToFee(st.feeRate, rep).
  double feeRateEff{0.0};

  // Bid/ask spread used when pricing station inventory consumption.
  double bidAskSpread{0.10};
};

// -----------------------------------------------------------------------------
// Hull repair
// -----------------------------------------------------------------------------

// Repair material recipe tuning.
// Units are in station commodity "units" per hull point.
inline constexpr double kRepairMetalsPerHull = 0.20;
inline constexpr double kRepairMachineryPerHull = 0.14;

struct HullRepairQuote {
  bool ok{false};
  const char* reason{nullptr};

  double hullMissing{0.0};
  double hullToRepair{0.0};

  // Planned consumption.
  double metalsUnits{0.0};
  double machineryUnits{0.0};

  // Pricing snapshot.
  double metalsAsk{0.0};
  double machineryAsk{0.0};

  double costCr{0.0};
  bool limitedByStock{false};
  bool limitedByCredits{false};
};

HullRepairQuote quoteHullRepairToFull(const econ::StationEconomyState& stEcon,
                                     const econ::StationEconomyModel& model,
                                     double hullCurrent,
                                     double hullMax,
                                     const StationServicePriceModel& pm,
                                     double creditsBudgetCr = -1.0);

struct HullRepairResult {
  bool ok{false};
  const char* reason{nullptr};

  double hullRepaired{0.0};
  double creditsPaid{0.0};
  double metalsTaken{0.0};
  double machineryTaken{0.0};
};

HullRepairResult applyHullRepairToFull(econ::StationEconomyState& stEcon,
                                      const econ::StationEconomyModel& model,
                                      double& ioCredits,
                                      double& ioHullCurrent,
                                      double hullMax,
                                      const StationServicePriceModel& pm);

// -----------------------------------------------------------------------------
// Refuel
// -----------------------------------------------------------------------------

struct RefuelQuote {
  bool ok{false};
  const char* reason{nullptr};

  double fuelMissing{0.0};
  double fuelToBuy{0.0};
  double fuelAsk{0.0};
  double fuelInventory{0.0};

  double costCr{0.0};
  bool limitedByStock{false};
  bool limitedByCredits{false};
};

RefuelQuote quoteRefuelToFull(const econ::StationEconomyState& stEcon,
                             const econ::StationEconomyModel& model,
                             double fuelCurrent,
                             double fuelMax,
                             const StationServicePriceModel& pm,
                             double creditsBudgetCr = -1.0);

struct RefuelResult {
  bool ok{false};
  const char* reason{nullptr};

  double fuelBought{0.0};
  double creditsPaid{0.0};
};

RefuelResult applyRefuelToFull(econ::StationEconomyState& stEcon,
                              const econ::StationEconomyModel& model,
                              double& ioCredits,
                              double& ioFuelCurrent,
                              double fuelMax,
                              const StationServicePriceModel& pm);

// -----------------------------------------------------------------------------
// Countermeasure restock
// -----------------------------------------------------------------------------

struct CountermeasureRestockQuote {
  bool ok{false};
  const char* reason{nullptr};

  CountermeasureCaps caps{};

  int haveFlares{0};
  int haveChaff{0};
  int haveHeatSinks{0};

  int needFlares{0};
  int needChaff{0};
  int needHeatSinks{0};

  int buyFlares{0};
  int buyChaff{0};
  int buyHeatSinks{0};

  // Commodity usage plan.
  double useWeapons{0.0};
  double useElectronics{0.0};
  double useMetals{0.0};
  double useMachinery{0.0};

  double costCr{0.0};
  bool limitedByStock{false};
  bool limitedByCredits{false};
};

// Quote a "best effort" restock plan toward the hull's caps.
//
// - If creditsBudgetCr < 0, assumes infinite credits (inventory-limited plan).
// - Uses a small brute-force search (caps are tiny) to maximize defensive utility
//   under commodity and (optional) credit constraints.
CountermeasureRestockQuote quoteCountermeasureRestockAll(const econ::StationEconomyState& stEcon,
                                                        const econ::StationEconomyModel& model,
                                                        ShipHullClass hullClass,
                                                        int haveFlares,
                                                        int haveChaff,
                                                        int haveHeatSinks,
                                                        const StationServicePriceModel& pm,
                                                        double creditsBudgetCr = -1.0);

CountermeasureRestockQuote quoteCountermeasureRestockFlaresToCap(const econ::StationEconomyState& stEcon,
                                                                 const econ::StationEconomyModel& model,
                                                                 ShipHullClass hullClass,
                                                                 int haveFlares,
                                                                 const StationServicePriceModel& pm,
                                                                 double creditsBudgetCr = -1.0);

CountermeasureRestockQuote quoteCountermeasureRestockChaffToCap(const econ::StationEconomyState& stEcon,
                                                                const econ::StationEconomyModel& model,
                                                                ShipHullClass hullClass,
                                                                int haveChaff,
                                                                const StationServicePriceModel& pm,
                                                                double creditsBudgetCr = -1.0);

CountermeasureRestockQuote quoteCountermeasureRestockHeatSinksToCap(const econ::StationEconomyState& stEcon,
                                                                    const econ::StationEconomyModel& model,
                                                                    ShipHullClass hullClass,
                                                                    int haveHeatSinks,
                                                                    const StationServicePriceModel& pm,
                                                                    double creditsBudgetCr = -1.0);

struct CountermeasureRestockResult {
  bool ok{false};
  const char* reason{nullptr};

  int flaresBought{0};
  int chaffBought{0};
  int heatSinksBought{0};

  double creditsPaid{0.0};

  // Actual commodity consumption.
  double weaponsTaken{0.0};
  double electronicsTaken{0.0};
  double metalsTaken{0.0};
  double machineryTaken{0.0};
};

CountermeasureRestockResult applyCountermeasureRestockAll(econ::StationEconomyState& stEcon,
                                                         const econ::StationEconomyModel& model,
                                                         double& ioCredits,
                                                         ShipHullClass hullClass,
                                                         int& ioFlares,
                                                         int& ioChaff,
                                                         int& ioHeatSinks,
                                                         const StationServicePriceModel& pm);

CountermeasureRestockResult applyCountermeasureRestockFlaresToCap(econ::StationEconomyState& stEcon,
                                                                  const econ::StationEconomyModel& model,
                                                                  double& ioCredits,
                                                                  ShipHullClass hullClass,
                                                                  int& ioFlares,
                                                                  const StationServicePriceModel& pm);

CountermeasureRestockResult applyCountermeasureRestockChaffToCap(econ::StationEconomyState& stEcon,
                                                                 const econ::StationEconomyModel& model,
                                                                 double& ioCredits,
                                                                 ShipHullClass hullClass,
                                                                 int& ioChaff,
                                                                 const StationServicePriceModel& pm);

CountermeasureRestockResult applyCountermeasureRestockHeatSinksToCap(econ::StationEconomyState& stEcon,
                                                                     const econ::StationEconomyModel& model,
                                                                     double& ioCredits,
                                                                     ShipHullClass hullClass,
                                                                     int& ioHeatSinks,
                                                                     const StationServicePriceModel& pm);

// -----------------------------------------------------------------------------
// Ordnance / missile rearm
// -----------------------------------------------------------------------------

struct OrdnanceRearmQuote {
  bool ok{false};
  const char* reason{nullptr};

  WeaponType weapon{WeaponType::BeamLaser};
  int ammoMax{0};
  int haveAmmo{0};
  int needAmmo{0};

  int buyAmmo{0};

  double useWeapons{0.0};
  double useElectronics{0.0};

  double costCr{0.0};
  bool limitedByStock{false};
  bool limitedByCredits{false};
};

OrdnanceRearmQuote quoteOrdnanceRearmToFull(const econ::StationEconomyState& stEcon,
                                           const econ::StationEconomyModel& model,
                                           ShipHullClass hullClass,
                                           WeaponType weapon,
                                           core::u8 ammoCurrent,
                                           const StationServicePriceModel& pm,
                                           double creditsBudgetCr = -1.0);

struct OrdnanceRearmResult {
  bool ok{false};
  const char* reason{nullptr};

  int ammoBought{0};
  double creditsPaid{0.0};
  double weaponsTaken{0.0};
  double electronicsTaken{0.0};
};

OrdnanceRearmResult applyOrdnanceRearmToFull(econ::StationEconomyState& stEcon,
                                            const econ::StationEconomyModel& model,
                                            double& ioCredits,
                                            ShipHullClass hullClass,
                                            WeaponType weapon,
                                            core::u8& ioAmmo,
                                            const StationServicePriceModel& pm);

struct OrdnanceRearmAllQuote {
  bool ok{false};
  const char* reason{nullptr};

  WeaponType weaponPrimary{WeaponType::BeamLaser};
  WeaponType weaponSecondary{WeaponType::BeamLaser};

  int buyPrimary{0};
  int buySecondary{0};

  double useWeapons{0.0};
  double useElectronics{0.0};

  double costCr{0.0};
  bool limitedByStock{false};
  bool limitedByCredits{false};
};

OrdnanceRearmAllQuote quoteOrdnanceRearmAllToFull(const econ::StationEconomyState& stEcon,
                                                 const econ::StationEconomyModel& model,
                                                 ShipHullClass hullClass,
                                                 WeaponType weaponPrimary,
                                                 core::u8 ammoPrimary,
                                                 WeaponType weaponSecondary,
                                                 core::u8 ammoSecondary,
                                                 const StationServicePriceModel& pm,
                                                 double creditsBudgetCr = -1.0);

struct OrdnanceRearmAllResult {
  bool ok{false};
  const char* reason{nullptr};

  int primaryBought{0};
  int secondaryBought{0};

  double creditsPaid{0.0};
  double weaponsTaken{0.0};
  double electronicsTaken{0.0};
};

OrdnanceRearmAllResult applyOrdnanceRearmAllToFull(econ::StationEconomyState& stEcon,
                                                  const econ::StationEconomyModel& model,
                                                  double& ioCredits,
                                                  ShipHullClass hullClass,
                                                  WeaponType weaponPrimary,
                                                  core::u8& ioAmmoPrimary,
                                                  WeaponType weaponSecondary,
                                                  core::u8& ioAmmoSecondary,
                                                  const StationServicePriceModel& pm);

} // namespace stellar::sim
