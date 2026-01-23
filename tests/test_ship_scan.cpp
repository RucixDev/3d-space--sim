#include "test_harness.h"

#include "stellar/sim/ShipScan.h"

#include <cmath>

int test_ship_scan() {
  using namespace stellar;
  using namespace stellar::sim;

  int failures = 0;

  // Quality should rise with sensor strength and fall with jammer power.
  const double qLow = shipScanQuality01(0.20, 0.0);
  const double qHigh = shipScanQuality01(0.90, 0.0);
  CHECK(qHigh > qLow);

  const double qJammed = shipScanQuality01(0.90, 1.0);
  CHECK(qJammed < qHigh);

  // Threat rating sanity: a well-equipped fighter should score above a mining hauler.
  ShipScanInput fighter{};
  fighter.targetId = 1;
  fighter.hullClass = ShipHullClass::Fighter;
  fighter.thrusterMk = 2;
  fighter.shieldMk = 3;
  fighter.distributorMk = 2;
  fighter.weapon = WeaponType::Railgun;
  fighter.hullFrac = 1.0;
  fighter.shieldFrac = 1.0;
  fighter.aiSkill01 = 0.85;

  ShipScanInput hauler{};
  hauler.targetId = 2;
  hauler.hullClass = ShipHullClass::Hauler;
  hauler.thrusterMk = 1;
  hauler.shieldMk = 1;
  hauler.distributorMk = 1;
  hauler.weapon = WeaponType::MiningLaser;
  hauler.hullFrac = 1.0;
  hauler.shieldFrac = 1.0;
  hauler.aiSkill01 = 0.35;

  const double tF = shipThreatRating01(fighter);
  const double tH = shipThreatRating01(hauler);
  CHECK(tF > tH);

  // Cargo estimate error band should shrink with quality.
  ShipScanInput lowScan = fighter;
  lowScan.targetId = 42;
  lowScan.cargoCommodity = econ::CommodityId::Metals;
  lowScan.cargoUnits = 20.0;
  lowScan.scanStrength01 = 0.18;

  ShipScanInput highScan = lowScan;
  highScan.scanStrength01 = 0.95;

  const auto repLow = computeShipScanReport(999u, lowScan);
  const auto repHigh = computeShipScanReport(999u, highScan);

  CHECK(repLow.cargoDetected);
  CHECK(repHigh.cargoDetected);
  CHECK(repHigh.quality01 > repLow.quality01);

  const double bandLow = repLow.cargoValueMaxCr - repLow.cargoValueMinCr;
  const double bandHigh = repHigh.cargoValueMaxCr - repHigh.cargoValueMinCr;
  CHECK(bandHigh < bandLow);

  CHECK(repHigh.cargoKnown);
  CHECK(repHigh.cargoCommodityKnown);
  CHECK(repHigh.cargoUnitsKnown);

  // Jammer detection should be gated by quality.
  ShipScanInput jamLow = lowScan;
  jamLow.targetId = 77;
  jamLow.jammerPower = 1.0;
  jamLow.scanStrength01 = 0.15;

  ShipScanInput jamHigh = jamLow;
  jamHigh.scanStrength01 = 0.90;

  const auto repJamLow = computeShipScanReport(123u, jamLow);
  const auto repJamHigh = computeShipScanReport(123u, jamHigh);

  CHECK(!repJamLow.jammerDetected);
  CHECK(repJamHigh.jammerDetected);
  CHECK(repJamHigh.jammerStrength01 >= 0.0 && repJamHigh.jammerStrength01 <= 1.0);

  return failures;
}
