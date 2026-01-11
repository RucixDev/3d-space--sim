#pragma once

#include "stellar/core/Types.h"

#include <algorithm>
#include <cstddef>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// ShipLoadout — shared gameplay balance definitions
// -----------------------------------------------------------------------------
//
// The SDL/OpenGL prototype previously kept hull/module/weapon tuning tables as
// local data inside apps/stellar_game/main.cpp.
//
// This header centralizes those *gameplay* definitions so:
//  - the meaning of SaveGame::shipHull / thrusterMk / weaponPrimary stays stable
//  - derived stat formulas can be unit tested
//  - future tools (e.g. stellar_sandbox) can reason about loadouts without
//    depending on the renderer app.
//
// IMPORTANT: numeric values here implicitly define save-file semantics.
// Keep enum values stable and do not reorder existing items.

enum class ShipHullClass : core::u8 {
  Scout  = 0,
  Hauler = 1,
  Fighter = 2,
};

enum class WeaponType : core::u8 {
  BeamLaser   = 0,
  PulseLaser  = 1,
  Cannon      = 2,
  Railgun     = 3,
  MiningLaser = 4,
  HomingMissile = 5,
  RadarMissile = 6,
};

struct HullDef {
  const char* name;
  double priceCr;
  double hullMax;
  double shieldBase;
  double linAccelKmS2;
  double angAccelRadS2;
  double cargoMult;
  double fuelMult;
};

struct MkDef {
  const char* name;
  double priceCr;
  double mult;
};

struct WeaponDef {
  WeaponType type;
  const char* name;
  double priceCr;
  double cooldownSimSec;
  double heatPerShot;
  double dmg;
  double rangeKm;
  double projSpeedKmS; // ignored for beam
  bool beam;
  bool guided;         // homing projectile (missile) if true
  double blastRadiusKm; // explosion/splash radius for guided munitions (0 for none)
  double turnRateRadS;  // max turn rate for guided munitions (0 for none)
  float r, g, b;       // HUD/beam tint (game-side convenience)
};

// ---- Static tuning tables ----
// Kept as inline variables so this file stays header-only.

inline constexpr HullDef kHullDefs[] = {
  // name,    price,   hull,  shield, linAccel, angAccel, cargoMult, fuelMult
  {"Scout",   0.0,     100.0, 100.0,  0.080,    1.20,     1.00,      1.00},
  {"Hauler",  14000.0, 125.0,  90.0,  0.070,    1.05,     1.35,      1.25},
  {"Fighter", 22000.0,  95.0, 120.0,  0.095,    1.35,     0.90,      1.05},
};

inline constexpr MkDef kThrusters[] = {
  {"(invalid)",        0.0,    1.0},
  {"Thrusters Mk1",    0.0,    1.00},
  {"Thrusters Mk2",    8500.0, 1.18},
  {"Thrusters Mk3",    17500.0, 1.35},
};

inline constexpr MkDef kShields[] = {
  {"(invalid)",     0.0,    1.0},
  {"Shields Mk1",   0.0,    1.00},
  {"Shields Mk2",   9000.0, 1.25},
  {"Shields Mk3",   18500.0, 1.55},
};

inline constexpr MkDef kDistributors[] = {
  {"(invalid)",         0.0,    1.0},
  {"Distributor Mk1",   0.0,    1.00},
  {"Distributor Mk2",   7500.0, 1.22},
  {"Distributor Mk3",   16000.0, 1.45},
};

inline constexpr WeaponDef kWeaponDefs[] = {
  // type,               name,          price,    cd,    heat, dmg,  range,    proj,   beam, guided, blast,  turn,   rgb
  {WeaponType::BeamLaser,   "Beam Laser",   0.0,     0.18,  2.2,  7.5, 210000.0,   0.0,  true,  false,  0.0,   0.0,  1.00f, 0.35f, 0.15f},
  {WeaponType::PulseLaser,  "Pulse Laser",  5200.0,  0.28,  1.9,  6.0, 230000.0,   0.0,  true,  false,  0.0,   0.0,  1.00f, 0.80f, 0.20f},
  {WeaponType::Cannon,      "Cannon",       0.0,     0.90,  4.5, 22.0, 260000.0, 120.0, false, false,  0.0,   0.0,  1.00f, 1.00f, 0.90f},
  {WeaponType::Railgun,     "Railgun",      9800.0,  1.65,  7.5, 45.0, 320000.0, 240.0, false, false,  0.0,   0.0,  0.60f, 0.90f, 1.00f},
  {WeaponType::MiningLaser, "Mining Laser", 6500.0,  0.22,  1.6,  3.0, 180000.0,   0.0,  true,  false,  0.0,   0.0,  0.30f, 1.00f, 0.35f},
  {WeaponType::HomingMissile, "Homing Missile", 12000.0, 2.80, 6.2, 60.0, 300000.0, 150.0, false, true,  5200.0, 0.55,  1.00f, 0.35f, 0.10f},
  {WeaponType::RadarMissile,  "Radar Missile",  14500.0, 3.20, 6.8, 58.0, 340000.0, 175.0, false, true,  5200.0, 0.70,  0.25f, 0.65f, 1.00f},
};

// ---- Derived stat helpers ----

struct ShipDerivedStats {
  double hullMax{100.0};
  double shieldMax{100.0};
  double shieldRegenPerSimMin{2.5};
  double heatCoolRate{10.0};

  // Baseline kinematic caps used by sim::Ship (km/s^2 and rad/s^2).
  double baseLinAccelKmS2{0.08};
  double baseAngAccelRadS2{1.2};
};

inline constexpr std::size_t hullDefCount() {
  return sizeof(kHullDefs) / sizeof(kHullDefs[0]);
}

inline constexpr std::size_t weaponDefCount() {
  return sizeof(kWeaponDefs) / sizeof(kWeaponDefs[0]);
}

inline const HullDef& hullDef(ShipHullClass cls) {
  const int idx = std::clamp((int)cls, 0, (int)hullDefCount() - 1);
  return kHullDefs[idx];
}

inline const WeaponDef& weaponDef(WeaponType t) {
  const int idx = std::clamp((int)t, 0, (int)weaponDefCount() - 1);
  return kWeaponDefs[idx];
}


inline bool weaponUsesAmmo(WeaponType t) {
  return weaponDef(t).guided;
}

// For now, only guided weapons consume ammunition.
// Capacity scales mildly with hull size so larger ships can carry more ordnance.
inline int weaponAmmoMax(WeaponType t, ShipHullClass hullClass) {
  if (!weaponUsesAmmo(t)) return 0;

  switch (t) {
    case WeaponType::HomingMissile:
      switch (hullClass) {
        case ShipHullClass::Scout:   return 8;
        case ShipHullClass::Hauler:  return 10;
        case ShipHullClass::Fighter: return 12;
      }
      return 8;

    case WeaponType::RadarMissile:
      switch (hullClass) {
        case ShipHullClass::Scout:   return 6;
        case ShipHullClass::Hauler:  return 8;
        case ShipHullClass::Fighter: return 10;
      }
      return 6;

    default:
      return 0;
  }
}

inline double weaponAmmoUnitPriceCr(WeaponType t) {
  if (!weaponUsesAmmo(t)) return 0.0;
  // Simple economy: ammo is a small fraction of the launcher purchase price.
  return std::max(30.0, weaponDef(t).priceCr * 0.012);
}


// Mirrors the gameplay formula historically used by stellar_game.
inline ShipDerivedStats computeShipDerivedStats(
  ShipHullClass hullClass,
  int thrusterMk,
  int shieldMk,
  int distributorMk
) {
  const HullDef& hull = hullDef(hullClass);

  const int tMk = std::clamp(thrusterMk, 1, 3);
  const int sMk = std::clamp(shieldMk, 1, 3);
  const int dMk = std::clamp(distributorMk, 1, 3);

  const double thrMult = kThrusters[tMk].mult;
  const double shMult = kShields[sMk].mult;
  const double distMult = kDistributors[dMk].mult;

  ShipDerivedStats out{};
  out.hullMax = hull.hullMax;
  out.shieldMax = hull.shieldBase * shMult;

  out.baseLinAccelKmS2 = hull.linAccelKmS2 * thrMult;
  // Angular tuning: slightly less sensitive than linear to keep higher Mk
  // upgrades from feeling too "twitchy".
  out.baseAngAccelRadS2 = hull.angAccelRadS2 * (0.92 + 0.08 * thrMult);

  // Regen & cooling scale mostly with distributor. Shields help slightly.
  out.shieldRegenPerSimMin = 2.5 * (0.85 + 0.15 * (double)sMk) * (0.85 + 0.15 * (double)dMk);
  out.heatCoolRate = 10.0 * distMult;

  return out;
}

} // namespace stellar::sim
