#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/PowerDistributor.h"
#include "stellar/sim/SecurityModel.h"
#include "stellar/sim/ShipLoadout.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// TrafficConvoyEncounter
// -----------------------------------------------------------------------------
//
// When a Traffic Convoy signal is "resolved" in the SDL prototype, the game can
// optionally "physicalize" the convoy into local-space contacts:
//   - a trader ship carrying the shipment
//   - optional police escorts (in lawful systems)
//   - optional pirate raiders (depending on piracy/security/cargo value)
//
// Historically this was implemented ad-hoc in apps/stellar_game/main.cpp.
// This module provides a deterministic, reusable encounter planner that is:
//   - headless (usable by stellar_sandbox and unit tests)
//   - stable (sub-stream RNG to avoid perturbing decisions when visuals change)
//   - compact (returns a small list of spawn specs for the app to instantiate)

enum class TrafficEncounterNpcRole : core::u8 {
  Trader = 0,
  Police = 1,
  Pirate = 2,
};

struct TrafficEncounterNpcSpec {
  TrafficEncounterNpcRole role{TrafficEncounterNpcRole::Trader};
  bool leader{false};

  // Spawn offsets relative to the Traffic Convoy signal position.
  math::Vec3d posOffsetKm{0, 0, 0};
  math::Vec3d velOffsetKmS{0, 0, 0};

  // Recommended loadout knobs (consumed by the app's NPC factory).
  ShipHullClass hullClass{ShipHullClass::Scout};
  int thrusterMk{1};
  int shieldMk{1};
  int distributorMk{1};
  WeaponType weapon{WeaponType::BeamLaser};

  // Multipliers for derived stats.
  double hullMul{1.0};
  double shieldMul{1.0};
  double regenMul{1.0};
  double accelMul{1.0};
  double aiSkill{0.55};

  // Recommended pips allocation.
  Pips pips{2, 2, 2};
};

struct TrafficConvoyEncounterParams {
  // Cargo value scaling used to compute value01.
  double valueScaleCr{25000.0};

  // Spawn distance bands (km).
  double convoyOffsetMinKm{3000.0};
  double convoyOffsetMaxKm{12000.0};

  double escortOffsetMinKm{14000.0};
  double escortOffsetMaxKm{24000.0};
  double escortVelJitterMaxKmS{0.6};

  double pirateOffsetMinKm{52000.0};
  double pirateOffsetMaxKm{90000.0};
  double pirateVelJitterMaxKmS{1.2};

  // Ambush risk model:
  //   risk = clamp(base + piracyW*piracy + valueW*value - securityW*security, 0..max)
  double pirateRiskBase{0.12};
  double pirateRiskPiracyWeight{0.65};
  double pirateRiskValueWeight{0.35};
  double pirateRiskSecurityWeight{0.55};
  double pirateRiskMax{0.90};

  // Escort scaling (lawful systems only).
  double escortSecurity2{0.70};
  double escortSecurity1{0.55};
  double escortExtraValueThresh{0.85};
  double escortExtraChance{0.55};
  int escortMax{3};

  // Pirate pack sizing (if an ambush occurs).
  double pirateExtra1Risk{0.60};
  double pirateExtra1Chance{0.45};
  double pirateExtra2Risk{0.82};
  double pirateExtra2Chance{0.25};
  int pirateMin{2};
  int pirateMax{5};

  // Soft cap: if the local space already has >= cap pirates alive, suppress
  // additional convoy raiders.
  int pirateAliveCap{10};

  // Baseline police regen tuning used to compute regenMul for escort ships.
  // (This matches the prototype's in-game tuning constant.)
  double npcShieldRegenPerSec{0.035};
};

struct TrafficConvoyEncounterPlan {
  bool ok{false};
  core::u64 seed{0};

  // Informational scalars (useful for UI/debug overlays).
  double cargoValueCr{0.0};
  double value01{0.0};
  double risk01{0.0};

  int escortCount{0};
  int pirateCount{0};
  bool ambush{false};

  // Ordering convention:
  //  - npc[0] is always the trader (convoy leader)
  //  - escorts (if any) follow
  //  - raiders (if any) come last
  std::vector<TrafficEncounterNpcSpec> npcs{};
};

// Plan the local-space contacts for a convoy encounter.
//
// The caller passes cargoValueCr (estimated by the app from a market), plus the
// effective security profile for the system.
//
// alivePiratesNow is a *local-space* count (not a system-wide metric). It is
// used to avoid exploding contact counts when the player chain-resolves many
// convoys while already fighting pirates.
TrafficConvoyEncounterPlan planTrafficConvoyEncounter(core::u64 universeSeed,
                                                      core::u64 convoyId,
                                                      double cargoValueCr,
                                                      const SystemSecurityProfile& sec,
                                                      bool lawfulSystem,
                                                      int alivePiratesNow,
                                                      const TrafficConvoyEncounterParams& params = {});

} // namespace stellar::sim
