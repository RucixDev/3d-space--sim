#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/ShipLoadout.h"

namespace stellar::sim {

// -----------------------------------------------------------------------------
// PowerDistributor
// -----------------------------------------------------------------------------
// A lightweight 3-channel capacitor system inspired by classic space sims.
//
// Channels:
//  - ENG: boost availability
//  - WEP: weapon firing budget
//  - SYS: shield regeneration budget
//
// A "pips" allocation controls how recharge is distributed between channels.

static constexpr int kPipMax   = 4;
static constexpr int kPipTotal = 6;

struct Pips {
  int eng{2};
  int wep{2};
  int sys{2};
};

// Clamp pips to [0..kPipMax] and normalize the sum to kPipTotal.
// If the sum is 0, resets to the default 2/2/2.
void normalizePips(Pips& p);

struct DistributorConfig {
  // Capacitor capacities (arbitrary "energy units").
  double capEng{90.0};
  double capWep{90.0};
  double capSys{100.0};

  // Total recharge rate (energy units per simulated second).
  // Pips determine how this recharge is distributed.
  double rechargePerSimSec{24.0};

  // ENG drain while boosting (energy units per simulated second).
  double boostCostPerSimSec{15.0};

  // SYS drain per shield point regenerated.
  double shieldRegenCostPerPoint{0.6};
};

struct DistributorState {
  double eng{0.0};
  double wep{0.0};
  double sys{0.0};
};

// Returns a gameplay tuning profile for a given hull + distributor Mk.
DistributorConfig distributorConfig(ShipHullClass hull, int distributorMk);

// Convenience: start with full capacitors.
DistributorState makeFull(const DistributorConfig& cfg);

// Applies boost drain for dtSim and returns the fraction of dtSim that can be
// boosted in [0..1]. If ENG runs out part-way, it will be fully consumed.
//
// This is intended to support "partial boost" integration in callers:
//  - run ship physics for dtSim * frac with boost=true
//  - run remaining dtSim * (1-frac) with boost=false
double consumeBoostFraction(DistributorState& st, const DistributorConfig& cfg, double dtSim);

// Recharge capacitors for dtSim based on pip allocation.
// If a channel is full, its share is redistributed to non-full channels.
void stepDistributor(DistributorState& st, const DistributorConfig& cfg, const Pips& pips, double dtSim);

// SYS pips effect helper: multiplier applied to base shield regen.
// Designed so 2 pips ~= 1.0, 4 pips ~= 1.4, 0 pips ~= 0.6.
double shieldRegenMultiplierFromPips(int sysPips);

// Simple heuristic for weapon capacitor cost (energy units per shot).
// Cost scales with cooldown and per-shot damage.
double weaponCapacitorCost(const WeaponDef& w);

} // namespace stellar::sim
