#include "stellar/sim/ShipScan.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

double shipScanQuality01(double scanStrength01, double jammerPower) {
  const double s = std::clamp(scanStrength01, 0.0, 1.0);

  // Base confidence is biased low so the outer edge of the radar feels "fuzzy".
  double q = 0.15 + 0.85 * s;

  // Jammers reduce scan reliability even when the target is visible.
  // We keep the effect bounded so scans remain usable (just noisy).
  const double jam01 = std::clamp(jammerPower / 1.0, 0.0, 1.0);

  // Jamming hurts more when you're already weakly tracking the target.
  const double jamMult = 0.65 + 0.35 * (1.0 - s);
  q *= (1.0 - 0.35 * jam01 * jamMult);

  return std::clamp(q, 0.0, 1.0);
}

static double tankiness01(ShipHullClass cls, int thrusterMk, int shieldMk, int distributorMk) {
  const ShipDerivedStats stats = computeShipDerivedStats(cls, thrusterMk, shieldMk, distributorMk);
  const double tank = stats.hullMax + stats.shieldMax;

  // Calibrated against the current tuning tables in ShipLoadout.
  // Keep the mapping generous so future hulls/modules still fit.
  return std::clamp((tank - 170.0) / (360.0 - 170.0), 0.0, 1.0);
}

static double weaponThreat01(WeaponType t) {
  const WeaponDef& w = weaponDef(t);
  const double cd = std::max(0.05, w.cooldownSimSec);
  const double dps = w.dmg / cd;

  // Normalize by a high-DPS laser baseline.
  double out = std::clamp(dps / 45.0, 0.0, 1.0);

  // Guided weapons are bursty and harder to evade in the prototype.
  if (w.guided) out = std::clamp(out + 0.18, 0.0, 1.0);
  if (w.blastRadiusKm > 1.0) out = std::clamp(out + 0.08, 0.0, 1.0);
  if (w.rangeKm > 280000.0) out = std::clamp(out + 0.04, 0.0, 1.0);

  return out;
}

static double mkTier01(int thrusterMk, int shieldMk, int distributorMk) {
  const double avg = ((double)std::clamp(thrusterMk, 1, 3) + (double)std::clamp(shieldMk, 1, 3) + (double)std::clamp(distributorMk, 1, 3)) / 3.0;
  return std::clamp((avg - 1.0) / 2.0, 0.0, 1.0);
}

double shipThreatRating01(const ShipScanInput& in) {
  const double hullFrac = std::clamp(in.hullFrac, 0.0, 1.0);
  const double shieldFrac = std::clamp(in.shieldFrac, 0.0, 1.0);
  const double ai = std::clamp(in.aiSkill01, 0.0, 1.0);

  const double tank01 = tankiness01(in.hullClass, in.thrusterMk, in.shieldMk, in.distributorMk);
  const double wep01  = weaponThreat01(in.weapon);
  const double mk01   = mkTier01(in.thrusterMk, in.shieldMk, in.distributorMk);

  // Weighted blend: weapons dominate perceived threat, but a tanky ship matters.
  double threat = 0.40 * wep01 + 0.35 * tank01 + 0.15 * mk01 + 0.10 * ai;

  // Current condition influences how scary the ship is *right now*.
  threat *= (0.55 + 0.45 * hullFrac);
  threat *= (0.70 + 0.30 * shieldFrac);

  return std::clamp(threat, 0.0, 1.0);
}

ShipScanReport computeShipScanReport(core::u64 universeSeed, const ShipScanInput& in) {
  ShipScanReport out{};
  out.targetId = in.targetId;
  out.quality01 = shipScanQuality01(in.scanStrength01, in.jammerPower);
  out.threat01 = shipThreatRating01(in);

  const double q = out.quality01;

  // Identification thresholds.
  out.hullKnown = (q >= 0.25);
  out.hullClass = in.hullClass;

  out.weaponKnown = (q >= 0.35);
  out.weapon = in.weapon;

  out.healthKnown = (q >= 0.55);
  out.hullFrac = std::clamp(in.hullFrac, 0.0, 1.0);
  out.shieldFrac = std::clamp(in.shieldFrac, 0.0, 1.0);

  // Cargo baseline: if we know units, estimate value from base price.
  const double baseUnits = std::max(0.0, in.cargoUnits);
  const double baseValue = (baseUnits > 0.0)
                             ? (baseUnits * econ::commodityDef(in.cargoCommodity).basePrice)
                             : std::max(0.0, in.cargoValueCr);

  out.cargoDetected = (baseValue > 1.0);
  out.cargoKnown = out.cargoDetected && (q >= 0.45);
  out.cargoCommodityKnown = out.cargoDetected && (q >= 0.55);
  out.cargoUnitsKnown = out.cargoDetected && (q >= 0.65) && (baseUnits > 0.0);

  if (out.cargoCommodityKnown) out.cargoCommodity = in.cargoCommodity;
  out.cargoConfidence01 = out.cargoKnown ? std::clamp((q - 0.45) / 0.55, 0.0, 1.0) : 0.0;

  // Deterministic bounded noise.
  const double jam01 = std::clamp(in.jammerPower / 1.0, 0.0, 1.0);
  double err = 0.60 * (1.0 - q);
  err += 0.22 * jam01 * (1.0 - 0.55 * q);
  err = std::clamp(err, 0.05, 0.90);

  if (out.cargoDetected) {
    // Seed the RNG using stable inputs so reports remain deterministic for replays/tests.
    const core::u64 s0 = core::hashCombine(universeSeed, in.targetId);
    const core::u64 s1 = core::hashCombine(s0, core::fnv1a64("ship_scan"));
    core::SplitMix64 rng(s1);

    const double signedUnit = rng.nextDouble() * 2.0 - 1.0;
    const double mag = 0.35 + 0.65 * rng.nextDouble();
    const double f = 1.0 + signedUnit * err * mag;

    const double minV = baseValue * (1.0 - err);
    const double maxV = baseValue * (1.0 + err);

    out.cargoValueMinCr = std::max(0.0, minV);
    out.cargoValueMaxCr = std::max(out.cargoValueMinCr, maxV);

    double est = baseValue * f;
    est = std::clamp(est, out.cargoValueMinCr, out.cargoValueMaxCr);
    out.cargoValueEstCr = std::max(0.0, est);

    if (out.cargoUnitsKnown) {
      // Units are slightly noisier than value.
      const double uErr = std::clamp(err * 0.85, 0.05, 0.85);
      const double uSigned = rng.nextDouble() * 2.0 - 1.0;
      const double uMag = 0.35 + 0.65 * rng.nextDouble();
      const double uF = 1.0 + uSigned * uErr * uMag;
      out.cargoUnitsEst = std::max(0.0, baseUnits * uF);
    }
  }

  // EW: jammer detection.
  if (in.jammerPower > 0.05) {
    out.jammerSuspected = (q >= 0.20);
    out.jammerDetected = (q >= 0.40);

    if (out.jammerDetected) {
      // Estimated jammer strength is coarse (you rarely get perfect numbers).
      const double strength = std::clamp(in.jammerPower / 1.0, 0.0, 1.0);

      const core::u64 s0 = core::hashCombine(universeSeed, in.targetId);
      const core::u64 s1 = core::hashCombine(s0, core::fnv1a64("ship_scan_jam"));
      core::SplitMix64 rng(s1);

      // +/- 10% estimation error at perfect scans; larger at worse quality.
      const double e = std::clamp(0.10 + 0.35 * (1.0 - q), 0.10, 0.45);
      const double sign = rng.nextDouble() * 2.0 - 1.0;
      const double mag = 0.35 + 0.65 * rng.nextDouble();
      out.jammerStrength01 = std::clamp(strength * (1.0 + sign * e * mag), 0.0, 1.0);
    }
  }

  return out;
}

} // namespace stellar::sim
