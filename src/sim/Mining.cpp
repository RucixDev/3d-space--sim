#include "stellar/sim/Mining.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {
namespace {

static double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

// Convert a 64-bit hash to a stable double in [0, 1).
static double unitFromHash(core::u64 h) {
  // Use the top 53 bits as a mantissa to avoid precision issues.
  const core::u64 mant = (h >> 11) & ((core::u64(1) << 53) - 1);
  return (double)mant / (double)(core::u64(1) << 53);
}

static double volatileChance(ResourceFieldKind k) {
  // Small, field-dependent chance. These are gameplay knobs.
  switch (k) {
    case ResourceFieldKind::OreBelt: return 0.04;
    case ResourceFieldKind::MetalPocket: return 0.07;
    case ResourceFieldKind::IceField: return 0.02;
  }
  return 0.04;
}

} // namespace

MiningAsteroidTraits miningAsteroidTraits(core::u64 universeSeed,
                                         core::u64 asteroidId,
                                         ResourceFieldKind fieldKind) {
  MiningAsteroidTraits out{};

  // Deterministic hash key: seed + asteroid id + field kind.
  core::u64 h = core::hashCombine(universeSeed, asteroidId);
  h = core::hashCombine(h, (core::u64)(core::u8)fieldKind);

  const double u0 = unitFromHash(core::hashCombine(h, core::fnv1a64("volatile")));
  out.volatilePocket = (u0 < volatileChance(fieldKind));

  if (out.volatilePocket) {
    // Fracture threshold: fraction of remaining units where the seam breaks.
    // Keep it in a readable band so it doesn't trigger too early/late.
    const double u1 = unitFromHash(core::hashCombine(h, core::fnv1a64("fracture")));
    out.fractureFrac = 0.22 + 0.26 * u1; // ~[0.22, 0.48]
  }

  return out;
}

double miningEfficiency(double distKm,
                        double rangeKm,
                        double fullEfficiencyFrac,
                        double minEfficiency) {
  if (rangeKm <= 1e-9) return 0.0;

  const double ff = std::clamp(fullEfficiencyFrac, 0.0, 1.0);
  const double minEff = std::clamp(minEfficiency, 0.0, 1.0);

  const double frac = distKm / rangeKm;
  if (frac <= ff) return 1.0;
  if (frac >= 1.0) return minEff;

  // Linear drop from 1.0 at ff -> minEff at 1.0.
  const double t = (frac - ff) / std::max(1e-9, 1.0 - ff);
  const double eff = 1.0 - t * (1.0 - minEff);
  return std::clamp(eff, minEff, 1.0);
}

MiningHitResult computeMiningHit(const MiningHitInput& in) {
  MiningHitResult out{};

  const MiningAsteroidTraits traits = miningAsteroidTraits(in.universeSeed, in.asteroidId, in.fieldKind);
  out.volatilePocket = traits.volatilePocket;
  out.fractureFrac = traits.fractureFrac;

  if (in.remainingUnits <= 1e-9 || in.baseUnits <= 1e-9) {
    out.extractedUnits = 0.0;
    out.efficiency = 0.0;
    out.fractureTriggered = false;
    return out;
  }

  const double eff = miningEfficiency(std::max(0.0, in.distKm),
                                      std::max(1e-6, in.rangeKm));
  out.efficiency = eff;

  double units = std::max(0.0, in.baseUnitsPerHit) * eff;
  if (in.prospected) units *= 1.20;

  units = std::min(units, in.remainingUnits);
  out.extractedUnits = units;

  // Fracture event: triggers once when remaining crosses a deterministic threshold.
  out.fractureTriggered = false;
  if (traits.volatilePocket && !in.fractureAlreadyTriggered && traits.fractureFrac > 0.0) {
    const double thresh = std::max(0.0, in.baseUnits) * traits.fractureFrac;
    if (in.remainingUnits > thresh + 1e-9 && (in.remainingUnits - units) <= thresh + 1e-9) {
      out.fractureTriggered = true;
    }
  }

  return out;
}

} // namespace stellar::sim
