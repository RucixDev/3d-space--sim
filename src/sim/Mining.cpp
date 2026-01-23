#include "stellar/sim/Mining.h"

#include "stellar/core/Hash.h"
#include "stellar/math/Math.h"

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

math::Vec3d miningSeamDir(core::u64 universeSeed,
                          core::u64 asteroidId,
                          ResourceFieldKind fieldKind) {
  // Deterministic hash key: seed + asteroid id + field kind.
  core::u64 h = core::hashCombine(universeSeed, asteroidId);
  h = core::hashCombine(h, (core::u64)(core::u8)fieldKind);
  h = core::hashCombine(h, core::fnv1a64("seamDir"));

  // Fixed-cost spherical sampling (no rejection sampling) for stable RNG usage.
  const double u0 = unitFromHash(core::hashCombine(h, core::fnv1a64("u0")));
  const double u1 = unitFromHash(core::hashCombine(h, core::fnv1a64("u1")));

  const double z = u0 * 2.0 - 1.0;
  const double a = u1 * (2.0 * math::kPi);
  const double r = std::sqrt(std::max(0.0, 1.0 - z * z));

  math::Vec3d v{std::cos(a) * r, std::sin(a) * r, z};
  const double lsq = v.lengthSq();
  if (lsq < 1e-12) return {1.0, 0.0, 0.0};
  return v / std::sqrt(lsq);
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

  out.seamMultiplier = 1.0;
  if (in.prospected) {
    // Prospecting baseline bonus.
    units *= 1.20;

    // Optional seam-aware yield: when the hit direction is known (non-zero),
    // scale yield based on how well the beam aligns with the asteroid's hidden
    // "rich seam" axis. The curve is centered so the average multiplier is 1.0
    // over random hit directions.
    if (in.hitDirUnit.lengthSq() > 1e-12) {
      const math::Vec3d hitDir = in.hitDirUnit.normalized();
      const math::Vec3d seamDir = miningSeamDir(in.universeSeed, in.asteroidId, in.fieldKind);
      const double u = std::clamp(std::abs(math::dot(hitDir, seamDir)), 0.0, 1.0);

      // Use u^2 to emphasize the poles and keep a broad "meh" band.
      // For u ~ U[0,1], E[u^2] = 1/3. We subtract the mean so E[mult]=1.
      const double u2 = u * u;
      constexpr double kMean = 1.0 / 3.0;
      constexpr double kStrength = 0.45; // gameplay knob (range ~[0.85, 1.30])
      out.seamMultiplier = 1.0 + kStrength * (u2 - kMean);
      out.seamMultiplier = std::clamp(out.seamMultiplier, 0.70, 1.40);

      units *= out.seamMultiplier;
    }
  }

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
