#include "stellar/sim/ResourceField.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/core/LowDiscrepancy.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

const char* resourceFieldKindName(ResourceFieldKind k) {
  switch (k) {
    case ResourceFieldKind::OreBelt: return "Ore Belt";
    case ResourceFieldKind::MetalPocket: return "Metal Pocket";
    case ResourceFieldKind::IceField: return "Ice Field";
    default: return "Resource Field";
  }
}

const char* resourceFieldLayoutName(ResourceFieldLayout l) {
  switch (l) {
    case ResourceFieldLayout::Cluster: return "Cluster";
    case ResourceFieldLayout::Torus: return "Torus";
    case ResourceFieldLayout::Sheet: return "Sheet";
    default: return "Layout";
  }
}

const char* resourceFieldFeatureKindName(ResourceFieldFeatureKind k) {
  switch (k) {
    case ResourceFieldFeatureKind::Hotspot: return "Hotspot";
    case ResourceFieldFeatureKind::Gap: return "Gap";
    case ResourceFieldFeatureKind::Streak: return "Streak";
    case ResourceFieldFeatureKind::Spokes: return "Spokes";
    default: return "Feature";
  }
}

static math::Vec3d randUnitVec(core::SplitMix64& rng) {
  // Sample a random vector in [-1,1]^3 and normalize.
  // If we hit a near-zero vector, fall back to +X.
  const double x = rng.range(-1.0, 1.0);
  const double y = rng.range(-1.0, 1.0);
  const double z = rng.range(-1.0, 1.0);

  const double lsq = x * x + y * y + z * z;
  if (lsq < 1e-12) return {1.0, 0.0, 0.0};
  const double invLen = 1.0 / std::sqrt(lsq);
  return {x * invLen, y * invLen, z * invLen};
}

static void buildOrthonormalBasis(const math::Vec3d& nIn,
                                 math::Vec3d& outX,
                                 math::Vec3d& outY,
                                 math::Vec3d& outZ) {
  // outY is the requested normal (or +Y if degenerate).
  math::Vec3d n = nIn.normalized();
  if (n.lengthSq() < 1e-12) n = {0.0, 1.0, 0.0};

  // Pick a helper axis that is not parallel to n.
  const math::Vec3d a = (std::abs(n.y) < 0.90) ? math::Vec3d{0.0, 1.0, 0.0} : math::Vec3d{1.0, 0.0, 0.0};
  math::Vec3d x = math::cross(a, n).normalized();
  if (x.lengthSq() < 1e-12) x = math::cross(math::Vec3d{0.0, 0.0, 1.0}, n).normalized();
  math::Vec3d z = math::cross(n, x).normalized();

  outX = x;
  outY = n;
  outZ = z;
}

static ResourceFieldLayout layoutFor(ResourceFieldKind k) {
  switch (k) {
    case ResourceFieldKind::OreBelt: return ResourceFieldLayout::Torus;
    case ResourceFieldKind::MetalPocket: return ResourceFieldLayout::Cluster;
    case ResourceFieldKind::IceField: return ResourceFieldLayout::Sheet;
    default: return ResourceFieldLayout::Cluster;
  }
}

// Wrap an angle to (-pi, pi].
static double wrapAnglePi(double a) {
  const double twoPi = 2.0 * math::kPi;
  a = std::fmod(a + math::kPi, twoPi);
  if (a < 0.0) a += twoPi;
  return a - math::kPi;
}

static double gauss01(double x, double sigma) {
  const double s = std::max(1e-9, sigma);
  const double t = x / s;
  return std::exp(-0.5 * t * t);
}

static double softSat01(double x) {
  // Smoothly map x>=0 to [0,1).
  return 1.0 - std::exp(-std::max(0.0, x));
}

// Generate a candidate asteroid position inside the site's layout from three
// uniform random numbers in [0,1).
static math::Vec3d sampleLayoutPos01(const ResourceFieldSite& site, double u0, double u1, double u2) {
  constexpr double kTwoPi = 6.28318530717958647692;

  const math::Vec3d X = site.basisX;
  const math::Vec3d Y = site.basisY;
  const math::Vec3d Z = site.basisZ;

  switch (site.layout) {
    case ResourceFieldLayout::Cluster: {
      // Ellipsoid-like cluster: radius in X/Z is major, radius in Y is minor.
      const double rx = std::max(1.0, site.majorRadiusKm);
      const double ry = std::max(1.0, site.minorRadiusKm);
      const double rz = std::max(1.0, site.majorRadiusKm);

      // Uniform volume in unit sphere:
      //  - pick a direction uniformly on the unit sphere from (u1,u2)
      //  - pick a radius ~ cbrt(u0)
      const double r = std::cbrt(std::clamp(u0, 0.0, 1.0));
      const double z0 = 1.0 - 2.0 * std::clamp(u1, 0.0, 1.0);
      const double phi = std::clamp(u2, 0.0, 1.0) * kTwoPi;
      const double s = std::sqrt(std::max(0.0, 1.0 - z0 * z0));

      const math::Vec3d dir{ s * std::cos(phi), z0, s * std::sin(phi) };
      return site.posKm + X * (dir.x * r * rx) + Y * (dir.y * r * ry) + Z * (dir.z * r * rz);
    }
    case ResourceFieldLayout::Sheet: {
      // Uniform disc in plane (X/Z) with a small thickness along Y.
      const double R = std::max(1.0, site.majorRadiusKm);
      const double T = std::max(0.0, site.minorRadiusKm);

      const double r = std::sqrt(std::max(0.0, std::clamp(u0, 0.0, 1.0))) * R;
      const double th = std::clamp(u1, 0.0, 1.0) * kTwoPi;
      const double y = (std::clamp(u2, 0.0, 1.0) * 2.0 - 1.0) * T;

      return site.posKm + X * (r * std::cos(th)) + Z * (r * std::sin(th)) + Y * y;
    }
    case ResourceFieldLayout::Torus: {
      // Torus around the origin: majorRadiusKm is the ring radius, minorRadiusKm is tube radius.
      const double R = std::max(1.0, site.majorRadiusKm);
      const double rTube = std::max(0.0, site.minorRadiusKm);

      // Arc coverage.
      double phi = std::clamp(u0, 0.0, 1.0) * kTwoPi;
      if (site.arcRad < kTwoPi - 1e-6) {
        const double half = 0.5 * std::max(0.0, site.arcRad);
        phi = site.arcCenterRad + (std::clamp(u0, 0.0, 1.0) * 2.0 - 1.0) * half;
      }

      const math::Vec3d radial = X * std::cos(phi) + Z * std::sin(phi);

      // Cross-section circle in the plane spanned by (radial, Y).
      const double rr = std::sqrt(std::max(0.0, std::clamp(u1, 0.0, 1.0))) * rTube;
      const double th = std::clamp(u2, 0.0, 1.0) * kTwoPi;
      const math::Vec3d off = radial * (rr * std::cos(th)) + Y * (rr * std::sin(th));

      return site.posKm + radial * R + off;
    }
  }

  return site.posKm;
}

// Stochastic fallback sampler.
static math::Vec3d sampleLayoutPos(const ResourceFieldSite& site, core::SplitMix64& rng) {
  return sampleLayoutPos01(site, rng.nextDouble(), rng.nextDouble(), rng.nextDouble());
}


struct KindProfile {
  ResourceFieldKind kind{ResourceFieldKind::OreBelt};
  econ::CommodityId primary{econ::CommodityId::Ore};
  econ::CommodityId secondary{econ::CommodityId::Metals};
  double secondaryChance{0.0};
  double richnessBase{1.0};
};

static KindProfile profileFor(ResourceFieldKind k) {
  switch (k) {
    case ResourceFieldKind::OreBelt:
      return {k, econ::CommodityId::Ore, econ::CommodityId::Metals, 0.18, 1.00};
    case ResourceFieldKind::MetalPocket:
      return {k, econ::CommodityId::Metals, econ::CommodityId::Ore, 0.28, 1.18};
    case ResourceFieldKind::IceField:
      return {k, econ::CommodityId::Water, econ::CommodityId::Ore, 0.08, 0.92};
    default:
      return {k, econ::CommodityId::Ore, econ::CommodityId::Metals, 0.15, 1.00};
  }
}

static ResourceFieldKind pickKind(core::SplitMix64& rng) {
  // Weighted pick.
  const double r = rng.nextDouble();
  if (r < 0.60) return ResourceFieldKind::OreBelt;
  if (r < 0.85) return ResourceFieldKind::MetalPocket;
  return ResourceFieldKind::IceField;
}

std::vector<ResourceFieldFeature> filterFeaturesForField(const std::vector<ResourceFieldFeature>& features,
                                                         core::u64 fieldId) {
  std::vector<ResourceFieldFeature> out;
  out.reserve(features.size());
  for (const auto& f : features) {
    if (f.fieldId == fieldId) out.push_back(f);
  }
  return out;
}

static math::Vec3d fieldLocalKm(const ResourceFieldSite& site, const math::Vec3d& worldPosKm) {
  const math::Vec3d d = worldPosKm - site.posKm;
  return {math::dot(d, site.basisX), math::dot(d, site.basisY), math::dot(d, site.basisZ)};
}

static double torusArcMask01(const ResourceFieldSite& site, double phiRad) {
  constexpr double kTwoPi = 6.28318530717958647692;
  if (site.layout != ResourceFieldLayout::Torus) return 1.0;
  if (site.arcRad >= kTwoPi - 1e-6) return 1.0;

  const double half = 0.5 * std::max(0.0, site.arcRad);
  const double d = wrapAnglePi(phiRad - site.arcCenterRad);
  // Smooth falloff at arc boundaries.
  const double a = std::abs(d);
  const double edge = std::max(1e-6, half);
  if (a <= edge) {
    // Fade over ~10 degrees to avoid hard steps in heatmaps.
    const double fade = math::degToRad(10.0);
    const double t = math::clamp((edge - a) / std::max(1e-6, fade), 0.0, 1.0);
    return t;
  }
  return 0.0;
}

static double density01Local(const ResourceFieldSite& site,
                             const std::vector<ResourceFieldFeature>& features,
                             const math::Vec3d& localKm) {
  const double x = localKm.x;
  const double y = localKm.y;
  const double z = localKm.z;

  const double major = std::max(1.0, site.majorRadiusKm);
  const double minor = std::max(1.0, site.minorRadiusKm);

  if (site.layout == ResourceFieldLayout::Cluster) {
    // Normalize to unit-ish ellipsoid.
    const math::Vec3d p{x / major, y / minor, z / major};
    const double r2 = p.lengthSq();
    const double base = std::exp(-r2 / (2.0 * 0.62 * 0.62));

    double peaks = 0.0;
    for (const auto& f : features) {
      if (f.fieldId != site.id) continue;
      if (f.kind != ResourceFieldFeatureKind::Hotspot) continue;
      const double sigma = std::max(0.05, f.width);
      const math::Vec3d dp = p - f.localPos;
      const double d2 = dp.lengthSq();
      peaks += std::clamp(f.strength01, 0.0, 1.0) * std::exp(-d2 / (2.0 * sigma * sigma));
    }

    const double peak01 = softSat01(peaks);
    const double dens = (0.25 + 0.55 * base) + 0.55 * peak01;
    return std::clamp(dens, 0.0, 1.0);
  }

  if (site.layout == ResourceFieldLayout::Sheet) {
    const double r = std::sqrt(x * x + z * z);
    const double rr = std::clamp(r / major, 0.0, 1.0);
    // Dense toward the center; fade at the outer edge.
    double base = 1.0 - rr * rr;
    base = base * base;

    // Thin the sheet in Y.
    const double thick = gauss01(y, minor * 0.85 + 1.0);

    double peaks = 0.0;
    for (const auto& f : features) {
      if (f.fieldId != site.id) continue;
      if (f.kind != ResourceFieldFeatureKind::Streak) continue;
      const double theta = f.angleRad;
      const double w = std::max(150.0, f.width);

      // Perpendicular distance to the streak line (through the origin) in local X/Z.
      const double d = std::abs(x * std::sin(theta) - z * std::cos(theta));
      double s = gauss01(d, w);

      // Fade streaks outward a bit so they don't dominate the whole disc.
      s *= (0.35 + 0.65 * (1.0 - rr));

      peaks += std::clamp(f.strength01, 0.0, 1.0) * s;
    }

    const double peak01 = softSat01(peaks);
    const double dens = (0.35 + 0.65 * base) * thick * std::clamp(0.60 + 0.70 * peak01, 0.0, 1.0);
    return std::clamp(dens, 0.0, 1.0);
  }

  // Torus (belt)
  const double r = std::sqrt(x * x + z * z);
  const double dr = r - major;
  const double tube = std::sqrt(dr * dr + y * y);
  const double tube01 = std::clamp(1.0 - tube / (minor * 1.05 + 1.0), 0.0, 1.0);

  const double phi = std::atan2(z, x);
  const double arcMask = torusArcMask01(site, phi);
  if (arcMask <= 1e-6) return 0.0;

  double dipMul = 1.0;
  double peaks = 0.0;
  double spokeMul = 1.0;

  for (const auto& f : features) {
    if (f.fieldId != site.id) continue;
    if (f.kind == ResourceFieldFeatureKind::Gap) {
      const double w = std::max(1e-4, f.width);
      const double d = wrapAnglePi(phi - f.angleRad);
      const double g = gauss01(d, w);
      dipMul *= (1.0 - std::clamp(f.strength01, 0.0, 1.0) * g);
    } else if (f.kind == ResourceFieldFeatureKind::Hotspot) {
      const double w = std::max(1e-4, f.width);
      const double d = wrapAnglePi(phi - f.angleRad);
      const double g = gauss01(d, w);
      peaks += std::clamp(f.strength01, 0.0, 1.0) * g;
    } else if (f.kind == ResourceFieldFeatureKind::Spokes) {
      // param encodes frequency; angleRad is phase.
      const int m = std::clamp((int)std::llround(f.param), 3, 24);
      const double amp = std::clamp(f.strength01, 0.0, 1.0) * std::clamp(f.width, 0.0, 1.0);
      const double t = 0.5 + 0.5 * std::cos((double)m * (phi - f.angleRad));
      // Blend between flat (1.0) and spokes (0.65..1.0).
      const double spoke = 0.65 + 0.35 * t;
      spokeMul *= math::lerp(1.0, spoke, amp);
    }
  }

  const double peak01 = softSat01(peaks);
  const double base = 0.22 + 0.78 * std::pow(tube01, 0.75);
  const double ang = std::clamp(0.70 + 0.60 * peak01, 0.0, 1.0) * std::clamp(dipMul, 0.0, 1.0) * spokeMul;

  return std::clamp(base * ang * arcMask, 0.0, 1.0);
}

double resourceFieldDensity01(const ResourceFieldSite& site,
                              const std::vector<ResourceFieldFeature>& features,
                              const math::Vec3d& worldPosKm) {
  return density01Local(site, features, fieldLocalKm(site, worldPosKm));
}

ResourceFieldPlan generateResourceFields(core::u64 universeSeed,
                                        SystemId systemId,
                                        const math::Vec3d& anchorPosKm,
                                        double anchorCommsRangeKm,
                                        int fieldCount,
                                        math::Vec3d preferredPlaneNormalKm) {
  ResourceFieldPlan plan{};
  if (fieldCount <= 0) return plan;

  const core::u64 sysKey = core::hashCombine(universeSeed, static_cast<core::u64>(systemId));

  // Plan RNG is only used for top-level field selection; individual fields/asteroids
  // use their ids as seeds so results stay stable even if the generation order changes.
  // "RESFIELD" as an ASCII 64-bit tag.
  core::SplitMix64 prng(core::hashCombine(sysKey, 0x5245534649454C44ull));

  plan.fields.reserve(static_cast<std::size_t>(fieldCount));
  // Features are small and sparse; reserve a bit to avoid churn.
  plan.features.reserve(static_cast<std::size_t>(fieldCount) * 12u);
  plan.features.reserve(static_cast<std::size_t>(fieldCount) * 12);

  for (int i = 0; i < fieldCount; ++i) {
    const ResourceFieldKind kind = pickKind(prng);
    const auto prof = profileFor(kind);

    // Keep id scheme aligned with the game prototype (typeCode=1 for resource fields)
    // so asteroid depletion state can be shared when the renderer integrates this module.
    const core::u64 typeCode = 1ull;
    const core::u64 fieldId = makeDeterministicWorldId(core::hashCombine(sysKey, typeCode), static_cast<core::u64>(i));

    core::SplitMix64 frng(core::hashCombine(fieldId, 0xA11CE5EEDull));
    const math::Vec3d dir = randUnitVec(frng);

    // Keep fields close enough to be discoverable from the station, but distinct.
    const double distKm = anchorCommsRangeKm * (1.25 + 0.12 * (double)i) + 110000.0 + frng.range(0.0, 45000.0);

    ResourceFieldSite site{};
    site.id = fieldId;
    site.kind = kind;
    site.posKm = anchorPosKm + dir * distKm;

    const double richnessJitter = frng.range(0.85, 1.35);
    site.richness = std::clamp(prof.richnessBase * richnessJitter, 0.70, 1.45);

    site.primary = prof.primary;
    site.secondary = prof.secondary;
    site.secondaryChance = std::clamp(prof.secondaryChance * frng.range(0.85, 1.15), 0.02, 0.45);

    // Layout geometry (stable per-field id).
    site.layout = layoutFor(kind);
    {
      // Pick a stable plane normal for sheet/torus.
      // If the caller provides a preferred plane normal (e.g., the anchor station's
      // orbital angular momentum direction), align belts/sheets to that plane.
      math::Vec3d n = preferredPlaneNormalKm;
      if (n.lengthSq() < 1e-12) {
        n = randUnitVec(frng);
      }
      buildOrthonormalBasis(n, site.basisX, site.basisY, site.basisZ);

      // Randomize the in-plane rotation (stable per-field) so multiple fields can share
      // the same plane without being axis-aligned to each other.
      constexpr double kTwoPi = 6.28318530717958647692;
      const double rot = frng.nextDouble() * kTwoPi;
      const double c = std::cos(rot);
      const double s = std::sin(rot);
      const math::Vec3d X0 = site.basisX;
      const math::Vec3d Z0 = site.basisZ;
      site.basisX = X0 * c + Z0 * s;
      site.basisZ = Z0 * c - X0 * s;

      // Radii tuned to sit in the same overall size band as the legacy generator
      // (so field discoverability / scanning doesn't change drastically).
      const double rich = std::clamp(site.richness, 0.70, 1.45);

      if (site.layout == ResourceFieldLayout::Torus) {
        // Belt radius and thickness.
        site.majorRadiusKm = frng.range(42000.0, 76000.0) * (0.88 + 0.26 * rich);
        site.minorRadiusKm = frng.range(7000.0, 19000.0) * (0.90 + 0.20 * rich);
        site.minorRadiusKm = std::min(site.minorRadiusKm, site.majorRadiusKm * 0.85);

        // Occasionally generate an arc segment instead of a full ring.
        constexpr double kTwoPi = 6.28318530717958647692;
        if (frng.chance(0.22)) {
          site.arcRad = frng.range(0.95 * math::kPi, 1.85 * math::kPi);
          site.arcCenterRad = frng.range(0.0, kTwoPi);
        } else {
          site.arcRad = kTwoPi;
          site.arcCenterRad = 0.0;
        }
      } else if (site.layout == ResourceFieldLayout::Sheet) {
        // Large disc with small thickness.
        site.majorRadiusKm = frng.range(52000.0, 90000.0) * (0.92 + 0.22 * rich);
        site.minorRadiusKm = frng.range(2600.0, 9800.0) * (0.85 + 0.20 * rich);
        site.arcRad = 0.0;
        site.arcCenterRad = 0.0;
      } else {
        // Cluster (ellipsoid-ish). Make pockets tighter than belts.
        site.majorRadiusKm = frng.range(36000.0, 72000.0) * (0.90 + 0.25 * rich);
        const double flat = frng.range(0.55, 0.95);
        site.minorRadiusKm = std::max(8000.0, site.majorRadiusKm * flat);
        site.arcRad = 0.0;
        site.arcCenterRad = 0.0;
      }
    }

    // --- Structural features (gaps, hotspots, streaks, spokes) ---
    //
    // These are used by the density function to make fields feel less uniform
    // while still being fully deterministic.
    {
      core::SplitMix64 srng(core::hashCombine(site.id, core::fnv1a64("struct")));
      constexpr double kTwoPi = 6.28318530717958647692;

      const auto add = [&](ResourceFieldFeatureKind kind) -> ResourceFieldFeature& {
        plan.features.push_back({});
        auto& f = plan.features.back();
        f.fieldId = site.id;
        f.kind = kind;
        return f;
      };

      if (site.layout == ResourceFieldLayout::Torus) {
        const double arcFrac = std::clamp(site.arcRad / kTwoPi, 0.0, 1.0);

        // Belt gaps (think resonant depletion bands / shepherding).
        int gapCount = srng.range(1, 3);
        if (arcFrac < 0.80) gapCount = std::min(gapCount, 1);
        for (int g = 0; g < gapCount; ++g) {
          auto& f = add(ResourceFieldFeatureKind::Gap);
          const double center = (site.arcRad < kTwoPi - 1e-6)
                                    ? (site.arcCenterRad + srng.range(-0.45, 0.45) * site.arcRad)
                                    : srng.range(0.0, kTwoPi);
          f.angleRad = center;
          f.width = srng.range(math::degToRad(8.0), math::degToRad(22.0));
          f.strength01 = srng.range(0.55, 0.92);
        }

        // Hotspots / clumps (collision families, rubble piles).
        int hotCount = srng.range(2, 5);
        if (site.richness > 1.15) hotCount += 1;
        hotCount = std::clamp(hotCount, 1, 6);
        for (int h = 0; h < hotCount; ++h) {
          auto& f = add(ResourceFieldFeatureKind::Hotspot);
          const double center = (site.arcRad < kTwoPi - 1e-6)
                                    ? (site.arcCenterRad + srng.range(-0.48, 0.48) * site.arcRad)
                                    : srng.range(0.0, kTwoPi);
          f.angleRad = center;
          f.width = srng.range(math::degToRad(6.0), math::degToRad(18.0));
          f.strength01 = srng.range(0.30, 0.95);
        }

        // Subtle azimuthal spokes (texture-like modulation).
        if (srng.chance(0.55)) {
          auto& f = add(ResourceFieldFeatureKind::Spokes);
          f.angleRad = srng.range(0.0, kTwoPi);  // phase
          f.param = (double)srng.range(6, 14);   // frequency
          f.width = srng.range(0.35, 0.95);      // amplitude control (clamped in density fn)
          f.strength01 = srng.range(0.25, 0.80);
        }
      } else if (site.layout == ResourceFieldLayout::Sheet) {
        // Filament-like streaks.
        int streakCount = srng.range(2, 5);
        if (site.kind == ResourceFieldKind::IceField) streakCount += 1;
        streakCount = std::clamp(streakCount, 1, 6);

        for (int s = 0; s < streakCount; ++s) {
          auto& f = add(ResourceFieldFeatureKind::Streak);
          f.angleRad = srng.range(0.0, math::kPi); // direction (pi-periodic)
          f.width = std::max(250.0, site.majorRadiusKm * srng.range(0.035, 0.095));
          f.strength01 = srng.range(0.25, 0.90);
        }
      } else {
        // Cluster pockets/sub-clumps.
        int hotCount = srng.range(2, 4);
        if (site.kind == ResourceFieldKind::MetalPocket) hotCount += 1;
        hotCount = std::clamp(hotCount, 1, 6);

        for (int h = 0; h < hotCount; ++h) {
          auto& f = add(ResourceFieldFeatureKind::Hotspot);
          // Sample a point inside a unit sphere and keep it away from the edge.
          const double u0 = srng.nextDouble();
          const double u1 = srng.nextDouble();
          const double u2 = srng.nextDouble();

          const double r = std::cbrt(std::clamp(u0, 0.0, 1.0)) * 0.78;
          const double z0 = 1.0 - 2.0 * std::clamp(u1, 0.0, 1.0);
          const double phi = std::clamp(u2, 0.0, 1.0) * kTwoPi;
          const double s = std::sqrt(std::max(0.0, 1.0 - z0 * z0));

          f.localPos = {s * std::cos(phi) * r, z0 * r, s * std::sin(phi) * r};
          f.width = srng.range(0.14, 0.38);
          f.strength01 = srng.range(0.35, 0.95);
        }
      }
    }

    plan.fields.push_back(site);

    // --- Asteroids ---
    const int baseCount = frng.range(22, 34);
    int count = (int)std::round(baseCount * site.richness);
    count = std::clamp(count, 18, 42);

    // Back-compat: the original renderer prototype spawned 28 asteroids per persistent field.
    // Keep a minimum of 28 so older saves keyed by asteroid id don't "lose" mined rocks after
    // generator tuning changes.
    count = std::max(count, 28);

    plan.asteroids.reserve(plan.asteroids.size() + static_cast<std::size_t>(count));

    // Quasi-random (low discrepancy) candidate set rotation for this field.
    // This reduces structured aliasing from pure Halton points while keeping determinism.
    core::SplitMix64 qmcRng(core::hashCombine(site.id, core::fnv1a64("qmc")));
    const double qmcShift0 = qmcRng.nextDouble();
    const double qmcShift1 = qmcRng.nextDouble();
    const double qmcShift2 = qmcRng.nextDouble();
    constexpr int kCandidates = 16;

    for (int j = 0; j < count; ++j) {
      const core::u64 asteroidId = makeDeterministicWorldId(site.id, static_cast<core::u64>(j));

      // Use separate RNG streams so asteroid yield/units don't depend on how many
      // position re-sampling attempts were needed.
      core::SplitMix64 propRng(core::hashCombine(site.id, static_cast<core::u64>(j)));
      core::SplitMix64 posRng(core::hashCombine(asteroidId, core::fnv1a64("pos")));

      ResourceAsteroid a{};
      a.id = asteroidId;
      a.fieldId = site.id;

      // Size first (keeps radius stable even if we need to resample position for spacing).
      a.radiusKm = propRng.range(2500.0, 8200.0);

      // --- Position (layout-aware, deterministic blue-noise spacing) ---
      //
      // We place asteroids using a fixed-size quasi-random candidate set (Halton bases 2/3/5)
      // and pick the candidate that maximizes clearance to previously placed asteroids in the
      // same field (Mitchell "best-candidate" style). This tends to avoid visible clumps
      // without relying on unbounded rejection loops.
      const double padKm = 1800.0;
      const double baseMinSep = 2.0 * a.radiusKm + padKm;

      auto minClearanceKm = [&](const math::Vec3d& candidate, double dens01) -> double {
        double minClear = 1.0e30;
        for (const auto& prev : plan.asteroids) {
          if (prev.fieldId != site.id) continue;
          const double d = (prev.posKm - candidate).length();
          // Variable-density spacing: in hotspots we allow tighter packing, while
          // voids have a slightly larger exclusion radius.
          const double densPair = 0.5 * (std::clamp(dens01, 0.0, 1.0) + std::clamp(prev.density01, 0.0, 1.0));
          const double sepMul = math::lerp(1.28, 0.78, densPair);
          const double baseReq = std::max(baseMinSep, (prev.radiusKm + a.radiusKm) * 1.85 + padKm);
          const double req = baseReq * sepMul;
          minClear = std::min(minClear, d - req);
        }
        return minClear;
      };

      math::Vec3d pos = site.posKm;
      double bestScore = -1.0e30;
      double bestClear = -1.0e30;
      double bestDens = 1.0;

      // Bias: density is in [0,1] so this is a modest tie-breaker.
      const double densityBiasKm = baseMinSep * 0.45;

      // Candidate set: unique slice per asteroid index so patterns don't repeat.
      for (int c = 0; c < kCandidates; ++c) {
        const std::uint32_t idx = static_cast<std::uint32_t>(j * kCandidates + c + 1);
        const auto h = core::halton3(idx);

        const double u0 = core::frac(h.x + qmcShift0);
        const double u1 = core::frac(h.y + qmcShift1);
        const double u2 = core::frac(h.z + qmcShift2);

        const math::Vec3d candidate = sampleLayoutPos01(site, u0, u1, u2);
        const double dens01 = resourceFieldDensity01(site, plan.features, candidate);
        const double clearance = minClearanceKm(candidate, dens01);
        const double score = clearance + dens01 * densityBiasKm;

        if (score > bestScore) {
          bestScore = score;
          bestClear = clearance;
          bestDens = dens01;
          pos = candidate;
        }
      }

      bool placed = (bestClear >= 0.0) || (j == 0);

      // Fallback stochastic search (rare; mostly triggers for unusually dense pockets).
      for (int attempt = 0; attempt < 96 && !placed; ++attempt) {
        const math::Vec3d candidate = sampleLayoutPos(site, posRng);
        const double dens01 = resourceFieldDensity01(site, plan.features, candidate);
        const double clearance = minClearanceKm(candidate, dens01);
        if (clearance >= 0.0) {
          pos = candidate;
          bestDens = dens01;
          placed = true;
        }
      }

      a.posKm = pos;
      a.density01 = std::clamp(bestDens, 0.0, 1.0);

      // Yield mix (primary dominates). In hotspots, secondary deposits are a bit
      // more likely to appear.
      double secondaryChance = site.secondaryChance;
      secondaryChance *= (0.72 + 0.70 * a.density01);
      secondaryChance = std::clamp(secondaryChance, 0.01, 0.80);
      a.yield = propRng.chance(secondaryChance) ? site.secondary : site.primary;

      // Units scale with radius, richness, and local density. Secondary goods are slightly scarcer.
      double baseUnits = propRng.range(90.0, 260.0) * (a.radiusKm / 5000.0) * site.richness;
      baseUnits *= (0.80 + 0.60 * a.density01);
      if (a.yield == site.secondary) baseUnits *= 0.68;
      if (a.yield == econ::CommodityId::Water) baseUnits *= 0.80;
      a.baseUnits = std::max(1.0, baseUnits);

      plan.asteroids.push_back(a);
    }
  }

  return plan;
}

std::vector<ResourceAsteroid> filterAsteroidsForField(const std::vector<ResourceAsteroid>& asteroids,
                                                      core::u64 fieldId) {
  std::vector<ResourceAsteroid> out;
  out.reserve(asteroids.size());
  for (const auto& a : asteroids) {
    if (a.fieldId == fieldId) out.push_back(a);
  }
  return out;
}

} // namespace stellar::sim
