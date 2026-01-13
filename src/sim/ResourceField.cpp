#include "stellar/sim/ResourceField.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
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

// Sample a point in the unit sphere (uniform volume) via rejection sampling.
static math::Vec3d randInUnitSphere(core::SplitMix64& rng) {
  for (int i = 0; i < 24; ++i) {
    const double x = rng.range(-1.0, 1.0);
    const double y = rng.range(-1.0, 1.0);
    const double z = rng.range(-1.0, 1.0);
    const double d2 = x * x + y * y + z * z;
    if (d2 > 1e-12 && d2 <= 1.0) return {x, y, z};
  }
  return {0.0, 0.0, 0.0};
}

// Generate a candidate asteroid position inside the site's layout.
static math::Vec3d sampleLayoutPos(const ResourceFieldSite& site, core::SplitMix64& rng) {
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
      const math::Vec3d u = randInUnitSphere(rng);
      return site.posKm + X * (u.x * rx) + Y * (u.y * ry) + Z * (u.z * rz);
    }
    case ResourceFieldLayout::Sheet: {
      // Uniform disc in plane (X/Z) with a small thickness along Y.
      const double R = std::max(1.0, site.majorRadiusKm);
      const double T = std::max(0.0, site.minorRadiusKm);
      const double u0 = rng.nextDouble();
      const double u1 = rng.nextDouble();
      const double r = std::sqrt(std::max(0.0, u0)) * R;
      const double th = u1 * kTwoPi;
      const double y = rng.range(-T, T);
      return site.posKm + X * (r * std::cos(th)) + Z * (r * std::sin(th)) + Y * y;
    }
    case ResourceFieldLayout::Torus: {
      // Torus around the origin: majorRadiusKm is the ring radius, minorRadiusKm is tube radius.
      const double R = std::max(1.0, site.majorRadiusKm);
      const double rTube = std::max(0.0, site.minorRadiusKm);

      // Arc coverage.
      double phi = rng.nextDouble() * kTwoPi;
      if (site.arcRad < kTwoPi - 1e-6) {
        const double half = 0.5 * std::max(0.0, site.arcRad);
        phi = site.arcCenterRad + rng.range(-half, half);
      }

      const math::Vec3d radial = X * std::cos(phi) + Z * std::sin(phi);

      // Cross-section circle in the plane spanned by (radial, Y).
      const double u0 = rng.nextDouble();
      const double u1 = rng.nextDouble();
      const double rr = std::sqrt(std::max(0.0, u0)) * rTube;
      const double th = u1 * kTwoPi;
      const math::Vec3d off = radial * (rr * std::cos(th)) + Y * (rr * std::sin(th));

      return site.posKm + radial * R + off;
    }
  }

  return site.posKm;
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

ResourceFieldPlan generateResourceFields(core::u64 universeSeed,
                                        SystemId systemId,
                                        const math::Vec3d& anchorPosKm,
                                        double anchorCommsRangeKm,
                                        int fieldCount) {
  ResourceFieldPlan plan{};
  if (fieldCount <= 0) return plan;

  const core::u64 sysKey = core::hashCombine(universeSeed, static_cast<core::u64>(systemId));

  // Plan RNG is only used for top-level field selection; individual fields/asteroids
  // use their ids as seeds so results stay stable even if the generation order changes.
  // "RESFIELD" as an ASCII 64-bit tag.
  core::SplitMix64 prng(core::hashCombine(sysKey, 0x5245534649454C44ull));

  plan.fields.reserve(static_cast<std::size_t>(fieldCount));

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
      // Pick a stable plane normal for sheet/torus (and an arbitrary orientation for clusters).
      const math::Vec3d n = randUnitVec(frng);
      buildOrthonormalBasis(n, site.basisX, site.basisY, site.basisZ);

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

      // --- Position (layout-aware, with simple blue-noise spacing) ---
      // We use a deterministic rejection sampler to avoid obvious clumps/overlaps.
      // Field sizes are small (<=42 rocks), so O(n^2) checks are fine.
      const double padKm = 1800.0;
      const double baseMinSep = 2.0 * a.radiusKm + padKm;

      math::Vec3d pos = site.posKm;
      bool placed = false;

      // Rebuild a small local list of already-placed asteroids for this field.
      // NOTE: This is done per asteroid so we don't need to keep an additional
      // container in the plan; count is tiny so the overhead is negligible.
      // We only compare against earlier asteroids in the same field.
      for (int attempt = 0; attempt < 48 && !placed; ++attempt) {
        pos = sampleLayoutPos(site, posRng);

        const double minSep = std::max(4500.0, baseMinSep * (0.96 - 0.006 * (double)attempt));
        placed = true;

        // Compare against earlier asteroids in this field that are already in the plan.
        // We can scan backwards a small amount; worst case is small.
        const std::size_t scanLimit = 96; // at most 42 per field, but keep a small cap
        std::size_t scanned = 0;
        for (std::size_t k = plan.asteroids.size(); k > 0 && scanned < scanLimit; --k, ++scanned) {
          const auto& prev = plan.asteroids[k - 1];
          if (prev.fieldId != site.id) continue;
          const double d2 = (prev.posKm - pos).lengthSq();
          const double req = (prev.radiusKm + a.radiusKm) * 1.85 + padKm;
          const double req2 = std::max(minSep, req);
          if (d2 < req2 * req2) { placed = false; break; }
        }
      }
      a.posKm = pos;

      // Yield mix (primary dominates).
      a.yield = propRng.chance(site.secondaryChance) ? site.secondary : site.primary;

      // Units scale with radius and richness; secondary goods are slightly scarcer.
      double baseUnits = propRng.range(90.0, 260.0) * (a.radiusKm / 5000.0) * site.richness;
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
