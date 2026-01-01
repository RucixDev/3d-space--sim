#include "stellar/sim/ResourceField.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

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

    const core::u64 typeCode = 0xF13D0000ull | static_cast<core::u64>(kind);
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

    plan.fields.push_back(site);

    // --- Asteroids ---
    const int baseCount = frng.range(22, 34);
    const int count = std::clamp((int)std::round(baseCount * site.richness), 18, 42);

    plan.asteroids.reserve(plan.asteroids.size() + static_cast<std::size_t>(count));

    for (int j = 0; j < count; ++j) {
      const core::u64 asteroidId = makeDeterministicWorldId(site.id, static_cast<core::u64>(j));

      core::SplitMix64 arng(core::hashCombine(site.id, static_cast<core::u64>(j)));

      ResourceAsteroid a{};
      a.id = asteroidId;
      a.fieldId = site.id;

      a.posKm = site.posKm + randUnitVec(arng) * arng.range(15000.0, 78000.0);
      a.radiusKm = arng.range(2500.0, 8200.0);

      // Yield mix (primary dominates).
      a.yield = arng.chance(site.secondaryChance) ? site.secondary : site.primary;

      // Units scale with radius and richness; secondary goods are slightly scarcer.
      double baseUnits = arng.range(90.0, 260.0) * (a.radiusKm / 5000.0) * site.richness;
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
