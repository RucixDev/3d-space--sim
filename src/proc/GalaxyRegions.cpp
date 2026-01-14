#include "stellar/proc/GalaxyRegions.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/proc/NameGenerator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace stellar::proc {
namespace {

static double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

struct FeaturePoint {
  core::u64 id{0};
  math::Vec3d p{0,0,0};
};

static FeaturePoint featurePointForCell(core::u64 universeSeed, core::i64 cx, core::i64 cy, core::i64 cz, double cellSizeLy) {
  // Combine seed + cell coord into a deterministic per-cell RNG.
  core::u64 h = core::hashCombine(universeSeed, core::fnv1a64("galaxy_region_cell_v1"));
  h = core::hashCombine(h, static_cast<core::u64>(cx));
  h = core::hashCombine(h, static_cast<core::u64>(cy));
  h = core::hashCombine(h, static_cast<core::u64>(cz));

  core::SplitMix64 rng(h);

  const math::Vec3d cellMin{(double)cx * cellSizeLy, (double)cy * cellSizeLy, (double)cz * cellSizeLy};
  const math::Vec3d p = cellMin + math::Vec3d{
    rng.range(0.0, cellSizeLy),
    rng.range(0.0, cellSizeLy),
    rng.range(0.0, cellSizeLy)
  };

  // Stable region id. We intentionally avoid using the random point coordinates.
  const core::u64 id = core::hashCombine(core::hashCombine(universeSeed, core::fnv1a64("galaxy_region_id_v1")), h);
  return FeaturePoint{id, p};
}

static GalaxyRegionKind classifyRegionKind(core::u64 regionSeed, const math::Vec3d& centerLy) {
  // A lightweight stylized classifier.
  // We bias by galactocentric radius so the galaxy has a recognizable "core / disc / rim" structure,
  // then sprinkle special regions (nebula/cluster/rift) using the regionSeed.

  const double rxy = std::sqrt(centerLy.x * centerLy.x + centerLy.y * centerLy.y);

  GalaxyRegionKind base = GalaxyRegionKind::OuterRim;
  if (rxy < 8000.0) base = GalaxyRegionKind::Core;
  else if (rxy < 22000.0) base = GalaxyRegionKind::InnerDisc;
  else base = GalaxyRegionKind::OuterRim;

  // Special region roll.
  core::SplitMix64 rng(core::hashCombine(regionSeed, core::fnv1a64("galaxy_region_kind_v1")));
  const double u = rng.nextDouble();

  // Nebulae are more common in the inner disc.
  const double nebulaP = (base == GalaxyRegionKind::InnerDisc) ? 0.14 : 0.07;
  // Clusters are relatively rare but can appear anywhere.
  const double clusterP = 0.06;
  // Rifts are rare and skew toward the outer rim.
  const double riftP = (base == GalaxyRegionKind::OuterRim) ? 0.06 : 0.03;

  if (u < nebulaP) return GalaxyRegionKind::Nebula;
  if (u < nebulaP + clusterP) return GalaxyRegionKind::Cluster;
  if (u < nebulaP + clusterP + riftP) return GalaxyRegionKind::Rift;

  return base;
}

static std::string regionName(core::u64 regionSeed, GalaxyRegionKind kind) {
  core::SplitMix64 rng(core::hashCombine(regionSeed, core::fnv1a64("galaxy_region_name_v1")));

  NameGenerator ng(regionSeed);
  std::string base = ng.systemName();

  // Pick a suffix table per kind.
  const char* suffix = "";
  switch (kind) {
    case GalaxyRegionKind::Core: {
      static const char* s[] = {" Core", " Heart", " Crucible", " Nexus"};
      suffix = s[rng.range<int>(0, 3)];
    } break;
    case GalaxyRegionKind::InnerDisc: {
      static const char* s[] = {" Belt", " Reach", " Expanse", " Marches"};
      suffix = s[rng.range<int>(0, 3)];
    } break;
    case GalaxyRegionKind::OuterRim: {
      static const char* s[] = {" Rim", " Fringe", " Expanse", " Wilds"};
      suffix = s[rng.range<int>(0, 3)];
    } break;
    case GalaxyRegionKind::Nebula: {
      static const char* s[] = {" Nebula", " Veil", " Cloud", " Mist"};
      suffix = s[rng.range<int>(0, 3)];
    } break;
    case GalaxyRegionKind::Cluster: {
      static const char* s[] = {" Cluster", " Swarm", " Knot", " Crown"};
      suffix = s[rng.range<int>(0, 3)];
    } break;
    case GalaxyRegionKind::Rift: {
      static const char* s[] = {" Rift", " Scar", " Maw", " Chasm"};
      suffix = s[rng.range<int>(0, 3)];
    } break;
    default: {
      suffix = " Expanse";
    } break;
  }

  // Occasionally add a leading article for flavor.
  const bool the = rng.nextDouble() < 0.18;
  if (the) {
    return std::string("The ") + base + suffix;
  }
  return base + suffix;
}

} // namespace

const char* galaxyRegionKindName(GalaxyRegionKind k) {
  switch (k) {
    case GalaxyRegionKind::Core: return "Core";
    case GalaxyRegionKind::InnerDisc: return "Inner Disc";
    case GalaxyRegionKind::OuterRim: return "Outer Rim";
    case GalaxyRegionKind::Nebula: return "Nebula";
    case GalaxyRegionKind::Cluster: return "Cluster";
    case GalaxyRegionKind::Rift: return "Rift";
    default: return "Unknown";
  }
}

GalaxyRegionSample sampleGalaxyRegion(core::u64 universeSeed, const math::Vec3d& posLy, double cellSizeLy) {
  GalaxyRegionSample out{};

  cellSizeLy = std::max(1.0, cellSizeLy);

  const core::i64 cx = static_cast<core::i64>(std::floor(posLy.x / cellSizeLy));
  const core::i64 cy = static_cast<core::i64>(std::floor(posLy.y / cellSizeLy));
  const core::i64 cz = static_cast<core::i64>(std::floor(posLy.z / cellSizeLy));

  FeaturePoint best{};
  FeaturePoint second{};
  double bestD2 = std::numeric_limits<double>::infinity();
  double secondD2 = std::numeric_limits<double>::infinity();

  // A 3x3x3 neighborhood is sufficient for classic Worley sampling when each
  // cell contains exactly one feature point.
  for (core::i64 dz = -1; dz <= 1; ++dz) {
    for (core::i64 dy = -1; dy <= 1; ++dy) {
      for (core::i64 dx = -1; dx <= 1; ++dx) {
        const FeaturePoint fp = featurePointForCell(universeSeed, cx + dx, cy + dy, cz + dz, cellSizeLy);
        const math::Vec3d d = fp.p - posLy;
        const double d2 = d.lengthSq();

        if (d2 < bestD2) {
          second = best;
          secondD2 = bestD2;
          best = fp;
          bestD2 = d2;
        } else if (d2 < secondD2) {
          second = fp;
          secondD2 = d2;
        }
      }
    }
  }

  if (!std::isfinite(bestD2)) {
    // Should never happen, but keep it robust.
    return out;
  }

  const double d1 = std::sqrt(bestD2);
  const double d2 = std::isfinite(secondD2) ? std::sqrt(secondD2) : (d1 + cellSizeLy);

  out.regionId = best.id;
  out.regionSeed = core::hashCombine(best.id, core::fnv1a64("galaxy_region_seed_v1"));
  out.centerLy = best.p;
  out.distanceToCenterLy = d1;

  // Heuristic radius: at the feature point (d1=0), d2 approximates distance to the
  // nearest competing region center, so half of that is a reasonable "radius".
  // Away from the center, using 0.5*(d1+d2) stays in the right ballpark.
  out.approxRadiusLy = 0.5 * (d1 + d2);

  // Edge factor using classic Worley f2-f1 logic:
  //  - Near boundary: d2 ~ d1 => edge01 ~ 1
  //  - Deep inside: d2 >> d1 => edge01 ~ 0
  out.edge01 = clamp01(1.0 - (d2 - d1) / std::max(1.0e-9, d2));

  out.kind = classifyRegionKind(out.regionSeed, out.centerLy);
  out.name = regionName(out.regionSeed, out.kind);
  return out;
}

} // namespace stellar::proc
