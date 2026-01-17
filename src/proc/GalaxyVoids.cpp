#include "stellar/proc/GalaxyVoids.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {
namespace {

static double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

static double smoothstep01(double t) {
  t = clamp01(t);
  return t * t * (3.0 - 2.0 * t);
}

struct VoidCenter {
  bool present{false};
  core::u64 id{0};
  core::u64 seed{0};
  math::Vec3d centerLy{0,0,0};
  double radiusLy{0.0};
  double intrinsicStrength{1.0};
};

static VoidCenter voidForCell(core::u64 universeSeed,
                              core::i64 cx,
                              core::i64 cy,
                              core::i64 cz,
                              const GalaxyVoidsParams& p) {
  VoidCenter c{};

  const double cell = std::max(1.0, p.cellSizeLy);
  const double chance = std::clamp(p.chancePerCell, 0.0, 1.0);
  if (chance <= 0.0 || !(p.radiusLy > 0.0)) return c;

  // Deterministic per-cell seed.
  core::u64 h = core::hashCombine(universeSeed, core::fnv1a64("galaxy_void_cell_v1"));
  h = core::hashCombine(h, static_cast<core::u64>(cx));
  h = core::hashCombine(h, static_cast<core::u64>(cy));
  h = core::hashCombine(h, static_cast<core::u64>(cz));
  c.seed = h;

  core::SplitMix64 rng(h);

  // Enable void with per-cell probability.
  if (rng.nextDouble() >= chance) return c;

  c.present = true;

  // Jittered center within the cell.
  const math::Vec3d cellMin{(double)cx * cell, (double)cy * cell, (double)cz * cell};
  c.centerLy = cellMin + math::Vec3d{rng.range(0.0, cell), rng.range(0.0, cell), rng.range(0.0, cell)};

  const double rJ = clamp01(p.radiusJitter01);
  const double sJ = clamp01(p.strengthJitter01);

  const double rMul = rng.range(1.0 - rJ, 1.0 + rJ);
  c.radiusLy = std::max(1.0, p.radiusLy * rMul);

  const double sMul = rng.range(1.0 - sJ, 1.0 + sJ);
  c.intrinsicStrength = std::max(0.0, sMul);

  // Stable id intentionally avoids the randomized center coordinates.
  c.id = core::hashCombine(core::hashCombine(universeSeed, core::fnv1a64("galaxy_void_id_v1")), h);
  return c;
}

static double falloff(double dist, double radius, double power) {
  if (!(radius > 0.0)) return 0.0;

  const double t = 1.0 - (dist / radius);
  if (t <= 0.0) return 0.0;

  // Smooth core so the field stays visually pleasing and numerically stable.
  double s = smoothstep01(t);

  // Exponent controls how tight the bubble is.
  const double p = std::max(0.05, power);
  return std::pow(s, p);
}

} // namespace

GalaxyVoidSample sampleGalaxyVoids(core::u64 universeSeed,
                                   const math::Vec3d& posLy,
                                   const GalaxyVoidsParams& params) {
  GalaxyVoidSample out{};

  if (!(params.cellSizeLy > 0.0)) return out;
  if (!(params.radiusLy > 0.0)) return out;
  if (!(params.chancePerCell > 0.0)) return out;

  const double cell = std::max(1.0, params.cellSizeLy);

  const core::i64 cx = static_cast<core::i64>(std::floor(posLy.x / cell));
  const core::i64 cy = static_cast<core::i64>(std::floor(posLy.y / cell));
  const core::i64 cz = static_cast<core::i64>(std::floor(posLy.z / cell));

  // Search enough neighboring cells to cover the maximum possible radius.
  const double rJ = clamp01(params.radiusJitter01);
  const double maxR = params.radiusLy * (1.0 + rJ);

  int search = static_cast<int>(std::ceil(maxR / cell)) + 1;
  search = std::clamp(search, 1, 4); // keep sampling bounded

  double best = 0.0;
  VoidCenter bestC{};

  for (core::i64 dz = -search; dz <= search; ++dz) {
    for (core::i64 dy = -search; dy <= search; ++dy) {
      for (core::i64 dx = -search; dx <= search; ++dx) {
        VoidCenter c = voidForCell(universeSeed, cx + dx, cy + dy, cz + dz, params);
        if (!c.present) continue;

        const math::Vec3d d = posLy - c.centerLy;
        const double dist = std::sqrt(d.lengthSq());
        if (dist > c.radiusLy) continue;

        const double f = falloff(dist, c.radiusLy, params.falloffPower);
        const double influ = f * c.intrinsicStrength;

        if (influ > best) {
          best = influ;
          bestC = c;
        }
      }
    }
  }

  if (best > 0.0 && bestC.present) {
    out.hasVoid = true;
    out.voidId = bestC.id;
    out.voidSeed = bestC.seed;
    out.centerLy = bestC.centerLy;
    out.radiusLy = bestC.radiusLy;
    out.intrinsicStrength = bestC.intrinsicStrength;
    out.void01 = clamp01(best);
  }

  return out;
}

} // namespace stellar::proc
