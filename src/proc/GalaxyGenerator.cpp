#include "stellar/proc/GalaxyGenerator.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/proc/GalaxyVoids.h"
#include "stellar/proc/NameGenerator.h"
#include "stellar/proc/Noise.h"
#include "stellar/proc/GalaxyClusters.h"
#include "stellar/proc/GalaxyMorphology.h"
#include "stellar/proc/TradeProfile.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {

static int tunedStationCount(int baseCount, const TradeProfile& p, core::u64 stubSeed) {
  baseCount = std::max(1, baseCount);

  // Proxy for "economic size" / "market pull".
  // Hub-ness is strongly amplified by population (a hub with no people isn't a hub).
  const double hubPop = std::clamp(p.hub, 0.0, 1.0) * (0.35 + 0.65 * std::clamp(p.population, 0.0, 1.0));
  const double dev = 0.45 * std::clamp(p.industry, 0.0, 1.0)
                   + 0.30 * std::clamp(p.wealth, 0.0, 1.0)
                   + 0.25 * std::clamp(p.technology, 0.0, 1.0);
  const double size01 = std::clamp(0.65 * hubPop + 0.35 * dev, 0.0, 1.0);

  // High lawlessness tends to suppress the number of large installations.
  const double law = std::clamp(p.lawlessness, 0.0, 1.0);

  // Deterministic micro-jitter per stub so not every "0.73 hub" system has
  // exactly the same station count.
  core::SplitMix64 srng(core::hashCombine(stubSeed, core::seedFromText("station_count")));

  int extra = static_cast<int>(std::floor(size01 * 7.0 + srng.nextDouble() * 0.9)); // 0..7
  extra -= static_cast<int>(std::floor(law * 2.2));                                 // 0..2-ish
  extra = std::max(0, extra);

  // Clamp to keep scan costs sane.
  return std::clamp(baseCount + extra, 1, 12);
}

std::size_t SectorCoordHash::operator()(const SectorCoord& c) const noexcept {
  core::u64 h = 0xcbf29ce484222325ull;
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(c.x)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(c.y)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(c.z)));
  return static_cast<std::size_t>(h);
}

GalaxyGenerator::GalaxyGenerator(core::u64 seed, GalaxyParams params)
: seed_(seed), params_(params) {}

SectorCoord GalaxyGenerator::sectorOf(const math::Vec3d& posLy) const {
  const double s = params_.sectorSizeLy;
  return {
    static_cast<core::i32>(std::floor(posLy.x / s)),
    static_cast<core::i32>(std::floor(posLy.y / s)),
    static_cast<core::i32>(std::floor(posLy.z / s)),
  };
}

sim::SystemId GalaxyGenerator::makeSystemId(const SectorCoord& coord, core::u32 localIndex) const {
  // Pack the sector coordinate + local index into 64 bits:
  //   [x:16][y:16][z:16][i:16]
  //
  // This allows Universe::getSystem(id) to decode the sector directly from the id
  // and stream/generate the correct sector on-demand (without requiring a hint stub).
  //
  // NOTE: Galaxy radius/sectorSize keep coords comfortably within +/- 32k sectors
  // (default: 50k ly / 10 ly = 5k sectors), so a 16-bit biased encoding is safe.
  if (coord.x < -32768 || coord.x > 32767 ||
      coord.y < -32768 || coord.y > 32767 ||
      coord.z < -32768 || coord.z > 32767) {
    // Extremely out-of-range query; fall back to a hashed id (won't be decodable).
    core::u64 h = seed_;
    h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.x)));
    h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.y)));
    h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.z)));
    h = core::hashCombine(h, static_cast<core::u64>(localIndex));
    return static_cast<sim::SystemId>(h);
  }

  const auto bias = [](core::i32 v) -> core::u16 {
    return static_cast<core::u16>(static_cast<core::i32>(v) + 32768);
  };

  const core::u64 bx = static_cast<core::u64>(bias(coord.x));
  const core::u64 by = static_cast<core::u64>(bias(coord.y));
  const core::u64 bz = static_cast<core::u64>(bias(coord.z));
  const core::u64 bi = static_cast<core::u64>(static_cast<core::u16>(localIndex & 0xFFFFu));

  const core::u64 id = (bx << 48) | (by << 32) | (bz << 16) | (bi << 0);
  return static_cast<sim::SystemId>(id);
}

static int poisson(core::SplitMix64& rng, double mean) {
  if (mean <= 0.0) return 0;
  if (mean < 30.0) {
    // Knuth
    const double L = std::exp(-mean);
    int k = 0;
    double p = 1.0;
    do {
      ++k;
      p *= rng.nextDouble();
    } while (p > L);
    return k - 1;
  }

  // Normal approximation for larger means (Box-Muller)
  const double u1 = std::max(1e-12, rng.nextDouble());
  const double u2 = rng.nextDouble();
  const double z0 = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * stellar::math::kPi * u2);
  const double k = mean + z0 * std::sqrt(mean);
  return std::max(0, static_cast<int>(std::round(k)));
}

static sim::StarClass pickStarClass(core::SplitMix64& rng) {
  const double r = rng.nextDouble();
  // Very rough main-sequence-ish distribution (not astrophysically accurate).
  if (r < 0.0003) return sim::StarClass::O;
  if (r < 0.0016) return sim::StarClass::B;
  if (r < 0.006)  return sim::StarClass::A;
  if (r < 0.03)   return sim::StarClass::F;
  if (r < 0.10)   return sim::StarClass::G;
  if (r < 0.30)   return sim::StarClass::K;
  return sim::StarClass::M;
}

static core::u32 pickFaction(const math::Vec3d& posLy, const std::vector<sim::Faction>& factions) {
  // 0 is Independent
  core::u32 best = 0;
  double bestD = 1e30;

  for (const auto& f : factions) {
    if (f.id == 0) continue;
    const double d = (posLy - f.homePosLy).length();
    if (d < f.influenceRadiusLy && d < bestD) {
      bestD = d;
      best = f.id;
    }
  }
  return best;
}

static double wrapAnglePi(double a) {
  const double twoPi = 2.0 * stellar::math::kPi;
  a = std::fmod(a + stellar::math::kPi, twoPi);
  if (a < 0.0) a += twoPi;
  return a - stellar::math::kPi;
}

static double fbmAmpSum(int octaves, double gain) {
  double amp = 1.0;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp;
    amp *= gain;
  }
  return sum;
}

static double spiralArmMask(const GalaxyParams& gp, core::u64 noiseSeed, double xLy, double yLy) {
  const int arms = gp.spiralArmCount;
  if (arms <= 0) return 0.0;

  const double rxy = std::sqrt(xLy * xLy + yLy * yLy);
  if (rxy < 1.0e-6) return 0.0;

  const double pitchDeg = std::clamp(gp.spiralPitchDeg, 1.0, 89.0);
  const double pitchRad = pitchDeg * (stellar::math::kPi / 180.0);
  const double k = 1.0 / std::tan(pitchRad);

  const double phaseRad = gp.spiralArmPhaseDeg * (stellar::math::kPi / 180.0);
  const double widthRad = std::max(1.0e-5, gp.spiralArmWidthDeg * (stellar::math::kPi / 180.0));

  // Reference radius for the log-spiral term. Using ~2% of disc radius yields
  // a reasonable number of turns across the galaxy at common pitch angles.
  const double rRef = std::max(1.0, gp.radiusLy * 0.02);
  const double lnTerm = std::log(std::max(1.0e-6, rxy / rRef));

  const double theta = std::atan2(yLy, xLy);

  // Find angular distance to nearest arm.
  double dMin = 1e30;
  const double twoPi = 2.0 * stellar::math::kPi;
  for (int i = 0; i < arms; ++i) {
    const double armTheta = k * lnTerm + phaseRad + twoPi * (static_cast<double>(i) / static_cast<double>(arms));
    const double d = std::abs(wrapAnglePi(theta - armTheta));
    if (d < dMin) dMin = d;
  }

  // Gaussian falloff from the arm centerline.
  const double t = dMin / widthRad;
  double mask = std::exp(-0.5 * t * t);

  // Fade arms out near the very center to avoid excessive winding.
  const double fadeStart = rRef * 0.25;
  const double fadeEnd = rRef * 0.60;
  const double fade = std::clamp((rxy - fadeStart) / std::max(1.0e-6, (fadeEnd - fadeStart)), 0.0, 1.0);
  mask *= fade;

  // Optional modulation along arms.
  const double ns = std::clamp(gp.spiralArmNoiseStrength, 0.0, 1.0);
  if (ns > 0.0 && gp.spiralArmNoiseFreq > 0.0) {
    const double n = smoothNoise2D(noiseSeed, xLy * gp.spiralArmNoiseFreq, yLy * gp.spiralArmNoiseFreq); // 0..1
    const double mul = (1.0 - ns) + (2.0 * ns) * n; // [1-ns, 1+ns]
    mask *= mul;
  }

  return mask;
}

static double densityNoiseMul(const GalaxyParams& gp, core::u64 seed, double xLy, double yLy) {
  const double s = std::clamp(gp.densityNoiseStrength, 0.0, 0.99);
  if (s <= 0.0 || gp.densityNoiseFreq <= 0.0) return 1.0;

  constexpr int kOctaves = 5;
  constexpr double kLacunarity = 2.0;
  constexpr double kGain = 0.5;
  const double ampSum = fbmAmpSum(kOctaves, kGain);

  double n = fbm2D(seed, xLy * gp.densityNoiseFreq, yLy * gp.densityNoiseFreq, kOctaves, kLacunarity, kGain);
  n /= std::max(1.0e-9, ampSum);
  n = std::clamp(n, 0.0, 1.0);

  // Map [0..1] -> [1-s .. 1+s]
  return (1.0 - s) + (2.0 * s) * n;
}

Sector GalaxyGenerator::generateSector(const SectorCoord& coord, const std::vector<sim::Faction>& factions) const {
  // Sector seed
  core::u64 h = seed_;
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.x)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.y)));
  h = core::hashCombine(h, static_cast<core::u64>(static_cast<core::i64>(coord.z)));

  core::SplitMix64 rng(h);

  Sector sector{};
  sector.coord = coord;

  // Density based on distance from galactic center + vertical falloff.
  const double s = params_.sectorSizeLy;
  const math::Vec3d centerLy{
    (coord.x + 0.5) * s,
    (coord.y + 0.5) * s,
    (coord.z + 0.5) * s
  };
  const double rxy = std::sqrt(centerLy.x * centerLy.x + centerLy.y * centerLy.y);

  const double radial = std::exp(-rxy / std::max(1.0, params_.radialScaleLengthLy));

  // Vertical falloff is evaluated relative to an optionally warped midplane and
  // optionally flared thickness.
  const double warpZCenter = galaxyWarpZLy(seed_, params_, centerLy.x, centerLy.y);
  const double halfThicknessCenter = galaxyThicknessHalfLy(params_, rxy);
  const double zRel = std::abs(centerLy.z - warpZCenter);
  const double vertical = std::exp(-zRel / std::max(1.0, halfThicknessCenter));

  const double mean = params_.baseMeanSystemsPerSector * radial * vertical;

  const bool useSpiral = (params_.spiralArmCount > 0) && (params_.spiralArmStrength > 0.0);
  const bool useDensityNoise = (params_.densityNoiseStrength > 0.0);

  const bool useClusters = (params_.clusterStrength > 0.0) &&
                           (params_.clusterCellSizeLy > 0.0) &&
                           (params_.clusterChancePerCell > 0.0) &&
                           (params_.clusterRadiusLy > 0.0);

  const bool useVoids = (params_.voidStrength > 0.0) &&
                        (params_.voidCellSizeLy > 0.0) &&
                        (params_.voidChancePerCell > 0.0) &&
                        (params_.voidRadiusLy > 0.0);

  const bool useBar = (params_.barStrength > 0.0);
  const bool useRing = (params_.ringStrength > 0.0);

  // Disc bounds helper: radius + (optional) warp and flare.
  const auto discContains = [&](const math::Vec3d& p) -> bool {
    const double rr = std::sqrt(p.x * p.x + p.y * p.y);
    if (rr > params_.radiusLy) return false;
    const double wz = galaxyWarpZLy(seed_, params_, p.x, p.y);
    const double halfT = galaxyThicknessHalfLy(params_, rr);
    return std::abs(p.z - wz) <= halfT;
  };

  // ---------------------------------------------------------------------------
  // Optional minimum-separation placement ("blue-noise-ish")
  // ---------------------------------------------------------------------------
  //
  // When enabled (minSystemSeparationLy > 0), we avoid tight clumps by enforcing
  // a minimum distance between neighboring systems.
  //
  // Design goals:
  //  - Deterministic from (seed, galaxy params, global cell coordinate)
  //  - Streaming safe: coherent across sector boundaries without multi-sector state
  //  - Compatible with existing density fields (radial/vertical falloff, spiral arms,
  //    density noise)
  //
  // Approach:
  //  - Define a global 3D grid of "placement cells" with size = r / sqrt(3).
  //  - Each cell contains one jittered candidate point.
  //  - The cell is "enabled" with probability derived from the local intensity.
  //  - Resolve close conflicts by accepting only the highest-priority candidate
  //    within radius r (priority is deterministic per cell).
  //
  // This approximates a Poisson disk / blue-noise distribution and guarantees
  // cross-sector minimum spacing because candidates live in global cell space.

  const double minSepReq = std::max(0.0, params_.minSystemSeparationLy);
  if (minSepReq > 0.0) {
    // Clamp to avoid pathological perf/memory when the separation is extremely small.
    // (At typical galaxy scales, <~0.25 ly separation is visually irrelevant.)
    const double minSep = std::max(0.25, minSepReq);

    const double cellSize = minSep / std::sqrt(3.0);
    const double cellVol = cellSize * cellSize * cellSize;
    const double sectorVol = s * s * s;

    const int nbr = std::max(1, static_cast<int>(std::ceil(minSep / std::max(1.0e-12, cellSize))));

    struct CellCand {
      math::Vec3d pos{};
      double priority{0.0};
      core::u64 tie{0};
      bool enabled{false};
    };

    // Sector bounds (half-open) in world-space ly.
    const double x0 = coord.x * s;
    const double y0 = coord.y * s;
    const double z0 = coord.z * s;
    const double x1 = x0 + s;
    const double y1 = y0 + s;
    const double z1 = z0 + s;

    // Independent noise seeds (so tweaks don't perturb sector seeding).
    const core::u64 spiralNoiseSeed = core::hashCombine(seed_, core::fnv1a64("galaxy_spiral_noise"));
    const core::u64 densityNoiseSeed = core::hashCombine(seed_, core::fnv1a64("galaxy_density_noise"));

    const double armStrength = std::max(0.0, params_.spiralArmStrength);
    const double clusterStrength = std::clamp(params_.clusterStrength, 0.0, 10.0);
    const double voidStrength = std::clamp(params_.voidStrength, 0.0, 10.0);
    const double barStrength = std::clamp(params_.barStrength, 0.0, 10.0);
    const double ringStrength = std::clamp(params_.ringStrength, 0.0, 10.0);

    GalaxyClustersParams clusterParams{};
    if (useClusters) {
      clusterParams.cellSizeLy = params_.clusterCellSizeLy;
      clusterParams.chancePerCell = params_.clusterChancePerCell;
      clusterParams.radiusLy = params_.clusterRadiusLy;
      clusterParams.radiusJitter01 = params_.clusterRadiusJitter;
      clusterParams.strengthJitter01 = params_.clusterStrengthJitter;
      clusterParams.falloffPower = params_.clusterFalloffPower;
    }

    GalaxyVoidsParams voidParams{};
    if (useVoids) {
      voidParams.cellSizeLy = params_.voidCellSizeLy;
      voidParams.chancePerCell = params_.voidChancePerCell;
      voidParams.radiusLy = params_.voidRadiusLy;
      voidParams.radiusJitter01 = params_.voidRadiusJitter;
      voidParams.strengthJitter01 = params_.voidStrengthJitter;
      voidParams.falloffPower = params_.voidFalloffPower;
    }


    auto cellMin = [&](double w0) -> long long {
      return static_cast<long long>(std::floor(w0 / std::max(1.0e-12, cellSize)));
    };
    auto cellMax = [&](double w1) -> long long {
      // Half-open [w0,w1): subtract epsilon so a boundary-aligned cell doesn't get included twice.
      return static_cast<long long>(std::floor((w1 - 1.0e-9) / std::max(1.0e-12, cellSize)));
    };

    const long long minCx0 = cellMin(x0);
    const long long maxCx0 = cellMax(x1);
    const long long minCy0 = cellMin(y0);
    const long long maxCy0 = cellMax(y1);
    const long long minCz0 = cellMin(z0);
    const long long maxCz0 = cellMax(z1);

    // Extend by neighbor range so boundary candidates see cross-sector competitors.
    const long long minCx = minCx0 - nbr;
    const long long maxCx = maxCx0 + nbr;
    const long long minCy = minCy0 - nbr;
    const long long maxCy = maxCy0 + nbr;
    const long long minCz = minCz0 - nbr;
    const long long maxCz = maxCz0 + nbr;

    const std::size_t nx = static_cast<std::size_t>(std::max(0ll, maxCx - minCx + 1));
    const std::size_t ny = static_cast<std::size_t>(std::max(0ll, maxCy - minCy + 1));
    const std::size_t nz = static_cast<std::size_t>(std::max(0ll, maxCz - minCz + 1));

    // Build candidate cache over the extended cell range.
    // This is faster than hashing each neighbor repeatedly.
    std::vector<CellCand> cells;
    cells.resize(nx * ny * nz);

    auto idx = [&](long long cx, long long cy, long long cz) -> std::size_t {
      const std::size_t ix = static_cast<std::size_t>(cx - minCx);
      const std::size_t iy = static_cast<std::size_t>(cy - minCy);
      const std::size_t iz = static_cast<std::size_t>(cz - minCz);
      return ix + nx * (iy + ny * iz);
    };

    auto makeCand = [&](long long cx, long long cy, long long cz) -> CellCand {
      CellCand c{};

      core::u64 cellSeed = seed_;
      cellSeed = core::hashCombine(cellSeed, static_cast<core::u64>(static_cast<core::i64>(cx)));
      cellSeed = core::hashCombine(cellSeed, static_cast<core::u64>(static_cast<core::i64>(cy)));
      cellSeed = core::hashCombine(cellSeed, static_cast<core::u64>(static_cast<core::i64>(cz)));

      c.tie = cellSeed;

      // Jittered candidate position within the cell.
      {
        core::SplitMix64 jrng(core::hashCombine(cellSeed, core::fnv1a64("galaxy_cell_jitter")));
        const double jx = jrng.nextDouble();
        const double jy = jrng.nextDouble();
        const double jz = jrng.nextDouble();
        c.pos = {
          (static_cast<double>(cx) + jx) * cellSize,
          (static_cast<double>(cy) + jy) * cellSize,
          (static_cast<double>(cz) + jz) * cellSize,
        };
      }

      // Priority used for local conflict resolution.
      {
        core::SplitMix64 prng(core::hashCombine(cellSeed, core::fnv1a64("galaxy_cell_priority")));
        c.priority = prng.nextDouble();
      }

      // Hard disc bounds check (includes optional warp + flare).
      const double rr = std::sqrt(c.pos.x * c.pos.x + c.pos.y * c.pos.y);
      if (rr > params_.radiusLy) {
        c.enabled = false;
        return c;
      }

      const double warpZ = galaxyWarpZLy(seed_, params_, c.pos.x, c.pos.y);
      const double halfThickness = galaxyThicknessHalfLy(params_, rr);
      if (std::abs(c.pos.z - warpZ) > halfThickness) {
        c.enabled = false;
        return c;
      }

      // Local intensity at candidate position.
      const double radialLocal = std::exp(-rr / std::max(1.0, params_.radialScaleLengthLy));
      const double verticalLocal = std::exp(-std::abs(c.pos.z - warpZ) / std::max(1.0, halfThickness));
      double meanLocal = params_.baseMeanSystemsPerSector * radialLocal * verticalLocal;

      double mul = 1.0;
      if (useSpiral) {
        mul *= (1.0 + armStrength * spiralArmMask(params_, spiralNoiseSeed, c.pos.x, c.pos.y));
      }
      if (useDensityNoise) {
        mul *= densityNoiseMul(params_, densityNoiseSeed, c.pos.x, c.pos.y);
      }
      if (useClusters) {
        const auto cs = sampleGalaxyClusters(seed_, c.pos, clusterParams);
        mul *= (1.0 + clusterStrength * cs.cluster01);
      }

      if (useVoids) {
        const auto vs = sampleGalaxyVoids(seed_, c.pos, voidParams);
        const double vMul = std::clamp(1.0 - voidStrength * vs.void01, 0.0, 1.0);
        mul *= vMul;
      }
      if (useBar) {
        mul *= (1.0 + barStrength * galaxyBarMask01(params_, c.pos.x, c.pos.y));
      }
      if (useRing) {
        mul *= (1.0 + ringStrength * galaxyRingMask01(params_, rr));
      }
      meanLocal *= mul;

      // Convert expected systems/sector into a per-cell Bernoulli probability.
      double pCell = meanLocal * (cellVol / std::max(1.0e-12, sectorVol));
      pCell = std::clamp(pCell, 0.0, 1.0);

      // Enable with probability pCell.
      core::SplitMix64 erng(core::hashCombine(cellSeed, core::fnv1a64("galaxy_cell_enable")));
      c.enabled = (erng.nextDouble() < pCell);
      return c;
    };

    for (long long cz = minCz; cz <= maxCz; ++cz) {
      for (long long cy = minCy; cy <= maxCy; ++cy) {
        for (long long cx = minCx; cx <= maxCx; ++cx) {
          cells[idx(cx, cy, cz)] = makeCand(cx, cy, cz);
        }
      }
    }

    struct Acc {
      core::u64 tie{0};
      math::Vec3d pos{};
    };

    std::vector<Acc> accepted;
    accepted.reserve(static_cast<std::size_t>(std::max(0.0, mean * 1.25)) + 8);

    const double minSep2 = minSep * minSep;
    const double epsP = 1.0e-12;

    // Accept only local maxima in a neighborhood of radius minSep.
    for (long long cz = minCz0; cz <= maxCz0; ++cz) {
      for (long long cy = minCy0; cy <= maxCy0; ++cy) {
        for (long long cx = minCx0; cx <= maxCx0; ++cx) {
          const CellCand& c = cells[idx(cx, cy, cz)];
          if (!c.enabled) continue;

          // Candidate must actually lie within this sector (cells can straddle sector boundaries).
          if (c.pos.x < x0 || c.pos.x >= x1 || c.pos.y < y0 || c.pos.y >= y1 || c.pos.z < z0 || c.pos.z >= z1) continue;

          bool ok = true;
          for (int dz = -nbr; dz <= nbr && ok; ++dz) {
            for (int dy = -nbr; dy <= nbr && ok; ++dy) {
              for (int dx = -nbr; dx <= nbr; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) continue;

                const long long nxC = cx + dx;
                const long long nyC = cy + dy;
                const long long nzC = cz + dz;

                const CellCand& n = cells[idx(nxC, nyC, nzC)];
                if (!n.enabled) continue;

                const double ddx = n.pos.x - c.pos.x;
                const double ddy = n.pos.y - c.pos.y;
                const double ddz = n.pos.z - c.pos.z;
                const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
                if (d2 >= minSep2) continue;

                // If a neighbor has a strictly higher (priority,tie), we lose.
                if (n.priority > c.priority + epsP) {
                  ok = false;
                  break;
                }
                if (std::abs(n.priority - c.priority) <= epsP && n.tie > c.tie) {
                  ok = false;
                  break;
                }
              }
            }
          }
          if (!ok) continue;

          accepted.push_back(Acc{c.tie, c.pos});

          // The id format encodes localIndex in 16 bits. This should never happen
          // under sane params, but guard anyway.
          if (accepted.size() >= 65535u) goto minsep_done;
        }
      }
    }

minsep_done:

    std::sort(accepted.begin(), accepted.end(), [](const Acc& a, const Acc& b) {
      return a.tie < b.tie;
    });

    sector.systems.reserve(accepted.size());

    NameGenerator ng{};

    for (std::size_t i = 0; i < accepted.size(); ++i) {
      sim::SystemStub stub{};
      stub.id = makeSystemId(coord, static_cast<core::u32>(i));
      stub.seed = core::hashCombine(seed_, static_cast<core::u64>(stub.id));

      ng.reseed(stub.seed);
      stub.name = ng.systemName();
      stub.posLy = accepted[i].pos;

      // Derive per-stub properties from the stub seed to avoid ordering sensitivity.
      core::SplitMix64 srng(core::hashCombine(stub.seed, core::seedFromText("stub_props")));
      stub.primaryClass = pickStarClass(srng);
      stub.planetCount = srng.range(0, 12);

      stub.stationCount = std::max(1, srng.range(0, 3));
      {
        const TradeProfile tp = generateTradeProfile(seed_, stub);
        stub.stationCount = tunedStationCount(stub.stationCount, tp, stub.seed);
      }

      stub.factionId = pickFaction(stub.posLy, factions);

      sector.systems.push_back(std::move(stub));
    }

    // Stable order for deterministic query results.
    std::sort(sector.systems.begin(), sector.systems.end(), [](const sim::SystemStub& a, const sim::SystemStub& b) {
      return a.id < b.id;
    });

    return sector;
  }

  if (!useSpiral && !useDensityNoise && !useClusters && !useVoids && !useBar && !useRing) {
    // ---- Legacy smooth-disc distribution path (keep deterministic signatures stable) ----
    const int n = poisson(rng, mean);

    sector.systems.reserve(static_cast<std::size_t>(std::max(0, n)));

    NameGenerator ng{};
    for (int i = 0; i < n; ++i) {
      // Try a few candidates to keep within disc bounds.
      bool placed = false;
      math::Vec3d pos{};
      for (int tries = 0; tries < 8 && !placed; ++tries) {
        const double ux = rng.nextDouble();
        const double uy = rng.nextDouble();
        const double uz = rng.nextDouble();
        pos = {
          (coord.x + ux) * s,
          (coord.y + uy) * s,
          (coord.z + uz) * s
        };

        if (!discContains(pos)) continue;
        placed = true;
      }
      if (!placed) continue;

      sim::SystemStub stub{};
      stub.id = makeSystemId(coord, static_cast<core::u32>(i));
      stub.seed = core::hashCombine(seed_, static_cast<core::u64>(stub.id));

      ng.reseed(stub.seed);
      stub.name = ng.systemName();
      stub.posLy = pos;
      stub.primaryClass = pickStarClass(rng);
      stub.planetCount = rng.range(0, 12);
      // Keep the legacy RNG draw for determinism, then tune based on the system's
      // macro trade profile.
      stub.stationCount = std::max(1, rng.range(0, 3));
      {
        const TradeProfile tp = generateTradeProfile(seed_, stub);
        stub.stationCount = tunedStationCount(stub.stationCount, tp, stub.seed);
      }
      stub.factionId = pickFaction(stub.posLy, factions);

      sector.systems.push_back(std::move(stub));
    }

    // Stable order for deterministic query results.
    std::sort(sector.systems.begin(), sector.systems.end(), [](const sim::SystemStub& a, const sim::SystemStub& b) {
      return a.id < b.id;
    });

    return sector;
  }

  // ---- Inhomogeneous galaxy path (spiral arms + density noise + clusters) ----
  const double armStrength = std::max(0.0, params_.spiralArmStrength);
  const double armNoiseStrength = std::clamp(params_.spiralArmNoiseStrength, 0.0, 1.0);
  const double densityStrength = std::clamp(params_.densityNoiseStrength, 0.0, 0.99);
  const double clusterStrength = std::clamp(params_.clusterStrength, 0.0, 10.0);
  const double voidStrength = std::clamp(params_.voidStrength, 0.0, 10.0);
  const double barStrength = std::clamp(params_.barStrength, 0.0, 10.0);
  const double ringStrength = std::clamp(params_.ringStrength, 0.0, 10.0);

  GalaxyClustersParams clusterParams{};
  if (useClusters) {
    clusterParams.cellSizeLy = params_.clusterCellSizeLy;
    clusterParams.chancePerCell = params_.clusterChancePerCell;
    clusterParams.radiusLy = params_.clusterRadiusLy;
    clusterParams.radiusJitter01 = params_.clusterRadiusJitter;
    clusterParams.strengthJitter01 = params_.clusterStrengthJitter;
    clusterParams.falloffPower = params_.clusterFalloffPower;
  }

  GalaxyVoidsParams voidParams{};
  if (useVoids) {
    voidParams.cellSizeLy = params_.voidCellSizeLy;
    voidParams.chancePerCell = params_.voidChancePerCell;
    voidParams.radiusLy = params_.voidRadiusLy;
    voidParams.radiusJitter01 = params_.voidRadiusJitter;
    voidParams.strengthJitter01 = params_.voidStrengthJitter;
    voidParams.falloffPower = params_.voidFalloffPower;
  }

  double mulMax = 1.0;
  if (useSpiral) {
    // armMask max is ~ (1+armNoiseStrength)
    mulMax *= (1.0 + armStrength * (1.0 + armNoiseStrength));
  }
  if (useDensityNoise) {
    // densityNoiseMul max is (1+densityStrength)
    mulMax *= (1.0 + densityStrength);
  }
  if (useClusters) {
    // cluster01 is clamped to [0..1], so the max multiplier is (1+clusterStrength)
    mulMax *= (1.0 + clusterStrength);
  }
  // Voids do not increase density (they only suppress), so they don't contribute
  // to mulMax.
  if (useBar) {
    // bar01 is clamped to [0..1], so the max multiplier is (1+barStrength)
    mulMax *= (1.0 + barStrength);
  }
  if (useRing) {
    // ring01 is clamped to [0..1], so the max multiplier is (1+ringStrength)
    mulMax *= (1.0 + ringStrength);
  }

  const double meanMax = mean * mulMax;
  const int nCand = poisson(rng, meanMax);

  sector.systems.reserve(static_cast<std::size_t>(std::max(0, nCand)));

  // Independent noise seeds (so tweaks don't perturb sector seeding).
  const core::u64 spiralNoiseSeed = core::hashCombine(seed_, core::fnv1a64("galaxy_spiral_noise"));
  const core::u64 densityNoiseSeed = core::hashCombine(seed_, core::fnv1a64("galaxy_density_noise"));

  NameGenerator ng{};

  for (int ci = 0; ci < nCand; ++ci) {
    // Try a few candidates to keep within disc bounds.
    bool placed = false;
    math::Vec3d pos{};
    for (int tries = 0; tries < 8 && !placed; ++tries) {
      const double ux = rng.nextDouble();
      const double uy = rng.nextDouble();
      const double uz = rng.nextDouble();
      pos = {
        (coord.x + ux) * s,
        (coord.y + uy) * s,
        (coord.z + uz) * s
      };

      if (!discContains(pos)) continue;
      placed = true;
    }
    if (!placed) continue;

    // Thinning: accept with probability lambda(pos) / lambda_max.
    double mul = 1.0;
    if (useSpiral) {
      mul *= (1.0 + armStrength * spiralArmMask(params_, spiralNoiseSeed, pos.x, pos.y));
    }
    if (useDensityNoise) {
      mul *= densityNoiseMul(params_, densityNoiseSeed, pos.x, pos.y);
    }
    if (useClusters) {
      const auto cs = sampleGalaxyClusters(seed_, pos, clusterParams);
      mul *= (1.0 + clusterStrength * cs.cluster01);
    }

    if (useVoids) {
      const auto vs = sampleGalaxyVoids(seed_, pos, voidParams);
      const double vMul = std::clamp(1.0 - voidStrength * vs.void01, 0.0, 1.0);
      mul *= vMul;
    }
    if (useBar) {
      mul *= (1.0 + barStrength * galaxyBarMask01(params_, pos.x, pos.y));
    }
    if (useRing) {
      const double rr = std::sqrt(pos.x * pos.x + pos.y * pos.y);
      mul *= (1.0 + ringStrength * galaxyRingMask01(params_, rr));
    }

    const double p = std::clamp(mul / std::max(1.0e-9, mulMax), 0.0, 1.0);
    if (rng.nextDouble() > p) continue;

    sim::SystemStub stub{};
    stub.id = makeSystemId(coord, static_cast<core::u32>(sector.systems.size()));
    stub.seed = core::hashCombine(seed_, static_cast<core::u64>(stub.id));

    ng.reseed(stub.seed);
    stub.name = ng.systemName();
    stub.posLy = pos;
    stub.primaryClass = pickStarClass(rng);
    stub.planetCount = rng.range(0, 12);
    // Keep the legacy RNG draw for determinism, then tune based on the system's
    // macro trade profile.
    stub.stationCount = std::max(1, rng.range(0, 3));
    {
      const TradeProfile tp = generateTradeProfile(seed_, stub);
      stub.stationCount = tunedStationCount(stub.stationCount, tp, stub.seed);
    }
    stub.factionId = pickFaction(stub.posLy, factions);

    sector.systems.push_back(std::move(stub));
  }

  // Stable order for deterministic query results.
  std::sort(sector.systems.begin(), sector.systems.end(), [](const sim::SystemStub& a, const sim::SystemStub& b) {
    return a.id < b.id;
  });

  return sector;
}

} // namespace stellar::proc
