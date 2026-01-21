#include "stellar/render/Nebula.h"

#include "stellar/core/LowDiscrepancy.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

namespace stellar::render {

static constexpr double kTwoPi = 6.283185307179586476925286766559;

static inline void hsvToRgb(float hDeg, float s, float v, float& r, float& g, float& b) {
  // h: [0,360)
  float h = std::fmod(hDeg, 360.0f);
  if (h < 0.0f) h += 360.0f;

  const float c = v * s;
  const float x = c * (1.0f - std::fabs(std::fmod(h / 60.0f, 2.0f) - 1.0f));
  const float m = v - c;

  float rr = 0.0f, gg = 0.0f, bb = 0.0f;
  if (h < 60.0f) {
    rr = c; gg = x; bb = 0.0f;
  } else if (h < 120.0f) {
    rr = x; gg = c; bb = 0.0f;
  } else if (h < 180.0f) {
    rr = 0.0f; gg = c; bb = x;
  } else if (h < 240.0f) {
    rr = 0.0f; gg = x; bb = c;
  } else if (h < 300.0f) {
    rr = x; gg = 0.0f; bb = c;
  } else {
    rr = c; gg = 0.0f; bb = x;
  }

  r = rr + m;
  g = gg + m;
  b = bb + m;
}

namespace {

static inline double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

static inline double smoothstep(double e0, double e1, double x) {
  const double t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3.0 - 2.0 * t);
}

static inline double fbmSumAmp(int octaves, double gain) {
  double amp = 0.5;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp;
    amp *= gain;
  }
  return std::max(1e-9, sum);
}

static inline math::Vec3d dirFromBandUV(double uLon01, double uY01, float bandPower) {
  const double ang = uLon01 * kTwoPi;

  // Start with uniform y in [-1,1], then squeeze toward 0 to create a "galactic plane" band.
  const double y0 = uY01 * 2.0 - 1.0;
  const double ySign = (y0 >= 0.0) ? 1.0 : -1.0;
  const double bp = std::max(1.0, (double)bandPower);
  const double y = ySign * std::pow(std::abs(y0), bp);

  const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
  const double x = r * std::cos(ang);
  const double z = r * std::sin(ang);
  return {x, y, z};
}

struct NebulaDensitySeeds {
  core::u64 warpX{0};
  core::u64 warpY{0};
  core::u64 warpZ{0};
  core::u64 filament{0};
  core::u64 blob{0};
  core::u64 cavity{0};
};

static NebulaDensitySeeds makeDensitySeeds(core::u64 seed) {
  NebulaDensitySeeds s;
  s.warpX = core::hashCombine(seed, core::fnv1a64("nebula_warp_x"));
  s.warpY = core::hashCombine(seed, core::fnv1a64("nebula_warp_y"));
  s.warpZ = core::hashCombine(seed, core::fnv1a64("nebula_warp_z"));
  s.filament = core::hashCombine(seed, core::fnv1a64("nebula_filament"));
  s.blob = core::hashCombine(seed, core::fnv1a64("nebula_blob"));
  s.cavity = core::hashCombine(seed, core::fnv1a64("nebula_cavity"));
  return s;
}

static double nebulaDensity01(const NebulaDensitySeeds& s, const math::Vec3d& dir, double depth01, float bandPower) {
  // Procedural density field for the nebula shell.
  //
  // Goals:
  //  - coherent, filament-like structure (ridged multifractal)
  //  - big soft blobs (low-frequency fBm)
  //  - void bubbles (cavity mask)
  //  - stable + deterministic (seeded), seam-free (evaluated in 3D)

  // Depth maps [0,1] -> [0.35,1.00] so inner/outer shell layers look different.
  const double r = 0.35 + 0.65 * clamp01(depth01);

  // Base sampling coordinates.
  double x = dir.x * r;
  double y = dir.y * r;
  double z = dir.z * r;

  // Domain warp (very low frequency).
  {
    constexpr int kOct = 3;
    constexpr double kGain = 0.55;
    constexpr double kLac = 2.1;
    const double norm = fbmSumAmp(kOct, kGain);

    const double wf = 2.2;     // warp frequency
    const double wa = 0.32;    // warp amplitude

    const double wx = proc::fbmPerlin3D(s.warpX, x * wf, y * wf, z * wf, kOct, kLac, kGain) / norm - 0.5;
    const double wy = proc::fbmPerlin3D(s.warpY, x * wf, y * wf, z * wf, kOct, kLac, kGain) / norm - 0.5;
    const double wz = proc::fbmPerlin3D(s.warpZ, x * wf, y * wf, z * wf, kOct, kLac, kGain) / norm - 0.5;

    x += wx * wa;
    y += wy * wa;
    z += wz * wa;
  }

  // Filament noise (ridged; higher frequency).
  double fil01 = 0.0;
  {
    constexpr int kOct = 5;
    constexpr double kGain = 0.55;
    constexpr double kLac = 2.2;
    const double norm = fbmSumAmp(kOct, kGain);

    const double ff = 8.0;
    const double f = proc::ridgedFbmPerlin3D(s.filament, x * ff, y * ff, z * ff, kOct, kLac, kGain) / norm;
    // Contrast a bit so filaments stand out.
    fil01 = std::pow(clamp01(f), 0.72);
  }

  // Large-scale blobs.
  double blob01 = 0.0;
  {
    constexpr int kOct = 4;
    constexpr double kGain = 0.55;
    constexpr double kLac = 2.0;
    const double norm = fbmSumAmp(kOct, kGain);

    const double bf = 1.65;
    const double b = proc::fbmPerlin3D(s.blob, x * bf, y * bf, z * bf, kOct, kLac, kGain) / norm;
    blob01 = clamp01(b);
    // Sharpen slightly to make big cloud masses.
    blob01 = std::pow(blob01, 1.15);
  }

  // Void bubbles / cavities.
  double cav01 = 0.0;
  {
    const double cf = 3.2;
    const double c = proc::perlin3D(s.cavity, x * cf, y * cf, z * cf);
    // Treat high values as "empty pockets".
    cav01 = smoothstep(0.58, 0.82, c);
  }

  // Plane factor (mild; bandPower already affects candidate distribution).
  const double plane = 1.0 - std::pow(std::abs(dir.y), 0.62);

  // Depth shaping: emphasize mid-shell slightly.
  double d = 1.0 - std::abs(clamp01(depth01) - 0.55) / 0.55;
  d = clamp01(d);
  d = d * d;

  // Combine.
  double dens = (0.16 + 0.84 * fil01) * (0.30 + 0.70 * blob01);
  dens *= (1.0 - 0.78 * cav01);
  dens *= (0.55 + 0.45 * plane);
  dens *= (0.78 + 0.22 * d);

  // Final contrast. (Keep gentle so it's still usable as a probability weight.)
  dens = std::pow(clamp01(dens), 1.25);

  // If caller squeezes the plane very hard, keep density slightly higher overall
  // so we still have enough strong candidates.
  const double bp = std::max(1.0, (double)bandPower);
  dens *= (bp > 2.5) ? 1.10 : 1.0;

  return clamp01(dens);
}

struct Candidate {
  math::Vec3d dir{0,0,1};
  float depth01{0.5f};
  float density01{0.5f};
  float tie01{0.0f};
  std::uint32_t key{0};
};

static inline std::uint32_t quant16(double x01) {
  const int q = (int)std::lround(clamp01(x01) * 65535.0);
  if (q < 0) return 0u;
  if (q > 65535) return 65535u;
  return (std::uint32_t)q;
}

} // namespace

void NebulaField::regenerate(core::u64 seed, int puffCount, float bandPower) {
  seed_ = seed;
  bandPower_ = std::max(1.0f, bandPower);

  if (puffCount < 0) puffCount = 0;

  puffs_.clear();
  puffs_.reserve((std::size_t)puffCount);

  if (puffCount == 0) {
    points_.clear();
    return;
  }

  // We generate a deterministic low-discrepancy candidate set on the unit sphere,
  // evaluate a coherent density field, then keep the best candidates.
  //
  // This produces more nebula-like filaments/voids than sampling puffs uniformly.
  const int candCount = std::max(puffCount * 4, puffCount + 256);

  const NebulaDensitySeeds densitySeeds = makeDensitySeeds(seed_);

  // Low-discrepancy sequence shifts (scramble) so the Halton lattice doesn't show.
  core::SplitMix64 scramble(core::hashCombine(seed_, core::fnv1a64("nebula_qmc")));
  const double sh0 = scramble.nextDouble();
  const double sh1 = scramble.nextDouble();
  const double sh2 = scramble.nextDouble();
  const double sh3 = scramble.nextDouble();

  std::vector<Candidate> cands;
  cands.reserve((std::size_t)candCount);

  for (int i = 0; i < candCount; ++i) {
    const std::uint32_t idx = (std::uint32_t)(i + 1);
    const auto h = core::halton3(idx);

    const double u0 = core::frac(h.x + sh0);
    const double u1 = core::frac(h.y + sh1);
    const double u2 = core::frac(h.z + sh2);
    const double u3 = core::frac(core::halton(idx, 7) + sh3);

    Candidate c{};
    c.dir = dirFromBandUV(u0, u1, bandPower_);

    // Depth selects where between inner/outer radius a puff lives.
    // Bias toward deeper (larger) clouds slightly.
    c.depth01 = (float)std::pow(u2, 0.75);
    c.tie01 = (float)u3;

    const double dens = nebulaDensity01(densitySeeds, c.dir, (double)c.depth01, bandPower_);
    c.density01 = (float)dens;

    // Integer sort key (stable-ish across platforms): 16-bit density + 16-bit tie-breaker.
    const std::uint32_t dq = quant16(dens);
    const std::uint32_t tq = quant16(u3);
    c.key = (dq << 16) | tq;

    cands.push_back(c);
  }

  std::vector<int> order(cands.size());
  std::iota(order.begin(), order.end(), 0);

  auto cmp = [&](int ia, int ib) {
    const auto& a = cands[(std::size_t)ia];
    const auto& b = cands[(std::size_t)ib];
    if (a.key != b.key) return a.key > b.key;
    return ia < ib;
  };

  if ((int)order.size() > puffCount) {
    std::nth_element(order.begin(), order.begin() + puffCount, order.end(), cmp);
    order.resize((std::size_t)puffCount);
  }
  std::sort(order.begin(), order.end(), cmp);

  // Seeded palette fields in HSV space.
  const core::u64 hueSeed = core::hashCombine(seed_, core::fnv1a64("nebula_hue"));
  const core::u64 satSeed = core::hashCombine(seed_, core::fnv1a64("nebula_sat"));
  const core::u64 valSeed = core::hashCombine(seed_, core::fnv1a64("nebula_val"));

  // Normalize the FBM sums to 0..1.
  const double hueNorm = fbmSumAmp(4, 0.55);
  const double satNorm = fbmSumAmp(3, 0.55);
  const double valNorm = fbmSumAmp(3, 0.55);

  for (int idx : order) {
    const Candidate& c = cands[(std::size_t)idx];

    Puff p{};
    p.dir = c.dir;
    p.depth01 = c.depth01;
    p.density01 = c.density01;

    // Color fields evaluated on the same sphere coords so nearby puffs tend to share hues.
    // The palette is biased toward blues/teals/purples with occasional warmer accents.
    const double sx = p.dir.x * (1.15 + 0.45 * (double)p.depth01);
    const double sy = p.dir.y * (1.15 + 0.45 * (double)p.depth01);
    const double sz = p.dir.z * (1.15 + 0.45 * (double)p.depth01);

    const double hn = proc::fbmPerlin3D(hueSeed, sx * 1.4, sy * 1.4, sz * 1.4, 4, 2.0, 0.55) / hueNorm;
    const double sn = proc::fbmPerlin3D(satSeed, sx * 2.2, sy * 2.2, sz * 2.2, 3, 2.1, 0.55) / satNorm;
    const double vn = proc::fbmPerlin3D(valSeed, sx * 2.0, sy * 2.0, sz * 2.0, 3, 2.1, 0.55) / valNorm;

    const float hue = (float)(190.0 + 165.0 * clamp01(hn) + 22.0 * ((double)c.tie01 - 0.5));

    // Saturation: reduce in cavities/void-ish regions; boost in filaments.
    float sat = (float)(0.35 + 0.50 * clamp01(sn));
    sat *= (float)(0.70 + 0.55 * (double)p.density01);
    sat = std::clamp(sat, 0.10f, 0.98f);

    // Value: brighter in denser regions, with some noise.
    float val = (float)(0.45 + 0.50 * clamp01(vn));
    val *= (float)(0.70 + 0.65 * (double)p.density01);
    val = std::clamp(val, 0.10f, 1.00f);

    hsvToRgb(hue, sat, val, p.r, p.g, p.b);

    // Base alpha is density-weighted so dense filaments read stronger.
    const double dens = (double)p.density01;
    p.alpha = (float)std::clamp((0.25 + 0.95 * dens) * (0.55 + 0.55 * (double)c.tie01), 0.08, 1.0);

    // Size: larger puffs in denser zones, but still a bit depth dependent.
    p.size01 = (float)std::clamp((0.25 + 0.55 * dens + 0.20 * (double)p.depth01) * (0.85 + 0.30 * (double)c.tie01), 0.0, 1.0);

    // Turbulence / twinkle.
    // Denser puffs should feel slower / heavier.
    p.twinkleSpeed = (float)(0.12 + 0.88 * (0.20 + 0.80 * (1.0 - dens)) * (0.35 + 0.65 * (double)c.tie01));
    p.phase = (float)((double)c.tie01 * kTwoPi);

    puffs_.push_back(p);
  }

  points_.clear();
  points_.reserve(puffs_.size());

  // Initialize the cached turbulence noise term.
  //
  // NOTE: We intentionally do this once during regeneration; update() will
  // refresh a small batch per frame.
  cachedNoise01_.clear();
  cachedNoise01_.reserve(puffs_.size());
  noiseCursor_ = 0;

  if (!puffs_.empty()) {
    const core::u64 noiseSeed = core::hashCombine(seed_, core::fnv1a64("nebula_update_v2"));
    const double nNorm = fbmSumAmp(3, 0.55);
    constexpr double tf = 4.1;
    constexpr double t0 = 0.0;

    for (const Puff& p : puffs_) {
      const double nx = p.dir.x * tf + t0 * 0.030;
      const double ny = p.dir.y * tf - t0 * 0.022;
      const double nz = p.dir.z * tf + t0 * 0.027;
      const double n = proc::fbmPerlin3D(noiseSeed, nx, ny, nz, 3, 2.1, 0.55) / nNorm;
      cachedNoise01_.push_back((float)n);
    }
  }
}

void NebulaField::update(const math::Vec3d& cameraPosU, double timeSeconds, const Settings& s) {
  points_.clear();
  points_.reserve(puffs_.size());

  const double par = std::clamp(s.parallax, 0.0, 1.0);
  const math::Vec3d anchor = cameraPosU * par;

  const double inner = std::max(0.0, std::min(s.innerRadiusU, s.outerRadiusU));
  const double outer = std::max(s.innerRadiusU, s.outerRadiusU);

  const float sizeMin = std::max(1.0f, s.sizeMinPx);
  const float sizeMax = std::max(sizeMin, s.sizeMaxPx);

  const float opacity = std::clamp(s.opacity, 0.0f, 1.0f);
  const float intensity = std::max(0.0f, s.intensity);

  const float turb = std::clamp(s.turbulence, 0.0f, 1.0f);
  const float turbSpeed = std::max(0.0f, s.turbulenceSpeed);

  const core::u64 noiseSeed = core::hashCombine(seed_, core::fnv1a64("nebula_update_v2"));
  const double nNorm = fbmSumAmp(3, 0.55);

  // Refresh a subset of the cached noise values each frame.
  //
  // This trades a tiny amount of temporal lag for much lower CPU cost vs
  // evaluating fBm noise for every puff every frame.
  if (cachedNoise01_.size() != puffs_.size()) {
    cachedNoise01_.assign(puffs_.size(), 0.5f);
    noiseCursor_ = 0;
  }

  const std::size_t nPuffs = puffs_.size();
  if (nPuffs > 0) {
    // Aim to refresh the whole field over ~24 frames (tunable).
    const int slices = 24;
    const int minBatch = 8;
    const int maxBatch = 256;
    int batch = (int)(nPuffs / (std::size_t)slices);
    batch = std::clamp(batch, minBatch, maxBatch);
    batch = std::min<int>(batch, (int)nPuffs);

    constexpr double tf = 4.1;
    const double t = timeSeconds;
    const float smooth = 0.18f; // exponential smoothing toward new samples

    for (int i = 0; i < batch; ++i) {
      const std::size_t idx = (noiseCursor_ + (std::size_t)i) % nPuffs;
      const Puff& p = puffs_[idx];

      const double nx = p.dir.x * tf + t * 0.030;
      const double ny = p.dir.y * tf - t * 0.022;
      const double nz = p.dir.z * tf + t * 0.027;
      const double ns = proc::fbmPerlin3D(noiseSeed, nx, ny, nz, 3, 2.1, 0.55) / nNorm;

      const float prev = cachedNoise01_[idx];
      const float target = (float)ns;
      cachedNoise01_[idx] = prev + (target - prev) * smooth;
    }

    noiseCursor_ = (noiseCursor_ + (std::size_t)batch) % nPuffs;
  }

  for (std::size_t i = 0; i < puffs_.size(); ++i) {
    const Puff& p = puffs_[i];
    const double r0 = inner + (outer - inner) * (double)p.depth01;

    // 3D turbulence noise sampled on the sphere, animated by drifting coords over time.
    const double n = (i < cachedNoise01_.size()) ? (double)cachedNoise01_[i] : 0.5;

    // Mild per-puff radial variation. Denser puffs jitter slightly less.
    const double jitter = (0.82 + 0.34 * n) * (0.96 - 0.08 * (double)p.density01);
    const double rr = r0 * jitter;

    const math::Vec3d pos = anchor + p.dir * rr;

    // Alpha wobble (turbulence) and density weighting.
    float wobble = 1.0f;
    if (turb > 1e-4f && turbSpeed > 1e-4f) {
      const float w = std::sin((float)(timeSeconds * (double)turbSpeed) * (float)(kTwoPi * (double)p.twinkleSpeed) + p.phase);
      // Denser puffs wobble less.
      wobble = 1.0f + turb * (0.44f - 0.26f * p.density01) * w;
    }

    float a = opacity * p.alpha * wobble;
    // Noise breakup.
    a *= (float)(0.72 + 0.55 * n);
    // Dense filaments appear a bit thicker.
    a *= (float)(0.65 + 0.65 * (double)p.density01);
    a = std::clamp(a, 0.0f, 1.0f);

    // Slightly boost size for dense puffs.
    const float size01 = std::clamp(p.size01 * (0.85f + 0.30f * p.density01), 0.0f, 1.0f);
    const float sz = sizeMin + (sizeMax - sizeMin) * size01;

    PointVertex v{};
    v.px = (float)pos.x;
    v.py = (float)pos.y;
    v.pz = (float)pos.z;
    v.cr = p.r * intensity;
    v.cg = p.g * intensity;
    v.cb = p.b * intensity;
    v.a = a;
    v.size = sz;

    points_.push_back(v);
  }
}

} // namespace stellar::render
