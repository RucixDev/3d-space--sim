#include "stellar/render/Starfield.h"

#include "stellar/math/Math.h"
#include "stellar/math/Quat.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace stellar::render {
namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;
// Golden angle (radians): pi * (3 - sqrt(5))
constexpr double kGoldenAngle = 2.3999632297286533222315555066336;

inline double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

inline float clamp01f(float v) {
  if (v < 0.0f) return 0.0f;
  if (v > 1.0f) return 1.0f;
  return v;
}

static math::Quatd randomUniformRotation(core::SplitMix64& rng) {
  // Shoemake (1992) uniform random quaternion.
  const double u1 = clamp01(rng.nextUnit());
  const double u2 = clamp01(rng.nextUnit());
  const double u3 = clamp01(rng.nextUnit());

  const double s1 = std::sqrt(std::max(0.0, 1.0 - u1));
  const double s2 = std::sqrt(std::max(0.0, u1));

  const double t1 = kTwoPi * u2;
  const double t2 = kTwoPi * u3;

  const double x = s1 * std::sin(t1);
  const double y = s1 * std::cos(t1);
  const double z = s2 * std::sin(t2);
  const double w = s2 * std::cos(t2);

  return math::Quatd{w, x, y, z}.normalized();
}

static math::Vec3d unitVecOnSphere(core::SplitMix64& rng, double bandPower) {
  // Uniform sampling on a sphere using y in [-1,1] and phi in [0,2pi).
  double y0 = rng.range(-1.0, 1.0);
  const double phi = rng.range(0.0, kTwoPi);

  // Optional band shaping: squeeze toward equator.
  if (bandPower != 1.0) {
    const double s = (y0 >= 0.0) ? 1.0 : -1.0;
    y0 = s * std::pow(std::abs(y0), bandPower);
  }

  const double r = std::sqrt(std::max(0.0, 1.0 - y0 * y0));
  return {r * std::cos(phi), y0, r * std::sin(phi)};
}

static math::Vec3d fibonacciDir(int i, int n, double bandPower) {
  // Deterministic quasi-uniform points on the sphere.
  // Use (i + 0.5)/n so we don't place points directly on the poles.
  const double u = (n > 0) ? ((double)i + 0.5) / (double)n : 0.5;
  double y0 = 1.0 - 2.0 * u;

  // Optional band shaping.
  if (bandPower != 1.0) {
    const double s = (y0 >= 0.0) ? 1.0 : -1.0;
    y0 = s * std::pow(std::abs(y0), bandPower);
  }

  const double r = std::sqrt(std::max(0.0, 1.0 - y0 * y0));
  const double phi = kGoldenAngle * (double)i;
  return {r * std::cos(phi), y0, r * std::sin(phi)};
}

static double expectedSpacingChord(int n) {
  // Heuristic: expected angular spacing for N points on sphere surface is
  // roughly sqrt(4*pi/N). Convert to chord length.
  const double nn = (double)std::max(1, n);
  const double ang = std::sqrt((4.0 * math::kPi) / nn);
  const double chord = 2.0 * std::sin(0.5 * ang);
  return std::clamp(chord, 1.0e-6, 2.0);
}

static void applyRotation(std::vector<math::Vec3d>& dirs, const math::Quatd& q) {
  for (auto& d : dirs) {
    d = q.rotate(d);
  }
}

static void applyTangentJitter(std::vector<math::Vec3d>& dirs,
                               core::SplitMix64& rng,
                               double amplitude) {
  if (amplitude <= 1.0e-12) return;

  for (auto& d : dirs) {
    const math::Vec3d n = d.normalized();

    // Build an orthonormal tangent basis at n.
    const math::Vec3d ref = (std::abs(n.y) < 0.9) ? math::Vec3d{0, 1, 0} : math::Vec3d{1, 0, 0};
    math::Vec3d t = math::cross(ref, n);
    if (t.lengthSq() <= 1.0e-12) {
      t = math::cross(math::Vec3d{0, 0, 1}, n);
    }
    t = t.normalized();
    const math::Vec3d b = math::cross(n, t).normalized();

    const double theta = rng.range(0.0, kTwoPi);
    const double rr = std::sqrt(clamp01(rng.nextUnit())); // uniform disk radius
    const double mag = amplitude * rr;

    const math::Vec3d jitter = t * std::cos(theta) + b * std::sin(theta);
    d = (n + jitter * mag).normalized();
  }
}

static int cellIndex(const math::Vec3d& p, int res) {
  auto toCell = [&](double v) -> int {
    const double t = (v + 1.0) * 0.5; // [-1,1] -> [0,1]
    int i = (int)std::floor(t * (double)res);
    if (i < 0) i = 0;
    if (i >= res) i = res - 1;
    return i;
  };

  const int ix = toCell(p.x);
  const int iy = toCell(p.y);
  const int iz = toCell(p.z);

  return (ix * res + iy) * res + iz;
}

static void relaxBlueNoise(std::vector<math::Vec3d>& dirs,
                           int iterations,
                           double strength,
                           double expectedSpacing,
                           double neighborRadius) {
  const int n = (int)dirs.size();
  if (n <= 1) return;
  if (iterations <= 0) return;
  if (!(strength > 0.0)) return;

  expectedSpacing = std::max(1.0e-6, expectedSpacing);
  neighborRadius = std::max(expectedSpacing, neighborRadius);

  const double neighSq = neighborRadius * neighborRadius;

  // Choose a coarse grid resolution so each bucket holds a small neighbor set.
  // The points live in [-1,1]^3.
  int res = (int)std::ceil(2.0 / neighborRadius);
  res = std::clamp(res, 8, 96);

  const std::size_t cellCount = (std::size_t)res * (std::size_t)res * (std::size_t)res;
  std::vector<int> head(cellCount, -1);
  std::vector<int> next((std::size_t)n, -1);

  std::vector<math::Vec3d> out((std::size_t)n);

  const double moveScale = strength * expectedSpacing * expectedSpacing;
  const double maxMove = expectedSpacing * 0.45;

  for (int it = 0; it < iterations; ++it) {
    std::fill(head.begin(), head.end(), -1);

    // Build linked lists per cell.
    for (int i = 0; i < n; ++i) {
      const int c = cellIndex(dirs[(std::size_t)i], res);
      next[(std::size_t)i] = head[(std::size_t)c];
      head[(std::size_t)c] = i;
    }

    for (int i = 0; i < n; ++i) {
      const math::Vec3d p = dirs[(std::size_t)i].normalized();
      const int base = cellIndex(p, res);

      const int bz = base % res;
      const int by = (base / res) % res;
      const int bx = base / (res * res);

      math::Vec3d force{0, 0, 0};
      int neighborCount = 0;

      // Visit 3x3x3 neighborhood.
      for (int dx = -1; dx <= 1; ++dx) {
        const int x = bx + dx;
        if (x < 0 || x >= res) continue;
        for (int dy = -1; dy <= 1; ++dy) {
          const int y = by + dy;
          if (y < 0 || y >= res) continue;
          for (int dz = -1; dz <= 1; ++dz) {
            const int z = bz + dz;
            if (z < 0 || z >= res) continue;

            const int cell = (x * res + y) * res + z;
            for (int j = head[(std::size_t)cell]; j != -1; j = next[(std::size_t)j]) {
              if (j == i) continue;
              const math::Vec3d q = dirs[(std::size_t)j];
              const math::Vec3d d = p - q;
              const double dsq = d.lengthSq();
              if (dsq <= 1.0e-12 || dsq >= neighSq) continue;

              // Repulsive force, stronger at small distances.
              force += d * (1.0 / (dsq + 1.0e-9));
              ++neighborCount;
            }
          }
        }
      }

      if (neighborCount > 0) {
        force *= (1.0 / (double)neighborCount);

        // Project to tangent plane of the sphere so we don't drift radially.
        force -= p * math::dot(force, p);

        math::Vec3d move = force * moveScale;
        const double ml = move.length();
        if (ml > maxMove && ml > 1.0e-12) {
          move *= (maxMove / ml);
        }

        out[(std::size_t)i] = (p + move).normalized();
      } else {
        out[(std::size_t)i] = p;
      }
    }

    dirs.swap(out);

    // Slightly reduce move scale over time for better convergence.
    // (We don't change `strength` itself to keep deterministic behavior stable.)
  }
}

static void kelvinToRgb(double kelvin, float& r, float& g, float& b) {
  // Approximate blackbody color temperature -> sRGB.
  // Input range clamped to common approximation domain.
  const double k = std::clamp(kelvin, 1000.0, 40000.0);
  const double t = k / 100.0;

  double rr = 0.0;
  double gg = 0.0;
  double bb = 0.0;

  if (t <= 66.0) {
    rr = 255.0;
    gg = 99.4708025861 * std::log(t) - 161.1195681661;
    if (t <= 19.0) {
      bb = 0.0;
    } else {
      bb = 138.5177312231 * std::log(t - 10.0) - 305.0447927307;
    }
  } else {
    rr = 329.698727446 * std::pow(t - 60.0, -0.1332047592);
    gg = 288.1221695283 * std::pow(t - 60.0, -0.0755148492);
    bb = 255.0;
  }

  r = (float)std::clamp(rr / 255.0, 0.0, 1.0);
  g = (float)std::clamp(gg / 255.0, 0.0, 1.0);
  b = (float)std::clamp(bb / 255.0, 0.0, 1.0);
}

static void pickStarColor(core::SplitMix64& rng,
                          double bright01,
                          float& r,
                          float& g,
                          float& b,
                          float& lumBoost) {
  // We bias spectral types based on brightness so the few bright stars are more
  // likely to be hotter/whiter (closer to what humans perceive in a night sky).
  const double t = std::pow(clamp01(bright01), 1.35);

  auto lerp = [](double a, double b, double t) { return a + (b - a) * t; };

  // Dim-heavy distribution (dominated by red dwarfs).
  const double dm = 0.72, dk = 0.18, dg = 0.07, df = 0.025, da = 0.004, dbt = 0.001, doo = 0.0002;
  // Bright-weighted distribution (more white/blue for visible bright points).
  const double bm = 0.25, bk = 0.20, bg = 0.18, bf = 0.14,  ba = 0.12,  bbt = 0.07,  boo = 0.04;

  double wM = lerp(dm, bm, t);
  double wK = lerp(dk, bk, t);
  double wG = lerp(dg, bg, t);
  double wF = lerp(df, bf, t);
  double wA = lerp(da, ba, t);
  double wB = lerp(dbt, bbt, t);
  double wO = lerp(doo, boo, t);

  const double sum = std::max(1.0e-12, wM + wK + wG + wF + wA + wB + wO);
  wM /= sum; wK /= sum; wG /= sum; wF /= sum; wA /= sum; wB /= sum; wO /= sum;

  const double u = clamp01(rng.nextUnit());

  enum class Type { M, K, G, F, A, B, O };
  Type type = Type::M;

  double c = 0.0;
  c += wM; if (u < c) type = Type::M;
  else { c += wK; if (u < c) type = Type::K;
  else { c += wG; if (u < c) type = Type::G;
  else { c += wF; if (u < c) type = Type::F;
  else { c += wA; if (u < c) type = Type::A;
  else { c += wB; if (u < c) type = Type::B;
  else type = Type::O; }}}}}

  double kelvin = 5200.0;
  lumBoost = 1.0f;

  switch (type) {
    case Type::M:
      kelvin = rng.range(2400.0, 3700.0);
      lumBoost = 0.90f;
      break;
    case Type::K:
      kelvin = rng.range(3700.0, 5200.0);
      lumBoost = 0.96f;
      break;
    case Type::G:
      kelvin = rng.range(5200.0, 6000.0);
      lumBoost = 1.00f;
      break;
    case Type::F:
      kelvin = rng.range(6000.0, 7500.0);
      lumBoost = 1.04f;
      break;
    case Type::A:
      kelvin = rng.range(7500.0, 10000.0);
      lumBoost = 1.10f;
      break;
    case Type::B:
      kelvin = rng.range(10000.0, 25000.0);
      lumBoost = 1.18f;
      break;
    case Type::O:
      kelvin = rng.range(25000.0, 40000.0);
      lumBoost = 1.25f;
      break;
  }

  kelvinToRgb(kelvin, r, g, b);

  // Gentle per-star variance.
  const float v = (float)rng.range(0.93, 1.07);
  r = clamp01f(r * v);
  g = clamp01f(g * v);
  b = clamp01f(b * v);

  // Subtle dust reddening for some stars (more noticeable on the dim ones).
  // This helps avoid an overly-clean palette.
  const float dust = (float)(std::pow(rng.nextUnit(), 1.6) * (0.10 + 0.20 * (1.0 - bright01)));
  b = clamp01f(b * (1.0f - 0.55f * dust));
  g = clamp01f(g * (1.0f - 0.20f * dust));
}

} // namespace

void Starfield::regenerate(core::u64 seed, int starCount) {
  regenerate(seed, starCount, settings_);
}

void Starfield::regenerate(core::u64 seed, int starCount, const Settings& settings) {
  settings_ = settings;
  seed_ = seed;

  stars_.clear();
  points_.clear();

  starCount = std::clamp(starCount, 0, 200000);
  if (starCount <= 0) return;

  stars_.reserve((std::size_t)starCount);

  // Independent streams keep style deterministic even when placement settings change.
  core::SplitMix64 rngDirs(core::hashCombine(seed, core::fnv1a64("starfield_dirs")));
  core::SplitMix64 rngStyle(core::hashCombine(seed, core::fnv1a64("starfield_style")));

  const double bandPower = std::max(0.05, settings_.bandPower);

  std::vector<math::Vec3d> dirs;
  dirs.reserve((std::size_t)starCount);

  switch (settings_.distribution) {
    case Distribution::UniformRandom:
      for (int i = 0; i < starCount; ++i) {
        dirs.push_back(unitVecOnSphere(rngDirs, bandPower));
      }
      break;
    case Distribution::Fibonacci:
    case Distribution::RelaxedFibonacci:
      for (int i = 0; i < starCount; ++i) {
        dirs.push_back(fibonacciDir(i, starCount, bandPower));
      }
      break;
  }

  // Randomize overall orientation to avoid obvious lattice alignment.
  if (settings_.randomRotation && starCount >= 3) {
    applyRotation(dirs, randomUniformRotation(rngDirs));
  }

  // Jitter is specified in terms of expected spacing.
  const double spacing = expectedSpacingChord(starCount);
  if (settings_.jitter01 > 0.0) {
    const double amp = std::clamp(settings_.jitter01, 0.0, 1.0) * spacing;
    applyTangentJitter(dirs, rngDirs, amp);
  }

  // Optional lightweight blue-noise relaxation.
  if (settings_.distribution == Distribution::RelaxedFibonacci) {
    const int iters = std::clamp(settings_.relaxIterations, 0, 16);
    const double str = std::clamp(settings_.relaxStrength, 0.0, 1.0);
    if (iters > 0 && str > 0.0 && starCount >= 16) {
      relaxBlueNoise(dirs, iters, str, spacing, spacing * 3.2);
    }
  }

  // Materialize star attributes.
  for (int i = 0; i < starCount; ++i) {
    Star s;
    s.dir = dirs[(std::size_t)i].normalized();

    // Magnitude-ish distribution: many dim, few bright.
    // We keep this separate from spectral type selection.
    const double m = std::pow(rngStyle.nextUnit(), 2.8); // skew toward 0
    const double bright01 = std::clamp(1.0 - m, 0.0, 1.0);

    float lumBoost = 1.0f;
    pickStarColor(rngStyle, bright01, s.r, s.g, s.b, lumBoost);

    // baseAlpha roughly [0.08 .. 1.0]
    s.baseAlpha = (float)std::clamp(0.08 + bright01 * 0.92, 0.02, 1.0);
    s.baseAlpha = std::clamp(s.baseAlpha * lumBoost, 0.02f, 1.0f);

    // Slight size correlation with brightness + stellar type.
    double baseSize = 1.0 + bright01 * 2.6;
    baseSize *= (double)std::sqrt(std::max(0.30f, lumBoost));
    s.sizePx = (float)std::clamp(baseSize + rngStyle.range(-0.25, 0.35), 0.75, 4.25);

    // Twinkle: slow for dim, slightly faster for bright.
    s.twinkleSpeed = (float)std::clamp(0.35 + bright01 * 1.8, 0.25, 2.5);
    s.phase = (float)rngStyle.range(0.0, kTwoPi);

    stars_.push_back(s);
  }
}

void Starfield::update(const math::Vec3d& cameraPosU, double timeSeconds) {
  points_.clear();
  points_.reserve(stars_.size());

  const double t = timeSeconds;

  for (const auto& s : stars_) {
    const math::Vec3d posU = cameraPosU + s.dir * radiusU_;

    const float tw = 0.85f + 0.15f * (float)std::sin(t * (double)s.twinkleSpeed + (double)s.phase);
    const float a = std::clamp(s.baseAlpha * tw, 0.0f, 1.0f);

    PointVertex v;
    v.px = (float)posU.x;
    v.py = (float)posU.y;
    v.pz = (float)posU.z;
    v.cr = s.r;
    v.cg = s.g;
    v.cb = s.b;
    v.a = a;
    v.size = s.sizePx;
    points_.push_back(v);
  }
}

void Starfield::exportGpuStars(std::vector<GpuStar>& out) const {
  out.clear();
  out.reserve(stars_.size());
  for (const auto& s : stars_) {
    GpuStar g;
    g.dx = (float)s.dir.x;
    g.dy = (float)s.dir.y;
    g.dz = (float)s.dir.z;
    g.r = s.r;
    g.g = s.g;
    g.b = s.b;
    g.baseAlpha = s.baseAlpha;
    g.sizePx = s.sizePx;
    g.twinkleSpeed = s.twinkleSpeed;
    g.phase = s.phase;
    out.push_back(g);
  }
}

} // namespace stellar::render
