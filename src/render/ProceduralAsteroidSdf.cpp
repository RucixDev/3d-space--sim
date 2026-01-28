#include "stellar/render/ProceduralAsteroidSdf.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace stellar::render {

namespace {

struct Crater {
  math::Vec3d c;
  double r{0.0};
  double smoothK{0.0};
};

struct CutPlane {
  math::Vec3d n;
  double h{0.0}; // keep dot(p,n) <= h
};

static inline double lerpD(double a, double b, double t) { return a + (b - a) * t; }

// Inigo Quilez smooth min.
static inline double smoothMin(double a, double b, double k) {
  if (k <= 1e-12) return std::min(a, b);
  const double h = std::clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
  return lerpD(b, a, h) - k * h * (1.0 - h);
}

static inline double smoothMax(double a, double b, double k) {
  // max(a,b) == -min(-a,-b)
  return -smoothMin(-a, -b, k);
}

static inline math::Vec3d randomUnitVec(core::SplitMix64& rng) {
  // Uniform direction on a sphere.
  const double u = rng.nextDouble();
  const double v = rng.nextDouble();
  const double z = 2.0 * u - 1.0;
  const double a = 2.0 * math::kPi * v;
  const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
  return math::Vec3d{r * std::cos(a), z, r * std::sin(a)};
}

static inline double clampd(double v, double lo, double hi) {
  return std::max(lo, std::min(v, hi));
}

static inline double degToRad(double deg) {
  return deg * (math::kPi / 180.0);
}

} // namespace


SdfMeshData generateAsteroidSdfMesh(core::u64 seed, const AsteroidSdfParams& params) {
  core::SplitMix64 rng(seed);

  const int craterCount = std::clamp(params.craterCount, 0, 48);
  const int cutCount = std::clamp(params.cutCount, 0, 16);

  const double baseR = std::max(1e-6, (double)params.baseRadius);

  // --- Feature pre-generation (deterministic) ---
  std::vector<Crater> craters;
  craters.reserve((std::size_t)craterCount);

  const double rMin = std::max(0.5, (double)params.craterRadiusMinDeg);
  const double rMax = std::max(rMin, (double)params.craterRadiusMaxDeg);
  const double smoothK = std::max(0.0, (double)params.craterSmoothK);
  const double depthFrac = clampd((double)params.craterDepth, 0.0, 1.0);

  for (int i = 0; i < craterCount; ++i) {
    const math::Vec3d dir = randomUnitVec(rng).normalized(1e-18);

    // Convert a crater angular radius to an approximate world-space radius.
    // For small angles: r ~= R * sin(theta).
    const double t = rng.nextDouble();
    const double ang = lerpD(rMin, rMax, t);
    const double theta = degToRad(ang);
    const double r = baseR * std::sin(theta);

    // Place the subtracting sphere slightly below the surface so we carve an indentation.
    const double depth = depthFrac * r;
    const math::Vec3d c = dir * (baseR - depth);

    Crater cr;
    cr.c = c;
    cr.r = std::max(1e-6, r);
    cr.smoothK = smoothK;
    craters.push_back(cr);
  }

  std::vector<CutPlane> cuts;
  cuts.reserve((std::size_t)cutCount);

  const double hMin = std::max(0.0, (double)params.cutOffsetMin);
  const double hMax = std::max(hMin, (double)params.cutOffsetMax);

  for (int i = 0; i < cutCount; ++i) {
    math::Vec3d n = randomUnitVec(rng).normalized(1e-18);

    // Bias planes to be roughly near the surface but not through the center.
    const double t = rng.nextDouble();
    const double h = lerpD(hMin, hMax, t) * baseR;
    cuts.push_back(CutPlane{n, h});
  }

  // --- Scalar field ---
  const double ax = std::max(1e-4, (double)params.axisScaleX);
  const double ay = std::max(1e-4, (double)params.axisScaleY);
  const double az = std::max(1e-4, (double)params.axisScaleZ);

  const int oct1 = std::clamp(params.noise1Octaves, 1, 12);
  const int oct2 = std::clamp(params.noise2Octaves, 1, 12);

  const double n1F = std::max(0.0, (double)params.noise1Frequency);
  const double n1A = (double)params.noise1Amplitude;
  const double n1L = std::max(1.01, (double)params.noise1Lacunarity);
  const double n1G = clampd((double)params.noise1Gain, 0.01, 0.99);

  const double n2F = std::max(0.0, (double)params.noise2Frequency);
  const double n2A = (double)params.noise2Amplitude;
  const double n2L = std::max(1.01, (double)params.noise2Lacunarity);
  const double n2G = clampd((double)params.noise2Gain, 0.01, 0.99);

  const double gStr = std::max(0.0, (double)params.grooveStrength);
  const double gFreq = std::max(0.0, (double)params.grooveFrequency);

  const core::u64 seed1 = core::hashCombine(seed, core::fnv1a64("ast_sdf_noise1"));
  const core::u64 seed2 = core::hashCombine(seed, core::fnv1a64("ast_sdf_noise2"));

  const auto field = [=](float x, float y, float z) -> float {
    // Domain scale for elongation. Keep the output mesh in roughly the same bounds.
    const math::Vec3d p{(double)x, (double)y, (double)z};
    const math::Vec3d q{p.x / ax, p.y / ay, p.z / az};

    // Base: sphere.
    double d = q.length() - baseR;

    // Crater carving: SDF subtraction (max(dA, -dB)).
    for (const Crater& c : craters) {
      const double ds = (q - c.c).length() - c.r;
      if (c.smoothK > 1e-9) {
        d = smoothMax(d, -ds, c.smoothK);
      } else {
        d = std::max(d, -ds);
      }
    }

    // Faceting: intersect with half-spaces.
    for (const CutPlane& cp : cuts) {
      const double plane = math::dot(q, cp.n) - cp.h;
      d = std::max(d, plane);
    }

    // Two-stage noise displacement.
    if (n1A != 0.0 && n1F > 0.0) {
      const double n = proc::fbm3D(seed1,
                                  q.x * n1F,
                                  q.y * n1F,
                                  q.z * n1F,
                                  oct1,
                                  n1L,
                                  n1G);
      d += n1A * (n * 2.0 - 1.0);
    }
    if (n2A != 0.0 && n2F > 0.0) {
      // Higher-frequency Perlin fbm tends to look less "grid" at small scales.
      const double n = proc::fbmPerlin3D(seed2,
                                        q.x * n2F,
                                        q.y * n2F,
                                        q.z * n2F,
                                        oct2,
                                        n2L,
                                        n2G);
      d += n2A * (n * 2.0 - 1.0);
    }

    // Subtle trig grooves (small engineered striations).
    if (gStr > 1e-8 && gFreq > 0.0) {
      const double gx = std::sin(q.x * gFreq);
      const double gy = std::sin(q.y * gFreq * 0.97);
      const double gz = std::sin(q.z * gFreq * 1.03);
      d += gStr * (gx * gy * gz);
    }

    // Undo axis scaling on distance to keep the shape close to the requested radius.
    // This is an approximation, but works well for small anisotropy.
    const double s = (ax + ay + az) / 3.0;
    return (float)(d * s);
  };

  SdfMesherParams mp{};
  mp.resolution = std::clamp(params.resolution, 12, 192);
  mp.bounds = std::max(0.25f, params.bounds);
  mp.iso = params.iso;
  mp.computeNormalsFromField = true;
  mp.normalEps = std::max(1e-5f, params.normalEps);
  mp.fixWindingFromNormals = true;

  return meshIsosurfaceMarchingTetrahedra(field, mp);
}

AsteroidMeshStats measureAsteroidSdfMesh(const SdfMeshData& mesh) {
  AsteroidMeshStats st;
  st.vertexCount = mesh.vertices.size();
  st.indexCount = mesh.indices.size();

  double rMin = std::numeric_limits<double>::infinity();
  double rMax = 0.0;

  for (const auto& v : mesh.vertices) {
    const double r = std::sqrt((double)v.px * (double)v.px + (double)v.py * (double)v.py + (double)v.pz * (double)v.pz);
    rMin = std::min(rMin, r);
    rMax = std::max(rMax, r);
  }

  if (!std::isfinite(rMin)) rMin = 0.0;
  st.minRadius = (float)rMin;
  st.maxRadius = (float)rMax;
  return st;
}

} // namespace stellar::render
