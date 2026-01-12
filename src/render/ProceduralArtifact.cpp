#include "stellar/render/ProceduralArtifact.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace stellar::render {

namespace {

struct Blob {
  math::Vec3d c;
  double r{0.0};
};

struct CutPlane {
  math::Vec3d n;
  double h{0.0}; // dot(p,n) <= h is kept
};

static inline double lerpD(double a, double b, double t) {
  return a + (b - a) * t;
}

static inline double smoothMin(double a, double b, double k) {
  // Inigo Quilez smooth min.
  if (k <= 1e-9) return std::min(a, b);
  const double h = std::clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
  return lerpD(b, a, h) - k * h * (1.0 - h);
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

} // namespace

SdfMeshData generateArtifactMesh(core::u64 seed, const ArtifactParams& params) {
  // Pre-generate features deterministically (blobs + cut planes).
  core::SplitMix64 rng(seed);

  const int blobCount = std::clamp(params.blobCount, 0, 64);
  const int cutCount = std::clamp(params.cutCount, 0, 32);

  std::vector<Blob> blobs;
  blobs.reserve(static_cast<std::size_t>(blobCount));

  for (int i = 0; i < blobCount; ++i) {
    const math::Vec3d dir = randomUnitVec(rng);

    const double t = rng.nextDouble();
    const double rad = lerpD((double)params.blobRadiusMin, (double)params.blobRadiusMax, t);

    // Bias centers toward the surface, but keep inside the domain.
    const double centerRadius = (double)params.blobCenterRadius * (0.5 + 0.5 * rng.nextDouble());
    const math::Vec3d c = dir * centerRadius;

    blobs.push_back(Blob{c, rad});
  }

  std::vector<CutPlane> cuts;
  cuts.reserve(static_cast<std::size_t>(cutCount));

  for (int i = 0; i < cutCount; ++i) {
    math::Vec3d n = randomUnitVec(rng);
    n = n.normalized(1e-18);

    const double tt = rng.nextDouble();
    const double h = lerpD((double)params.cutOffsetMin, (double)params.cutOffsetMax, tt);

    cuts.push_back(CutPlane{n, h});
  }

  // Main field.
  const auto field = [=](float x, float y, float z) -> float {
    // Symmetry fold: makes the artifact feel more "engineered".
    const math::Vec3d p{(double)x, (double)y, (double)z};
    const math::Vec3d q{std::abs((double)x), std::abs((double)y), std::abs((double)z)};

    const double r = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z);

    // Base sphere.
    double d = r - (double)params.baseRadius;

    // Noise displacement.
    {
      const double fx = q.x * (double)params.noiseFrequency;
      const double fy = q.y * (double)params.noiseFrequency;
      const double fz = q.z * (double)params.noiseFrequency;

      const double n = proc::fbm3D(seed ^ 0xA11E4F4EULL,
                                  fx,
                                  fy,
                                  fz,
                                  std::clamp(params.noiseOctaves, 1, 10),
                                  std::max(1.01, (double)params.noiseLacunarity),
                                  std::clamp((double)params.noiseGain, 0.01, 0.99));

      // fbm3D is roughly [0,1), remap to [-1,1] so displacement is centered.
      const double nn = n * 2.0 - 1.0;
      d += (double)params.noiseAmplitude * nn;
    }

    // Subtle grooves (trig patterns) layered on top.
    if (params.grooveStrength > 0.0f) {
      const double f = (double)params.grooveFrequency;
      const double gx = std::sin(q.x * f);
      const double gy = std::sin(q.y * f * 0.97);
      const double gz = std::sin(q.z * f * 1.03);
      const double g = gx * gy * gz;
      d += (double)params.grooveStrength * g;
    }

    // Metaball lumps (smooth union with the base).
    {
      const double k = std::max(1e-6, (double)params.blobSmoothK);
      for (const Blob& b : blobs) {
        // Union with a sphere centered at b.c.
        const math::Vec3d dp = p - b.c;
        const double br = std::sqrt(dp.x * dp.x + dp.y * dp.y + dp.z * dp.z) - b.r;
        d = smoothMin(d, br, k);
      }
    }

    // Planar cuts (intersection with half-spaces) -> faceting.
    for (const CutPlane& cp : cuts) {
      const double plane = math::dot(p, cp.n) - cp.h;
      d = std::max(d, plane);
    }

    return (float)d;
  };

  SdfMesherParams mp{};
  mp.resolution = params.resolution;
  mp.bounds = params.bounds;
  mp.iso = params.iso;
  mp.computeNormalsFromField = true;
  mp.normalEps = params.normalEps;
  mp.fixWindingFromNormals = true;

  return meshIsosurfaceMarchingTetrahedra(field, mp);
}

} // namespace stellar::render
