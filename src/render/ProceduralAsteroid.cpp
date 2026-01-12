#include "stellar/render/ProceduralAsteroid.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>

namespace stellar::render {

static inline double clamp01(double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); }
static inline double lerp(double a, double b, double t) { return a + (b - a) * t; }

static inline double smoothstep(double edge0, double edge1, double x) {
  if (edge0 == edge1) return (x < edge0) ? 0.0 : 1.0;
  const double t = clamp01((x - edge0) / (edge1 - edge0));
  return t * t * (3.0 - 2.0 * t);
}

static inline math::Vec3d randUnitVec(core::SplitMix64& rng) {
  math::Vec3d v{
      rng.nextDouble() * 2.0 - 1.0,
      rng.nextDouble() * 2.0 - 1.0,
      rng.nextDouble() * 2.0 - 1.0
  };
  if (v.lengthSq() < 1e-12) v = {0.0, 1.0, 0.0};
  return v.normalized();
}

static inline double ampSum(int octaves, double gain) {
  double amp = 0.5;
  double s = 0.0;
  for (int i = 0; i < octaves; ++i) {
    s += amp;
    amp *= gain;
  }
  return s;
}

struct Crater {
  math::Vec3d dir;
  double cosRadius{0.0};
  double invSpan{1.0};
  double depth{0.0};
  double rim{0.0};
};

static std::vector<Crater> makeCraters(core::u64 seed, const AsteroidParams& p) {
  std::vector<Crater> cr;
  if (p.craterCount <= 0 || p.craterDepth <= 0.0f) return cr;
  cr.reserve(static_cast<std::size_t>(p.craterCount));

  core::SplitMix64 rng(core::hashCombine(seed, core::fnv1a64("asteroid_craters")));

  const double rMin = math::degToRad((double)std::max(0.1f, p.craterRadiusMinDeg));
  const double rMax = math::degToRad((double)std::max(p.craterRadiusMinDeg, p.craterRadiusMaxDeg));
  const double depthBase = std::max(0.0f, p.craterDepth);
  const double rim = std::clamp((double)p.craterRim, 0.0, 1.0);

  for (int i = 0; i < p.craterCount; ++i) {
    const math::Vec3d d = randUnitVec(rng);
    const double rad = rng.range(rMin, rMax);
    const double cosR = std::cos(rad);
    const double invSpan = 1.0 / std::max(1e-6, (1.0 - cosR));

    // Slight per-crater depth variation keeps the field from looking tiled.
    const double depth = depthBase * rng.range(0.65, 1.15);

    Crater c;
    c.dir = d;
    c.cosRadius = cosR;
    c.invSpan = invSpan;
    c.depth = depth;
    c.rim = rim;
    cr.push_back(c);
  }

  return cr;
}

static void computeNormals(std::vector<VertexPNUT>& v, const std::vector<std::uint32_t>& idx) {
  std::vector<math::Vec3d> acc(v.size(), {0.0, 0.0, 0.0});
  for (std::size_t i = 0; i + 2 < idx.size(); i += 3) {
    const std::uint32_t i0 = idx[i + 0];
    const std::uint32_t i1 = idx[i + 1];
    const std::uint32_t i2 = idx[i + 2];
    if (i0 >= v.size() || i1 >= v.size() || i2 >= v.size()) continue;

    const math::Vec3d p0{v[i0].px, v[i0].py, v[i0].pz};
    const math::Vec3d p1{v[i1].px, v[i1].py, v[i1].pz};
    const math::Vec3d p2{v[i2].px, v[i2].py, v[i2].pz};
    const math::Vec3d n = math::cross(p1 - p0, p2 - p0);

    acc[i0] += n;
    acc[i1] += n;
    acc[i2] += n;
  }

  for (std::size_t i = 0; i < v.size(); ++i) {
    math::Vec3d n = acc[i].normalized();
    if (n.lengthSq() < 1e-12) {
      const math::Vec3d p{v[i].px, v[i].py, v[i].pz};
      n = p.normalized();
    }
    v[i].nx = (float)n.x;
    v[i].ny = (float)n.y;
    v[i].nz = (float)n.z;
  }
}

AsteroidMeshData generateAsteroidMesh(core::u64 seed, const AsteroidParams& params) {
  AsteroidParams p = params;
  p.slices = std::clamp(p.slices, 8, 128);
  p.stacks = std::clamp(p.stacks, 6, 96);
  p.noiseOctaves = std::clamp(p.noiseOctaves, 0, 12);
  p.noiseFrequency = std::max(0.0f, p.noiseFrequency);
  p.noiseAmplitude = std::max(0.0f, p.noiseAmplitude);
  p.baseRadius = std::max(0.01f, p.baseRadius);
  p.noiseGain = std::clamp(p.noiseGain, 0.05f, 0.95f);
  p.noiseLacunarity = std::clamp(p.noiseLacunarity, 1.05f, 4.0f);
  p.minRadius = std::max(0.01f, p.minRadius);
  p.maxRadius = std::max(p.minRadius, p.maxRadius);

  AsteroidMeshData out;
  const std::size_t vCount = static_cast<std::size_t>((p.stacks + 1) * (p.slices + 1));
  out.vertices.resize(vCount);
  out.indices.reserve(static_cast<std::size_t>(p.stacks * p.slices * 6));

  const double twoPi = 2.0 * math::kPi;
  const double pi = math::kPi;

  const core::u64 s0 = core::hashCombine(seed, core::fnv1a64("asteroid_noise_0"));
  const core::u64 s1 = core::hashCombine(seed, core::fnv1a64("asteroid_noise_1"));
  const double freq = (double)p.noiseFrequency;
  const double amp = (double)p.noiseAmplitude;
  const double baseR = (double)p.baseRadius;
  const int oct = p.noiseOctaves;
  const double lac = (double)p.noiseLacunarity;
  const double gain = (double)p.noiseGain;
  const double aSum = std::max(1e-9, ampSum(std::max(1, oct), gain));

  const std::vector<Crater> craters = makeCraters(seed, p);

  auto sample01 = [&](core::u64 s, const math::Vec3d& d, double fMul, double ox, double oy, double oz) {
    if (oct <= 0 || freq <= 0.0) return 0.5;
    const double n = proc::fbm3D(s,
                                d.x * freq * fMul + ox,
                                d.y * freq * fMul + oy,
                                d.z * freq * fMul + oz,
                                oct,
                                lac,
                                gain);
    return clamp01(n / aSum);
  };

  // Build displaced UV sphere vertices.
  for (int i = 0; i <= p.stacks; ++i) {
    const double v = (p.stacks == 0) ? 0.0 : (double)i / (double)p.stacks;
    const double phi = v * pi; // 0..pi
    const double y = std::cos(phi);
    const double r = std::sin(phi);

    for (int j = 0; j <= p.slices; ++j) {
      const double u = (p.slices == 0) ? 0.0 : (double)j / (double)p.slices;
      const double theta = u * twoPi;
      const double x = r * std::cos(theta);
      const double z = r * std::sin(theta);

      math::Vec3d dir{x, y, z};

      // Multi-layer fBm (two seeds) to avoid obvious repetition.
      const double n0 = sample01(s0, dir, 1.00, 0.0, 0.0, 0.0);
      const double n1 = sample01(s1, dir, 0.55, 19.3, -7.1, 3.7);
      const double nSigned = std::clamp((n0 * 2.0 - 1.0) * 0.75 + (n1 * 2.0 - 1.0) * 0.25, -1.0, 1.0);

      double radius = baseR * (1.0 + amp * nSigned);

      // Carve craters (spherical caps).
      if (!craters.empty()) {
        for (const auto& c : craters) {
          const double d = math::dot(dir, c.dir); // cos(angle)
          if (d <= c.cosRadius) continue;
          const double t = clamp01((d - c.cosRadius) * c.invSpan); // 0 at edge -> 1 at center
          const double w = smoothstep(0.0, 1.0, t);

          // Depression at the center.
          radius -= c.depth * w;

          // Small ridge ring near the edge.
          if (c.rim > 0.0) {
            const double ring = smoothstep(0.06, 0.22, t) * (1.0 - smoothstep(0.22, 0.55, t));
            radius += c.depth * c.rim * ring;
          }
        }
      }

      radius = std::clamp(radius, (double)p.minRadius, (double)p.maxRadius);

      const std::size_t idx = static_cast<std::size_t>(i * (p.slices + 1) + j);
      VertexPNUT vert{};
      vert.px = (float)(dir.x * radius);
      vert.py = (float)(dir.y * radius);
      vert.pz = (float)(dir.z * radius);
      // Normals are recomputed after building indices.
      vert.nx = (float)dir.x;
      vert.ny = (float)dir.y;
      vert.nz = (float)dir.z;
      vert.u = (float)u;
      vert.v = (float)v;
      out.vertices[idx] = vert;
    }
  }

  // Indices (two triangles per quad).
  for (int i = 0; i < p.stacks; ++i) {
    for (int j = 0; j < p.slices; ++j) {
      const std::uint32_t a = (std::uint32_t)(i * (p.slices + 1) + j);
      const std::uint32_t b = (std::uint32_t)((i + 1) * (p.slices + 1) + j);
      const std::uint32_t c = (std::uint32_t)((i + 1) * (p.slices + 1) + (j + 1));
      const std::uint32_t d = (std::uint32_t)(i * (p.slices + 1) + (j + 1));

      out.indices.push_back(a);
      out.indices.push_back(b);
      out.indices.push_back(c);

      out.indices.push_back(a);
      out.indices.push_back(c);
      out.indices.push_back(d);
    }
  }

  computeNormals(out.vertices, out.indices);
  return out;
}

AsteroidMeshStats measureAsteroidMesh(const AsteroidMeshData& mesh) {
  AsteroidMeshStats s;
  s.vertexCount = mesh.vertices.size();
  s.indexCount = mesh.indices.size();
  float rMin = 1e30f;
  float rMax = -1e30f;
  for (const auto& v : mesh.vertices) {
    const float r = std::sqrt(v.px * v.px + v.py * v.py + v.pz * v.pz);
    rMin = std::min(rMin, r);
    rMax = std::max(rMax, r);
  }
  if (mesh.vertices.empty()) {
    rMin = 0.0f;
    rMax = 0.0f;
  }
  s.minRadius = rMin;
  s.maxRadius = rMax;
  return s;
}

} // namespace stellar::render
