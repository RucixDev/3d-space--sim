#include "stellar/core/Hash.h"
#include "stellar/render/ProceduralAsteroid.h"

#include <cmath>
#include <cstdint>
#include <iostream>

static stellar::core::u64 hashMeshSignature(const stellar::render::AsteroidMeshData& m) {
  stellar::core::u64 h = stellar::core::fnv1a64("asteroid_signature");

  auto mixQ = [&](double v) {
    // Quantize to reduce platform-specific FP drift (while still catching logic changes).
    const std::int64_t q = (std::int64_t)std::llround(v * 100000.0);
    h = stellar::core::hashCombine(h, stellar::core::hashBytes(&q, sizeof(q)));
  };

  const std::size_t n = std::min<std::size_t>(m.vertices.size(), 96);
  for (std::size_t i = 0; i < n; ++i) {
    const auto& v = m.vertices[i];
    mixQ(v.px); mixQ(v.py); mixQ(v.pz);
    mixQ(v.nx); mixQ(v.ny); mixQ(v.nz);
    mixQ(v.u);  mixQ(v.v);
  }

  const std::size_t ni = std::min<std::size_t>(m.indices.size(), 96);
  for (std::size_t i = 0; i < ni; ++i) {
    const std::uint32_t ix = m.indices[i];
    h = stellar::core::hashCombine(h, (stellar::core::u64)ix);
  }
  return h;
}

int test_procedural_asteroid() {
  int fails = 0;

  const stellar::core::u64 seed = 424242;

  stellar::render::AsteroidParams p;
  p.slices = 16;
  p.stacks = 10;
  p.noiseFrequency = 2.8f;
  p.noiseAmplitude = 0.25f;
  p.noiseOctaves = 4;
  p.craterCount = 10;
  p.craterDepth = 0.11f;
  p.craterRim = 0.22f;
  p.minRadius = 0.60f;
  p.maxRadius = 1.40f;

  const auto m0 = stellar::render::generateAsteroidMesh(seed, p);
  const auto m1 = stellar::render::generateAsteroidMesh(seed, p);

  const std::size_t expectedV = (std::size_t)((p.stacks + 1) * (p.slices + 1));
  const std::size_t expectedI = (std::size_t)(p.stacks * p.slices * 6);
  if (m0.vertices.size() != expectedV) {
    std::cerr << "[test_procedural_asteroid] vertex count mismatch " << m0.vertices.size() << " vs " << expectedV << "\n";
    ++fails;
  }
  if (m0.indices.size() != expectedI) {
    std::cerr << "[test_procedural_asteroid] index count mismatch " << m0.indices.size() << " vs " << expectedI << "\n";
    ++fails;
  }

  // Determinism signature.
  const auto h0 = hashMeshSignature(m0);
  const auto h1 = hashMeshSignature(m1);
  if (h0 != h1) {
    std::cerr << "[test_procedural_asteroid] signature mismatch (non-deterministic)\n";
    ++fails;
  }

  // Normals should be finite and close to unit length.
  for (std::size_t i = 0; i < m0.vertices.size(); i += 17) {
    const auto& v = m0.vertices[i];
    const double nl2 = (double)v.nx * v.nx + (double)v.ny * v.ny + (double)v.nz * v.nz;
    if (!std::isfinite(nl2) || nl2 < 0.25 || nl2 > 2.25) {
      std::cerr << "[test_procedural_asteroid] bad normal length^2=" << nl2 << " at i=" << i << "\n";
      ++fails;
      break;
    }
  }

  // Radius clamp sanity.
  const auto st = stellar::render::measureAsteroidMesh(m0);
  if (st.minRadius < p.minRadius - 1e-3f || st.maxRadius > p.maxRadius + 1e-3f) {
    std::cerr << "[test_procedural_asteroid] radius clamp violated "
              << st.minRadius << ".." << st.maxRadius
              << " expected " << p.minRadius << ".." << p.maxRadius << "\n";
    ++fails;
  }

  if (fails == 0) std::cout << "[test_procedural_asteroid] pass\n";
  return fails;
}
