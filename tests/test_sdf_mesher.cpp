#include "test_harness.h"

#include "stellar/render/SdfMesher.h"
#include "stellar/render/ProceduralArtifact.h"

#include <algorithm>
#include <cmath>

using namespace stellar;

int test_sdf_mesher() {
  int failures = 0;

  // --- Marching tetrahedra: sphere SDF sanity ---
  {
    const render::ScalarField3D sphere = [](float x, float y, float z) -> float {
      const double r = std::sqrt((double)x * x + (double)y * y + (double)z * z);
      return (float)(r - 0.50); // radius 0.5
    };

    render::SdfMesherParams mp{};
    mp.resolution = 24;
    mp.bounds = 0.75f;
    mp.iso = 0.0f;
    mp.normalEps = 0.0025f;

    const auto m1 = render::meshIsosurfaceMarchingTetrahedra(sphere, mp);
    const auto s1 = render::measureSdfMesh(m1);

    CHECK(!m1.vertices.empty());
    CHECK(!m1.indices.empty());
    CHECK(m1.indices.size() % 3u == 0u);
    CHECK(s1.vertexCount == m1.vertices.size());
    CHECK(s1.indexCount == m1.indices.size());

    // Bounding box should roughly match the sphere radius.
    CHECK(s1.minX < -0.45f);
    CHECK(s1.maxX > 0.45f);
    CHECK(s1.minY < -0.45f);
    CHECK(s1.maxY > 0.45f);
    CHECK(s1.minZ < -0.45f);
    CHECK(s1.maxZ > 0.45f);

    // Normals should be reasonably normalized.
    const std::size_t sampleN = std::min<std::size_t>(m1.vertices.size(), 128);
    for (std::size_t i = 0; i < sampleN; ++i) {
      const auto& v = m1.vertices[i];
      const double len = std::sqrt((double)v.nx * v.nx + (double)v.ny * v.ny + (double)v.nz * v.nz);
      CHECK(len > 0.5);
      CHECK(len < 1.5);
    }

    // Determinism: same field + params => same topology sizes.
    const auto m2 = render::meshIsosurfaceMarchingTetrahedra(sphere, mp);
    CHECK(m2.vertices.size() == m1.vertices.size());
    CHECK(m2.indices.size() == m1.indices.size());
  }

  // --- Procedural artifact generator: determinism + variety ---
  {
    render::ArtifactParams ap{};
    ap.resolution = 32;
    ap.bounds = 1.35f;
    ap.baseRadius = 0.90f;
    ap.noiseAmplitude = 0.25f;
    ap.blobCount = 5;
    ap.cutCount = 4;
    ap.grooveStrength = 0.08f;

    const core::u64 seedA = 1234567ull;
    const core::u64 seedB = 1234568ull;

    const auto a1 = render::generateArtifactMesh(seedA, ap);
    const auto a2 = render::generateArtifactMesh(seedA, ap);
    const auto b1 = render::generateArtifactMesh(seedB, ap);

    CHECK(!a1.vertices.empty());
    CHECK(!a1.indices.empty());
    CHECK(a1.indices.size() % 3u == 0u);

    // Determinism: exact counts should match for the same seed.
    CHECK(a1.vertices.size() == a2.vertices.size());
    CHECK(a1.indices.size() == a2.indices.size());

    // Variety: different seed should not produce an identical mesh.
    const std::size_t n = std::min<std::size_t>(a1.vertices.size(), 128);
    const std::size_t nb = std::min<std::size_t>(b1.vertices.size(), 128);

    double sumA = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      sumA += (double)a1.vertices[i].px + (double)a1.vertices[i].py + (double)a1.vertices[i].pz;
    }

    double sumB = 0.0;
    for (std::size_t i = 0; i < nb; ++i) {
      sumB += (double)b1.vertices[i].px + (double)b1.vertices[i].py + (double)b1.vertices[i].pz;
    }

    CHECK(std::abs(sumA - sumB) > 1e-3);
  }

  return failures;
}
