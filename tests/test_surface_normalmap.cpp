#include "stellar/render/ProceduralPlanet.h"

#include "test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>

using namespace stellar;

static double decodeNormal01To11(std::uint8_t c) {
  return (double)c / 255.0 * 2.0 - 1.0;
}

int test_surface_normalmap() {
  using namespace stellar::render;

  int failures = 0;

  const core::u64 seed = 0x123456789ABCDEF0ull;
  const int w = 128;

  const auto a = generateSurfaceNormalMap(SurfaceKind::Rocky, seed, w);
  const auto b = generateSurfaceNormalMap(SurfaceKind::Rocky, seed, w);

  CHECK(a.w == w);
  CHECK(a.h == w / 2);
  CHECK(a.w == b.w);
  CHECK(a.h == b.h);
  CHECK(a.rgba == b.rgba);

  CHECK(a.rgba.size() == (std::size_t)a.w * (std::size_t)a.h * 4);

  double minLen = 1e9;
  double maxLen = 0.0;
  double meanLen = 0.0;

  double minX = 1e9, maxX = -1e9;
  double minY = 1e9, maxY = -1e9;
  double minZ = 1e9, maxZ = -1e9;

  const std::size_t pxCount = (std::size_t)a.w * (std::size_t)a.h;
  for (std::size_t i = 0; i < pxCount; ++i) {
    const std::size_t idx = i * 4;
    const double nx = decodeNormal01To11(a.rgba[idx + 0]);
    const double ny = decodeNormal01To11(a.rgba[idx + 1]);
    const double nz = decodeNormal01To11(a.rgba[idx + 2]);

    const double len = std::sqrt(nx * nx + ny * ny + nz * nz);
    meanLen += len;
    minLen = std::min(minLen, len);
    maxLen = std::max(maxLen, len);

    minX = std::min(minX, nx);
    maxX = std::max(maxX, nx);
    minY = std::min(minY, ny);
    maxY = std::max(maxY, ny);
    minZ = std::min(minZ, nz);
    maxZ = std::max(maxZ, nz);
  }
  if (pxCount > 0) meanLen /= (double)pxCount;

  // Encoded normals should be close to unit length, though quantization introduces error.
  CHECK(meanLen > 0.95);
  CHECK(meanLen < 1.05);
  CHECK(minLen > 0.80);
  CHECK(maxLen < 1.20);

  // Should not be entirely flat.
  CHECK((maxX - minX) > 0.02);
  CHECK((maxY - minY) > 0.02);

  // Z should generally face outward.
  CHECK(maxZ > 0.65);

  if (failures == 0) {
    std::cout << "[test_surface_normalmap] PASS\n";
  }
  return failures;
}
