#include "stellar/render/ProceduralRings.h"

#include "stellar/core/Hash.h"

#include "test_harness.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>

using namespace stellar;

static core::u64 imageHash(const render::RingImage& img) {
  // Include dimensions to avoid accidental collisions when testing multiple res.
  core::u64 h = core::hashCombine(core::fnv1a64("ring_img"), (core::u64)img.w);
  h = core::hashCombine(h, (core::u64)img.h);
  if (!img.rgba.empty()) {
    h = core::hashCombine(h, core::hashBytes(img.rgba.data(), img.rgba.size()));
  }
  return h;
}

int test_rings() {
  int failures = 0;

  const core::u64 seedA = 0x123456789abcdef0ull;
  const core::u64 seedB = 0x0fedcba987654321ull;

  const int w = 256;
  const int h = 80;

  const render::RingImage a0 = render::generateRingTexture(seedA, w, h);
  const render::RingImage a1 = render::generateRingTexture(seedA, w, h);
  const render::RingImage b0 = render::generateRingTexture(seedB, w, h);

  CHECK(a0.w == w && a0.h == h);
  CHECK((int)a0.rgba.size() == w * h * 4);

  // Deterministic for the same seed/res.
  CHECK(imageHash(a0) == imageHash(a1));

  // Different seeds should generally differ.
  CHECK(imageHash(a0) != imageHash(b0));

  // Seam continuity: u=0 and u=1 columns should match (or be extremely close).
  int seamMax = 0;
  for (int y = 0; y < h; ++y) {
    const std::size_t i0 = (std::size_t)(y * w + 0) * 4;
    const std::size_t i1 = (std::size_t)(y * w + (w - 1)) * 4;
    for (int c = 0; c < 4; ++c) {
      const int d = std::abs((int)a0.rgba[i0 + (std::size_t)c] - (int)a0.rgba[i1 + (std::size_t)c]);
      seamMax = std::max(seamMax, d);
    }
  }
  CHECK(seamMax <= 2);

  // Basic signal sanity: alpha should have both sparse and dense regions.
  std::uint8_t minA = 255;
  std::uint8_t maxA = 0;
  int countSparse = 0;
  int countDense = 0;
  const int nPix = w * h;

  for (int i = 0; i < nPix; ++i) {
    const std::uint8_t A = a0.rgba[(std::size_t)i * 4 + 3];
    minA = std::min(minA, A);
    maxA = std::max(maxA, A);
    if (A < 30) ++countSparse;
    if (A > 170) ++countDense;
  }

  CHECK(minA < 30);
  CHECK(maxA > 170);
  CHECK(countSparse > (int)(0.02 * (double)nPix));
  CHECK(countDense > (int)(0.05 * (double)nPix));

  // Not monochrome: some RGB variation should exist.
  std::uint8_t minR = 255, maxR = 0;
  std::uint8_t minG = 255, maxG = 0;
  std::uint8_t minB = 255, maxB = 0;
  for (int i = 0; i < nPix; ++i) {
    const std::uint8_t R = a0.rgba[(std::size_t)i * 4 + 0];
    const std::uint8_t G = a0.rgba[(std::size_t)i * 4 + 1];
    const std::uint8_t B = a0.rgba[(std::size_t)i * 4 + 2];
    minR = std::min(minR, R); maxR = std::max(maxR, R);
    minG = std::min(minG, G); maxG = std::max(maxG, G);
    minB = std::min(minB, B); maxB = std::max(maxB, B);
  }

  CHECK((int)(maxR - minR) > 8);
  CHECK((int)(maxG - minG) > 8);
  CHECK((int)(maxB - minB) > 8);

  if (failures) {
    std::cerr << "[stellar_tests] test_rings failures=" << failures << "\n";
  }
  return failures;
}