#include "stellar/math/Bvh.h"

#include "test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

using namespace stellar;

namespace {

static math::Aabb3d randomBox(std::mt19937_64& rng,
                             double centerRange,
                             double minExtent,
                             double maxExtent) {
  std::uniform_real_distribution<double> cdist(-centerRange, centerRange);
  std::uniform_real_distribution<double> edist(minExtent, maxExtent);
  const math::Vec3d c{cdist(rng), cdist(rng), cdist(rng)};
  const math::Vec3d e{edist(rng), edist(rng), edist(rng)};
  return math::Aabb3d::fromCenterExtents(c, e);
}

static bool sameIdSet(std::vector<int> a, std::vector<int> b) {
  auto norm = [](std::vector<int>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  };
  norm(a);
  norm(b);
  return a == b;
}

static bool approx(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

} // namespace

int test_geometry_bvh() {
  int failures = 0;

  std::mt19937_64 rng(0xC0FFEEULL);

  // Build a medium-sized random scene.
  std::vector<math::BvhItem3d> items;
  items.reserve(256);
  for (int i = 0; i < 256; ++i) {
    math::BvhItem3d it{};
    it.bounds = randomBox(rng, /*centerRange=*/50.0, /*minExtent=*/0.2, /*maxExtent=*/4.0);
    it.id = i;
    items.push_back(it);
  }

  math::Bvh3d bvhMedian;
  bvhMedian.build(items, /*leafSize=*/6);

  math::Bvh3d bvhMorton;
  bvhMorton.build(items, /*leafSize=*/6, math::BvhBuildMode::Morton);

  CHECK(!bvhMedian.empty());
  CHECK(!bvhMorton.empty());
  CHECK(bvhMedian.itemCount() == items.size());
  CHECK(bvhMorton.itemCount() == items.size());
  CHECK(bvhMedian.nodeCount() > 0);
  CHECK(bvhMorton.nodeCount() > 0);
  CHECK(bvhMedian.buildMode() == math::BvhBuildMode::Median);
  CHECK(bvhMorton.buildMode() == math::BvhBuildMode::Morton);

  // ---- AABB overlap query equivalence (BVH vs brute force). ----
  for (int q = 0; q < 100; ++q) {
    const math::Aabb3d query = randomBox(rng, /*centerRange=*/60.0, /*minExtent=*/0.5, /*maxExtent=*/12.0);

    std::vector<int> hitsMedian;
    bvhMedian.queryAabb(query, [&](int id) { hitsMedian.push_back(id); });

    std::vector<int> hitsMorton;
    bvhMorton.queryAabb(query, [&](int id) { hitsMorton.push_back(id); });

    std::vector<int> brute;
    brute.reserve(items.size());
    for (const auto& it : items) {
      if (it.bounds.intersectsAabb(query)) brute.push_back(it.id);
    }

    CHECK(sameIdSet(hitsMedian, brute));
    CHECK(sameIdSet(hitsMorton, brute));
    CHECK(sameIdSet(hitsMedian, hitsMorton));
  }

  // ---- Frustum query equivalence. ----
  {
    const double pi = 3.14159265358979323846;
    const math::Mat4d proj = math::Mat4d::perspective(0.5 * pi, /*aspect=*/1.2, /*znear=*/0.5, /*zfar=*/200.0);
    const math::Mat4d view = math::Mat4d::lookAt({15.0, 8.0, 25.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0});
    const math::Frustumd fr = math::frustumFromViewProjection(view, proj);
    CHECK(fr.isFinite());

    std::vector<int> hitsMedian;
    bvhMedian.queryFrustum(fr, [&](int id) { hitsMedian.push_back(id); });

    std::vector<int> hitsMorton;
    bvhMorton.queryFrustum(fr, [&](int id) { hitsMorton.push_back(id); });

    std::vector<int> brute;
    brute.reserve(items.size());
    for (const auto& it : items) {
      if (fr.intersectsAabb(it.bounds)) brute.push_back(it.id);
    }

    CHECK(sameIdSet(hitsMedian, brute));
    CHECK(sameIdSet(hitsMorton, brute));
    CHECK(sameIdSet(hitsMedian, hitsMorton));
  }

  // ---- Raycast equivalence (nearest hit). ----
  {
    std::uniform_real_distribution<double> odist(-80.0, 80.0);
    std::uniform_real_distribution<double> ddist(-1.0, 1.0);

    for (int r = 0; r < 200; ++r) {
      const math::Vec3d o{odist(rng), odist(rng), odist(rng)};
      math::Vec3d d{ddist(rng), ddist(rng), ddist(rng)};
      if (d.lengthSq() < 1e-8) {
        d = {0.1, 0.2, -0.3};
      }

      double bruteT = std::numeric_limits<double>::infinity();
      std::vector<int> bestIds;
      for (const auto& it : items) {
        double t0 = 0.0, t1 = 0.0;
        if (!it.bounds.rayIntersectionT(o, d, t0, t1)) continue;
        if (t0 + 1e-10 < bruteT) {
          bruteT = t0;
          bestIds.clear();
          bestIds.push_back(it.id);
        } else if (approx(t0, bruteT, 1e-10)) {
          bestIds.push_back(it.id);
        }
      }

      // Median BVH
      {
        double t = 0.0;
        int id = -1;
        const bool hit = bvhMedian.raycast(o, d, &t, &id);
        const bool bruteHit = std::isfinite(bruteT);
        CHECK(hit == bruteHit);

        if (hit && bruteHit) {
          CHECK(approx(t, bruteT, 1e-8));
          const bool idOk = std::find(bestIds.begin(), bestIds.end(), id) != bestIds.end();
          CHECK(idOk);
        }
      }

      // Morton BVH
      {
        double t = 0.0;
        int id = -1;
        const bool hit = bvhMorton.raycast(o, d, &t, &id);
        const bool bruteHit = std::isfinite(bruteT);
        CHECK(hit == bruteHit);

        if (hit && bruteHit) {
          CHECK(approx(t, bruteT, 1e-8));
          const bool idOk = std::find(bestIds.begin(), bestIds.end(), id) != bestIds.end();
          CHECK(idOk);
        }
      }
    }
  }

  // ---- Degenerate centroid bounds: all centers equal should not explode. ----
  {
    std::vector<math::BvhItem3d> flat;
    flat.reserve(32);

    for (int i = 0; i < 32; ++i) {
      math::BvhItem3d it{};
      const double e = 0.1 + 0.02 * static_cast<double>(i);
      it.bounds = math::Aabb3d::fromCenterExtents({0.0, 0.0, 0.0}, {e, e * 1.1, e * 0.9});
      it.id = i;
      flat.push_back(it);
    }

    math::Bvh3d lbvh;
    lbvh.build(flat, /*leafSize=*/4, math::BvhBuildMode::Morton);
    CHECK(!lbvh.empty());

    const math::Aabb3d q = math::Aabb3d::fromCenterExtents({0.0, 0.0, 0.0}, {10.0, 10.0, 10.0});
    std::vector<int> hits;
    lbvh.queryAabb(q, [&](int id) { hits.push_back(id); });

    std::vector<int> brute;
    brute.reserve(flat.size());
    for (const auto& it : flat) {
      if (it.bounds.intersectsAabb(q)) brute.push_back(it.id);
    }

    CHECK(sameIdSet(hits, brute));
  }

  return failures;
}
