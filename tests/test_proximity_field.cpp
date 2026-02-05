#include "stellar/sim/ProximityField.h"

#include "test_harness.h"

#include <cmath>
#include <algorithm>
#include <vector>

using namespace stellar;

namespace {

static bool approx(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

} // namespace

int test_proximity_field() {
  int failures = 0;

  // ---- Raycast hit / ignoreId. ----
  {
    sim::ProximityFieldKm f;
    std::vector<sim::SphereObstacleKm> o;
    o.push_back(sim::SphereObstacleKm{{0.0, 0.0, 5.0}, 1.0, 1.0}); // id 0
    o.push_back(sim::SphereObstacleKm{{0.0, 0.0, 9.0}, 1.0, 1.0}); // id 1
    f.build(o, /*leafSize=*/4);

    const auto hit = f.raycastFirstHit({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 10.0);
    CHECK(hit.hit);
    CHECK(hit.id == 0);
    CHECK(approx(hit.tKm, 4.0, 1e-9));

    const auto hitIgnore = f.raycastFirstHit({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 10.0, /*padKm=*/0.0, /*ignoreId=*/0);
    CHECK(hitIgnore.hit);
    CHECK(hitIgnore.id == 1);
    CHECK(approx(hitIgnore.tKm, 8.0, 1e-9));
  }

  // ---- Linear time-to-impact along velocity. ----
  {
    sim::ProximityFieldKm f;
    std::vector<sim::SphereObstacleKm> o;
    o.push_back(sim::SphereObstacleKm{{0.0, 0.0, 5.0}, 1.0, 1.0});
    f.build(o, /*leafSize=*/2);

    const auto hit = f.predictLinearImpact({0.0, 0.0, 0.0}, {0.0, 0.0, 2.0}, /*maxTimeSec=*/10.0);
    CHECK(hit.hit);
    CHECK(approx(hit.tKm, 4.0, 1e-9));
    CHECK(approx(hit.tSec, 2.0, 1e-9));
  }


  // ---- Frustum query broadphase equivalence. ----
  {
    sim::ProximityFieldKm f;
    std::vector<sim::SphereObstacleKm> o;
    o.push_back(sim::SphereObstacleKm{{0.0, 0.0, -6.0}, 1.0, 1.0});  // id 0 (in front)
    o.push_back(sim::SphereObstacleKm{{25.0, 0.0, -6.0}, 1.0, 1.0}); // id 1 (far right)
    o.push_back(sim::SphereObstacleKm{{0.0, 0.0, 6.0}, 1.0, 1.0});   // id 2 (behind camera)
    f.build(o, /*leafSize=*/2);

    const double pi = 3.14159265358979323846;
    const math::Mat4d proj = math::Mat4d::perspective(0.5 * pi, /*aspect=*/1.0, /*znear=*/0.5, /*zfar=*/50.0);
    const math::Mat4d view = math::Mat4d::lookAt({0.0, 0.0, 0.0}, {0.0, 0.0, -1.0}, {0.0, 1.0, 0.0});
    const math::Frustumd fr = math::frustumFromViewProjection(view, proj);
    CHECK(fr.isFinite());

    std::vector<int> hits;
    f.queryFrustum(fr, [&](int id) { hits.push_back(id); });

    std::vector<int> brute;
    brute.reserve(o.size());
    for (int i = 0; i < (int)o.size(); ++i) {
      const auto& s = o[(std::size_t)i];
      const math::Aabb3d aabb = math::Aabb3d::fromCenterExtents(
        s.centerKm, {s.radiusKm, s.radiusKm, s.radiusKm});
      if (fr.intersectsAabb(aabb)) brute.push_back(i);
    }

    std::sort(hits.begin(), hits.end());
    hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
    std::sort(brute.begin(), brute.end());
    brute.erase(std::unique(brute.begin(), brute.end()), brute.end());
    CHECK(hits == brute);
  }
  // ---- Avoidance steering: obstacle centered on the desired path. ----
  {
    sim::ProximityFieldKm f;
    std::vector<sim::SphereObstacleKm> o;
    o.push_back(sim::SphereObstacleKm{{0.0, 0.0, 8.0}, 1.0, 1.0});
    f.build(o, /*leafSize=*/2);

    sim::AvoidanceParamsKm p;
    p.enabled = true;
    p.shipRadiusKm = 0.0;
    p.padKm = 0.0;
    p.lookaheadTimeSec = 10.0;
    p.lookaheadBaseKm = 0.0;
    p.minSpeedForLookaheadKmS = 1.0;
    p.nearMissExtraKm = 2.0;
    p.strength = 1.5;
    p.maxSteerDeg = 45.0;

    const auto ar = sim::steerAvoidObstacles(f,
                                           /*posKm=*/{0.0, 0.0, 0.0},
                                           /*velKmS=*/{0.0, 0.0, 1.0},
                                           /*desiredDirUnitWorld=*/{0.0, 0.0, 1.0},
                                           /*desiredSpeedKmS=*/1.0,
                                           p,
                                           /*ignoreId=*/-1);

    CHECK(ar.lookaheadKm > 0.0);
    CHECK(ar.threatId == 0);
    CHECK(ar.steering);
    CHECK(ar.safeDirUnit.lengthSq() > 0.99);
    CHECK(math::dot(ar.desiredDirUnit, ar.safeDirUnit) < 0.9999);

    const auto arIgnore = sim::steerAvoidObstacles(f,
                                                 /*posKm=*/{0.0, 0.0, 0.0},
                                                 /*velKmS=*/{0.0, 0.0, 1.0},
                                                 /*desiredDirUnitWorld=*/{0.0, 0.0, 1.0},
                                                 /*desiredSpeedKmS=*/1.0,
                                                 p,
                                                 /*ignoreId=*/0);
    CHECK(!arIgnore.steering);
    CHECK(math::dot(arIgnore.desiredDirUnit, arIgnore.safeDirUnit) > 0.999999);
  }



  // ---- Tangent bypass waypoint: single obstacle blocks the direct segment. ----
  {
    sim::ProximityFieldKm f;
    std::vector<sim::SphereObstacleKm> o;
    o.push_back(sim::SphereObstacleKm{{0.0, 0.0, 5.0}, 1.0, 1.0}); // id 0
    f.build(o, /*leafSize=*/2);

    sim::TangentBypassParamsKm bp;
    bp.enabled = true;
    bp.shipRadiusKm = 0.0;
    bp.padKm = 0.0;
    bp.nearMissExtraKm = 0.5;
    bp.aheadClearanceMult = 1.0;
    bp.aheadExtraKm = 0.0;

    const math::Vec3d start{0.0, 0.0, 0.0};
    const math::Vec3d goal{0.0, 0.0, 10.0};

    // Sanity: direct segment should hit when inflated by the near-miss bubble.
    const double queryPad = bp.shipRadiusKm + bp.padKm + bp.nearMissExtraKm;
    const auto direct = f.raycastFirstHit(start, {0.0, 0.0, 1.0}, 10.0, queryPad);
    CHECK(direct.hit);
    CHECK(direct.id == 0);

    const auto br = sim::planTangentBypassWaypoint(f, start, goal, bp);
    CHECK(br.used);
    CHECK(br.obstacleId == 0);
    CHECK(br.sideSign == -1 || br.sideSign == 1);

    // The two legs should be clear with the same inflation used for planning.
    const auto toW = br.waypointKm - start;
    const double d1 = toW.length();
    CHECK(d1 > 0.1);
    const auto hit1 = f.raycastFirstHit(start, toW.normalized(), d1, queryPad);
    CHECK(!hit1.hit || hit1.tKm >= d1 - 1e-6);

    const auto toG = goal - br.waypointKm;
    const double d2 = toG.length();
    CHECK(d2 > 0.1);
    const auto hit2 = f.raycastFirstHit(br.waypointKm, toG.normalized(), d2, queryPad);
    CHECK(!hit2.hit || hit2.tKm >= d2 - 1e-6);

    // Preferred side should be respected when valid.
    const auto brPos = sim::planTangentBypassWaypoint(f, start, goal, bp, /*ignoreId=*/-1, /*preferredSideSign=*/1);
    CHECK(brPos.used);
    CHECK(brPos.sideSign == 1);

    const auto brNeg = sim::planTangentBypassWaypoint(f, start, goal, bp, /*ignoreId=*/-1, /*preferredSideSign=*/-1);
    CHECK(brNeg.used);
    CHECK(brNeg.sideSign == -1);
  }

  return failures;
}
