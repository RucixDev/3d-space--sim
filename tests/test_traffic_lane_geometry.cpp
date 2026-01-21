#include "test_harness.h"

#include "stellar/sim/TrafficLanes.h"
#include "stellar/sim/Universe.h"
#include "stellar/sim/WorldIds.h"

#include <cmath>
#include <algorithm>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

static bool nearlyVec(const stellar::math::Vec3d& a, const stellar::math::Vec3d& b, double eps = 1e-6) {
  return nearly(a.x, b.x, eps) && nearly(a.y, b.y, eps) && nearly(a.z, b.z, eps);
}

int test_traffic_lane_geometry() {
  int failures = 0;

  using namespace stellar;
  using namespace stellar::sim;

  Universe u(1337ull);

  // Find a system with at least 2 stations.
  const auto stubs = u.queryNearby({0, 0, 0}, 80.0, 64);
  const SystemStub* chosen = nullptr;
  for (const auto& s : stubs) {
    if (s.stationCount >= 2) { chosen = &s; break; }
  }
  CHECK(chosen != nullptr);
  if (!chosen) return failures;

  const auto& sys = u.getSystem(chosen->id, chosen);
  CHECK(sys.stations.size() >= 2);
  if (sys.stations.size() < 2) return failures;

  const StationId A = sys.stations[0].id;
  const StationId B = sys.stations[1].id;
  CHECK(A != B);

  TrafficLaneParams p{};
  p.bundleByStationPair = true;
  p.dualCarriageway = true;
  // Force a non-trivial arc so mid-flight offsets are measurable.
  p.arcMinKm = 6000.0;
  p.arcMaxKm = 90000.0;
  p.arcMaxFracOfDistance = 0.35;

  TrafficConvoy c1{};
  c1.id = makeDeterministicWorldId(123, 456);
  c1.systemId = sys.stub.id;
  c1.fromStation = A;
  c1.toStation = B;
  c1.factionId = sys.stations[0].factionId;
  c1.commodity = econ::CommodityId::Food;
  c1.units = 10.0;
  c1.departDay = 100.25;
  c1.arriveDay = 100.45;

  TrafficConvoy c2 = c1;
  c2.id = makeDeterministicWorldId(123, 457);
  c2.commodity = econ::CommodityId::Ore;

  TrafficConvoy c3 = c1;
  c3.id = makeDeterministicWorldId(123, 458);
  c3.fromStation = B;
  c3.toStation = A;

  // --- Bundled corridors: same station pair => same laneKey/arc/side ---------
  const auto g1 = computeTrafficLaneGeometry(c1, sys, p);
  const auto g2 = computeTrafficLaneGeometry(c2, sys, p);

  CHECK(g1.laneKey != 0);
  CHECK(g2.laneKey != 0);
  CHECK(g1.laneKey == g2.laneKey);

  CHECK(nearly(g1.arcMagKm, g2.arcMagKm, 1e-9));
  CHECK(nearlyVec(g1.side, g2.side, 1e-6));
  CHECK(g1.directionSign == g2.directionSign);

  // --- Dual carriageway: reversed direction => same corridor key, sign flip ---
  const auto g3 = computeTrafficLaneGeometry(c3, sys, p);
  CHECK(g3.laneKey == g1.laneKey);
  CHECK(nearly(g3.arcMagKm, g1.arcMagKm, 1e-9));
  CHECK(g1.directionSign == -g3.directionSign);

  // --- Mid-flight arc sanity: offset is perpendicular to chord and matches arcMag
  const double mid = 0.5 * (c1.departDay + c1.arriveDay);
  const auto stMid = evaluateTrafficConvoy(c1, sys, mid, p);

  const math::Vec3d chord = g1.p1Km - g1.p0Km;
  const math::Vec3d base = g1.p0Km + chord * 0.5;
  const math::Vec3d off = stMid.posKm - base;

  // At phase=0.5, sin^2(pi*phase)=1 so dot(off, side) ~= arcMag (+ jitter).
  const double expectedArc = std::max(0.0, g1.arcMagKm + g1.arcJitterKm);
  const double alongSide = math::dot(off, g1.side);

  CHECK(std::abs(alongSide - expectedArc) < 1e-3);
  CHECK(std::abs(math::dot(off, g1.dir)) < 1e-3);

  // --- Legacy mode: unbundled => laneKey differs per convoy id ---------------
  TrafficLaneParams legacy = p;
  legacy.bundleByStationPair = false;

  const auto lg1 = computeTrafficLaneGeometry(c1, sys, legacy);
  const auto lg2 = computeTrafficLaneGeometry(c2, sys, legacy);

  CHECK(lg1.laneKey != 0);
  CHECK(lg2.laneKey != 0);
  CHECK(lg1.laneKey != lg2.laneKey);

  return failures;
}
