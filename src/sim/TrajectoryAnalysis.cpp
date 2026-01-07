#include "stellar/sim/TrajectoryAnalysis.h"

#include "stellar/math/Math.h"
#include "stellar/sim/Units.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::sim {

namespace {

struct BodyDesc {
  GravityBody::Kind kind{GravityBody::Kind::Planet};
  core::u64 id{0};
  std::string name{};
  double muKm3S2{0.0};
  double radiusKm{0.0};
  // For planets: index into StarSystem::planets.
  int planetIndex{-1};
};

static math::Vec3d lerp(const math::Vec3d& a, const math::Vec3d& b, double t) {
  return a + (b - a) * t;
}

static GravityBody evalBodyAtTime(const StarSystem& sys, const BodyDesc& d, double timeDays) {
  GravityBody b;
  b.kind = d.kind;
  b.id = d.id;
  b.name = d.name;
  b.muKm3S2 = d.muKm3S2;
  b.radiusKm = d.radiusKm;

  if (d.kind == GravityBody::Kind::Star) {
    b.posKm = {0, 0, 0};
    b.velKmS = {0, 0, 0};
    return b;
  }

  const int i = d.planetIndex;
  if (i >= 0 && i < (int)sys.planets.size()) {
    const auto& p = sys.planets[(size_t)i];
    b.posKm = planetPosKm(p, timeDays);
    b.velKmS = planetVelKmS(p, timeDays);
  } else {
    b.posKm = {0, 0, 0};
    b.velKmS = {0, 0, 0};
  }
  return b;
}

static void buildBodyDescs(const StarSystem& sys,
                           const GravityParams& params,
                           std::vector<BodyDesc>& out) {
  out.clear();

  if (params.includeStar) {
    BodyDesc d;
    d.kind = GravityBody::Kind::Star;
    d.id = sys.stub.id;
    d.name = sys.star.cls == StarClass::Count ? std::string{} : sys.stub.name;
    // Prefer StarSystem stub name (if set), otherwise leave empty.
    if (d.name.empty()) d.name = "Star";
    d.muKm3S2 = muStarKm3S2(sys.star);
    d.radiusKm = radiusStarKm(sys.star);
    d.planetIndex = -1;
    out.push_back(std::move(d));
  }

  if (params.includePlanets) {
    for (int i = 0; i < (int)sys.planets.size(); ++i) {
      const auto& p = sys.planets[(size_t)i];
      BodyDesc d;
      d.kind = GravityBody::Kind::Planet;
      d.id = (core::u64)i;
      d.name = p.name;
      d.muKm3S2 = muPlanetKm3S2(p);
      d.radiusKm = radiusPlanetKm(p);
      d.planetIndex = i;
      out.push_back(std::move(d));
    }
  }
}

static bool solveFirstIntersectionU(const math::Vec3d& rel0,
                                    const math::Vec3d& rel1,
                                    double radiusKm,
                                    double& outU) {
  outU = 0.0;
  if (!(radiusKm > 0.0)) return false;

  const double r0 = rel0.length();
  const double r1 = rel1.length();

  // We only care about the first outside->inside transition.
  if (!(r0 > radiusKm && r1 <= radiusKm)) return false;

  const math::Vec3d d = rel1 - rel0;
  const double a = d.lengthSq();
  if (!(a > 1e-18)) return false;

  const double b = 2.0 * math::dot(rel0, d);
  const double c = rel0.lengthSq() - radiusKm * radiusKm;
  const double disc = b * b - 4.0 * a * c;
  if (disc < 0.0) return false;

  const double s = std::sqrt(std::max(0.0, disc));
  const double inv2a = 1.0 / (2.0 * a);
  const double u1 = (-b - s) * inv2a;
  const double u2 = (-b + s) * inv2a;

  // Smallest root in [0,1].
  double u = std::numeric_limits<double>::infinity();
  if (u1 >= 0.0 && u1 <= 1.0) u = std::min(u, u1);
  if (u2 >= 0.0 && u2 <= 1.0) u = std::min(u, u2);
  if (!std::isfinite(u)) return false;

  outU = u;
  return true;
}

} // namespace

TrajectoryAnalysisResult analyzeTrajectory(const StarSystem& sys,
                                          double startTimeDays,
                                          const std::vector<TrajectorySample>& samples,
                                          const GravityParams& params) {
  TrajectoryAnalysisResult out{};
  if (samples.empty()) return out;

  std::vector<BodyDesc> bodies;
  buildBodyDescs(sys, params, bodies);
  if (bodies.empty()) return out;

  out.closestApproachByBody.resize(bodies.size());
  for (auto& a : out.closestApproachByBody) {
    a.valid = false;
    a.distanceKm = std::numeric_limits<double>::infinity();
    a.altitudeKm = std::numeric_limits<double>::infinity();
  }

  out.closestOverall.valid = false;
  out.closestOverall.distanceKm = std::numeric_limits<double>::infinity();
  out.closestOverall.altitudeKm = std::numeric_limits<double>::infinity();

  // Handle the degenerate case of a single sample: treat it as both the closest
  // approach and potential impact state.
  if (samples.size() == 1) {
    const auto& s = samples[0];
    const double timeDays = startTimeDays + s.tSec / kSecondsPerDay;
    for (size_t bi = 0; bi < bodies.size(); ++bi) {
      const auto b = evalBodyAtTime(sys, bodies[bi], timeDays);
      const auto rel = s.posKm - b.posKm;
      const double dist = rel.length();
      const double alt = dist - b.radiusKm;

      TrajectoryApproach a;
      a.valid = true;
      a.body = b;
      a.tSec = s.tSec;
      a.distanceKm = dist;
      a.altitudeKm = alt;
      a.relPosKm = rel;
      a.relVelKmS = s.velKmS - b.velKmS;
      out.closestApproachByBody[bi] = a;

      if (dist < out.closestOverall.distanceKm) {
        out.closestOverall = a;
      }

      if (!out.firstImpact.valid && alt <= 0.0) {
        out.firstImpact.valid = true;
        out.firstImpact.body = b;
        out.firstImpact.tSec = s.tSec;
        out.firstImpact.posKm = s.posKm;
        out.firstImpact.velKmS = s.velKmS;
      }
    }
    return out;
  }

  // Segment-based scan (more accurate than sampling-only).
  for (size_t si = 1; si < samples.size(); ++si) {
    const auto& s0 = samples[si - 1];
    const auto& s1 = samples[si];
    const double dt = s1.tSec - s0.tSec;
    if (!(dt > 0.0)) continue;

    const double tDays0 = startTimeDays + s0.tSec / kSecondsPerDay;
    const double tDays1 = startTimeDays + s1.tSec / kSecondsPerDay;

    for (size_t bi = 0; bi < bodies.size(); ++bi) {
      const auto b0 = evalBodyAtTime(sys, bodies[bi], tDays0);
      const auto b1 = evalBodyAtTime(sys, bodies[bi], tDays1);

      const math::Vec3d rel0 = s0.posKm - b0.posKm;
      const math::Vec3d rel1 = s1.posKm - b1.posKm;
      const math::Vec3d d = rel1 - rel0;
      const double dd = d.lengthSq();

      // Closest approach along the segment in body-centric space.
      double uMin = 0.0;
      if (dd > 1e-18) {
        uMin = std::clamp(-math::dot(rel0, d) / dd, 0.0, 1.0);
      }

      const math::Vec3d relMin = rel0 + d * uMin;
      const double distMin = std::sqrt(std::max(0.0, relMin.lengthSq()));
      const double altMin = distMin - b0.radiusKm;
      const double tMinSec = s0.tSec + uMin * dt;
      const double tMinDays = startTimeDays + tMinSec / kSecondsPerDay;
      const auto bAtMin = evalBodyAtTime(sys, bodies[bi], tMinDays);

      if (distMin < out.closestApproachByBody[bi].distanceKm) {
        TrajectoryApproach a;
        a.valid = true;
        a.body = bAtMin;
        a.tSec = tMinSec;
        a.distanceKm = distMin;
        a.altitudeKm = distMin - bAtMin.radiusKm;
        a.relPosKm = relMin;

        const math::Vec3d shipVelMin = lerp(s0.velKmS, s1.velKmS, uMin);
        const math::Vec3d bodyVelMin = lerp(b0.velKmS, b1.velKmS, uMin);
        a.relVelKmS = shipVelMin - bodyVelMin;

        out.closestApproachByBody[bi] = a;
      }

      if (distMin < out.closestOverall.distanceKm) {
        out.closestOverall = out.closestApproachByBody[bi];
      }

      // First impact detection (surface intersection).
      if (!out.firstImpact.valid) {
        double uHit = 0.0;
        if (solveFirstIntersectionU(rel0, rel1, b0.radiusKm, uHit)) {
          const double tHitSec = s0.tSec + uHit * dt;
          const math::Vec3d pHit = lerp(s0.posKm, s1.posKm, uHit);
          const math::Vec3d vHit = lerp(s0.velKmS, s1.velKmS, uHit);
          const double tHitDays = startTimeDays + tHitSec / kSecondsPerDay;
          const auto bAtHit = evalBodyAtTime(sys, bodies[bi], tHitDays);

          out.firstImpact.valid = true;
          out.firstImpact.body = bAtHit;
          out.firstImpact.tSec = tHitSec;
          out.firstImpact.posKm = pHit;
          out.firstImpact.velKmS = vHit;
        }
      }
    }
  }

  return out;
}

} // namespace stellar::sim
