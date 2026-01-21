#include "stellar/sim/TimeTrial.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static math::Vec3d safeNormalize(const math::Vec3d& v, const math::Vec3d& fallback) {
  const double ls = v.lengthSq();
  if (ls < 1e-18) return fallback;
  return v / std::sqrt(ls);
}

static math::Vec3d anyPerpendicularUnit(const math::Vec3d& n) {
  // Pick a reference vector that is unlikely to be parallel to n.
  const math::Vec3d ref = (std::abs(n.x) < 0.6) ? math::Vec3d{1,0,0} : math::Vec3d{0,1,0};
  math::Vec3d u = math::cross(n, ref);
  if (u.lengthSq() < 1e-18) {
    u = math::cross(n, {0,0,1});
  }
  return safeNormalize(u, {1,0,0});
}

core::u64 hashPosKmQuantized(const math::Vec3d& pKm, double quantumKm) {
  const double q = (quantumKm <= 0.0) ? 1.0 : quantumKm;
  const auto qi = [&](double v) -> core::i64 {
    return static_cast<core::i64>(std::llround(v / q));
  };
  core::u64 h = core::fnv1a64("posKm");
  h = core::hashCombine(h, static_cast<core::u64>(qi(pKm.x)));
  h = core::hashCombine(h, static_cast<core::u64>(qi(pKm.y)));
  h = core::hashCombine(h, static_cast<core::u64>(qi(pKm.z)));
  return h;
}

TimeTrialCourse generateTimeTrialCourseStationSlalomKm(const math::Vec3d& anchorPosKm,
                                                       const math::Quatd& anchorOrient,
                                                       core::u64 seed,
                                                       const TimeTrialCourseParams& p) {
  TimeTrialCourse out;
  out.seed = seed;
  out.name = "Station Slalom";

  const int n = std::max(2, p.gateCount);
  const int loops = std::max(1, p.loops);

  // Stable key: seed + params + quantized anchor.
  core::u64 key = core::hashCombine(core::fnv1a64("TimeTrialCourse"), seed);
  key = core::hashCombine(key, static_cast<core::u64>(n));
  key = core::hashCombine(key, static_cast<core::u64>(loops));
  key = core::hashCombine(key, hashPosKmQuantized(anchorPosKm, 10.0));
  out.key = key;

  core::SplitMix64 rng(seed ^ 0xA1B2C3D4u);

  const double baseR = std::max(1000.0, p.courseRadiusKm);
  const double height = std::max(0.0, p.courseHeightKm);
  const double jitter = std::max(0.0, p.jitterKm);
  const double gateR = std::max(250.0, p.gateRadiusKm);

  const math::Vec3d right = anchorOrient.rotate({1,0,0});
  const math::Vec3d up    = anchorOrient.rotate({0,1,0});
  const math::Vec3d fwd   = anchorOrient.rotate({0,0,1});

  std::vector<math::Vec3d> centers;
  centers.resize(n);

  // Use a smooth-ish random walk for radius/height to avoid pure jitter noise.
  double rWalk = 0.0;
  double hWalk = 0.0;

  for (int i = 0; i < n; ++i) {
    const double t = (n > 1) ? (double)i / (double)(n - 1) : 0.0;
    const double ang = (2.0 * math::kPi) * (double)loops * t + rng.range(-0.15, 0.15);

    // Random walk (clamped) so neighboring gates have some continuity.
    rWalk = std::clamp(rWalk + rng.range(-1.0, 1.0) * 0.25, -1.0, 1.0);
    hWalk = std::clamp(hWalk + rng.range(-1.0, 1.0) * 0.35, -1.0, 1.0);

    const double r = baseR + rWalk * jitter;
    const double h = hWalk * height + std::sin(ang * 0.5) * (0.15 * height);

    const double ca = std::cos(ang);
    const double sa = std::sin(ang);

    centers[i] = anchorPosKm + right * (ca * r) + fwd * (sa * r) + up * h;
  }

  out.gates.resize(n);
  for (int i = 0; i < n; ++i) {
    const math::Vec3d a = centers[i];

    math::Vec3d dir{};
    if (i + 1 < n) {
      dir = centers[i + 1] - centers[i];
    } else if (p.closedLoop && n > 1) {
      dir = centers[0] - centers[i];
    } else if (i > 0) {
      dir = centers[i] - centers[i - 1];
    } else {
      dir = fwd;
    }

    out.gates[i].posKm = a;
    out.gates[i].normal = safeNormalize(dir, fwd);
    out.gates[i].radiusKm = gateR;
  }

  // If the last gate isn't part of a loop, give it a reasonable normal so the
  // player has a direction to fly through.
  if (!p.closedLoop && n >= 2) {
    out.gates[n - 1].normal = safeNormalize(out.gates[n - 1].normal, out.gates[n - 2].normal);
  }

  return out;
}

bool timeTrialGatePassed(const TimeTrialGate& g,
                         const math::Vec3d& prevPosKm,
                         const math::Vec3d& posKm,
                         const math::Vec3d& velKmS,
                         double* outTHit) {
  const math::Vec3d n = safeNormalize(g.normal, {0,0,1});

  const math::Vec3d a = prevPosKm - g.posKm;
  const math::Vec3d b = posKm - g.posKm;

  const double d0 = math::dot(a, n);
  const double d1 = math::dot(b, n);

  // Must cross from negative to positive side (forward direction).
  if (!(d0 < 0.0 && d1 >= 0.0)) return false;

  // Segment intersection with plane: prev + t*(pos-prev)
  const double denom = (d0 - d1);
  double t = 0.0;
  if (std::abs(denom) > 1e-12) {
    t = d0 / denom; // in [0,1] for segment
  }
  t = std::clamp(t, 0.0, 1.0);

  if (outTHit) *outTHit = t;
  const math::Vec3d p = prevPosKm + (posKm - prevPosKm) * t;

  const double r2 = g.radiusKm * g.radiusKm;
  if ((p - g.posKm).lengthSq() > r2) return false;

  // Extra robustness: if we have a meaningful velocity, require it to be at
  // least loosely aligned with the gate forward direction.
  if (velKmS.lengthSq() > 1e-12) {
    const double fwd = math::dot(velKmS, n);
    if (fwd < -1e-4) return false;
  }

  return true;
}

} // namespace stellar::sim
