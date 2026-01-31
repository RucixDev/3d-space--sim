#include "stellar/sim/TrajectoryEvents.h"

#include "stellar/math/Vec3.h"
#include "stellar/sim/Units.h"

#include <algorithm>
#include <cstddef>

namespace stellar::sim {

namespace {

inline bool sameBodyId(const GravityBody& a, const GravityBody& b) {
  return a.kind == b.kind && a.id == b.id;
}

inline math::Vec3d lerp(const math::Vec3d& a, const math::Vec3d& b, double t) {
  return (a * (1.0 - t)) + (b * t);
}

void refineDominanceTransitions(const StarSystem& sys,
                               double startTimeDays,
                               double t0,
                               const math::Vec3d& p0,
                               const GravityBody& b0,
                               double t1,
                               const math::Vec3d& p1,
                               const GravityBody& b1,
                               int depth,
                               int maxDepth,
                               const GravityParams& gravityParams,
                               std::vector<DominantBodyTransition>& out) {
  if (sameBodyId(b0, b1)) return;
  if (!(t1 > t0)) {
    DominantBodyTransition e;
    e.tSec = t0;
    e.from = b0;
    e.to = b1;
    out.push_back(e);
    return;
  }

  if (depth >= maxDepth) {
    DominantBodyTransition e;
    e.tSec = 0.5 * (t0 + t1);
    e.from = b0;
    e.to = b1;
    out.push_back(e);
    return;
  }

  const double tMid = 0.5 * (t0 + t1);
  const math::Vec3d pMid = lerp(p0, p1, 0.5);

  const auto domMid = dominantGravityBody(sys, startTimeDays + (tMid / kSecondsPerDay), pMid, gravityParams);
  const GravityBody bMid = domMid.valid ? domMid.body : b0;

  if (sameBodyId(bMid, b0)) {
    refineDominanceTransitions(sys,
                              startTimeDays,
                              tMid,
                              pMid,
                              bMid,
                              t1,
                              p1,
                              b1,
                              depth + 1,
                              maxDepth,
                              gravityParams,
                              out);
    return;
  }

  if (sameBodyId(bMid, b1)) {
    refineDominanceTransitions(sys,
                              startTimeDays,
                              t0,
                              p0,
                              b0,
                              tMid,
                              pMid,
                              bMid,
                              depth + 1,
                              maxDepth,
                              gravityParams,
                              out);
    return;
  }

  // A third body becomes dominant inside this segment. Split and search both halves.
  refineDominanceTransitions(sys,
                            startTimeDays,
                            t0,
                            p0,
                            b0,
                            tMid,
                            pMid,
                            bMid,
                            depth + 1,
                            maxDepth,
                            gravityParams,
                            out);
  refineDominanceTransitions(sys,
                            startTimeDays,
                            tMid,
                            pMid,
                            bMid,
                            t1,
                            p1,
                            b1,
                            depth + 1,
                            maxDepth,
                            gravityParams,
                            out);
}

} // namespace

std::vector<DominantBodyTransition> detectDominantBodyTransitions(const StarSystem& sys,
                                                                 double startTimeDays,
                                                                 const std::vector<TrajectorySample>& samples,
                                                                 const GravityParams& gravityParams,
                                                                 DominantBodyTransitionParams params) {
  std::vector<DominantBodyTransition> out;
  if (samples.size() < 2) return out;

  const auto dom0 = dominantGravityBody(sys,
                                       startTimeDays + (samples.front().tSec / kSecondsPerDay),
                                       samples.front().posKm,
                                       gravityParams);
  if (!dom0.valid) return out;

  const GravityBody startBody = dom0.body;

  std::vector<DominantBodyTransition> raw;
  raw.reserve(16);

  GravityBody prevBody = startBody;
  double prevTSec = samples.front().tSec;
  math::Vec3d prevPos = samples.front().posKm;

  const int refineDepth = std::max(0, params.refineDepth);

  for (std::size_t i = 1; i < samples.size(); ++i) {
    const auto dom = dominantGravityBody(sys,
                                        startTimeDays + (samples[i].tSec / kSecondsPerDay),
                                        samples[i].posKm,
                                        gravityParams);
    if (!dom.valid) continue;

    const GravityBody curBody = dom.body;
    if (!sameBodyId(prevBody, curBody)) {
      refineDominanceTransitions(sys,
                                startTimeDays,
                                prevTSec,
                                prevPos,
                                prevBody,
                                samples[i].tSec,
                                samples[i].posKm,
                                curBody,
                                0,
                                refineDepth,
                                gravityParams,
                                raw);
    }

    prevBody = curBody;
    prevTSec = samples[i].tSec;
    prevPos = samples[i].posKm;
  }

  if (raw.empty()) return out;

  std::sort(raw.begin(), raw.end(), [](const DominantBodyTransition& a, const DominantBodyTransition& b) {
    return a.tSec < b.tSec;
  });

  // Compress into a consistent transition timeline.
  GravityBody current = startBody;
  for (const auto& e : raw) {
    if (sameBodyId(e.to, current)) continue;
    if (!out.empty() && params.minSeparationSec > 0.0 && (e.tSec - out.back().tSec) < params.minSeparationSec) {
      continue;
    }

    DominantBodyTransition t;
    t.tSec = e.tSec;
    t.from = current;
    t.to = e.to;
    out.push_back(t);
    current = e.to;
  }

  return out;
}

} // namespace stellar::sim
