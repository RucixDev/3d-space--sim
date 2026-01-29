#include "stellar/sim/FieldRepair.h"

#include <cmath>

namespace stellar::sim {

static inline double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

FieldRepairQuote quoteFieldRepairToCap(const FieldRepairInventory& inv,
                                       double hullCurrent,
                                       double hullMax,
                                       const FieldRepairParams& p) {
  FieldRepairQuote q{};
  q.hullCurrent = hullCurrent;
  q.hullMax = hullMax;

  const double eps = std::max(0.0, p.eps);

  if (!(hullMax > eps)) {
    q.ok = false;
    q.reason = "invalid hullMax";
    return q;
  }

  const double capFrac = clamp01(p.maxRepairFrac);
  q.capHull = hullMax * capFrac;
  q.limitedByCap = (q.capHull + eps < hullMax);

  q.hullMissingToCap = std::max(0.0, q.capHull - hullCurrent);

  if (q.hullMissingToCap <= eps) {
    q.ok = false;
    q.reason = "hull at/above field repair cap";
    return q;
  }

  const double mPer = std::max(0.0, p.metalsPerHull);
  const double wPer = std::max(0.0, p.machineryPerHull);

  // Compute max repairable hull from each material. If recipe is zero, treat as unlimited.
  double maxByMetals = q.hullMissingToCap;
  double maxByMach = q.hullMissingToCap;

  if (mPer > eps) maxByMetals = std::max(0.0, inv.metals) / mPer;
  if (wPer > eps) maxByMach = std::max(0.0, inv.machinery) / wPer;

  const double maxByStock = std::max(0.0, std::min(maxByMetals, maxByMach));

  q.hullToRepair = std::max(0.0, std::min(q.hullMissingToCap, maxByStock));
  q.limitedByStock = (q.hullToRepair + eps < q.hullMissingToCap);

  if (q.hullToRepair <= eps) {
    q.ok = false;
    q.reason = "insufficient materials";
    return q;
  }

  q.metalsNeeded = q.hullToRepair * mPer;
  q.machineryNeeded = q.hullToRepair * wPer;

  const double rate = std::max(p.hullPerRealSec, eps);
  q.timeSec = q.hullToRepair / rate;

  q.ok = true;
  q.reason = nullptr;
  return q;
}

FieldRepairStepResult stepFieldRepair(double dtRealSec,
                                      double hullMax,
                                      FieldRepairInventory& inv,
                                      double& hullCurrent,
                                      const FieldRepairParams& p) {
  FieldRepairStepResult r{};

  const double eps = std::max(0.0, p.eps);

  if (!(dtRealSec > eps) || !(hullMax > eps)) {
    r.done = true;
    return r;
  }

  const double capFrac = clamp01(p.maxRepairFrac);
  r.capHull = hullMax * capFrac;

  if (hullCurrent + eps >= r.capHull) {
    r.reachedCap = true;
    r.done = true;
    return r;
  }

  const double missing = std::max(0.0, r.capHull - hullCurrent);
  const double rate = std::max(0.0, p.hullPerRealSec);
  double want = rate * dtRealSec;
  if (want <= eps) {
    r.done = true;
    return r;
  }
  want = std::min(want, missing);

  const double mPer = std::max(0.0, p.metalsPerHull);
  const double wPer = std::max(0.0, p.machineryPerHull);

  // Compute how much hull we can afford this step.
  double canByMetals = want;
  double canByMach = want;

  if (mPer > eps) canByMetals = std::max(0.0, inv.metals) / mPer;
  if (wPer > eps) canByMach = std::max(0.0, inv.machinery) / wPer;

  const double canByStock = std::max(0.0, std::min(canByMetals, canByMach));
  const double doRepair = std::max(0.0, std::min(want, canByStock));

  if (doRepair <= eps) {
    r.outOfSupplies = true;
    r.done = true;
    return r;
  }

  r.hullRepaired = doRepair;
  r.metalsTaken = doRepair * mPer;
  r.machineryTaken = doRepair * wPer;
  r.heatAdded = doRepair * std::max(0.0, p.heatPerHull);

  // Apply.
  inv.metals = std::max(0.0, inv.metals - r.metalsTaken);
  inv.machinery = std::max(0.0, inv.machinery - r.machineryTaken);

  hullCurrent = std::min(r.capHull, hullCurrent + doRepair);

  r.progressed = true;

  // Determine whether we can progress further.
  if (hullCurrent + eps >= r.capHull) {
    r.reachedCap = true;
    r.done = true;
    return r;
  }

  // If either material is effectively exhausted, we're done.
  double remByMetals = 1.0;
  double remByMach = 1.0;
  if (mPer > eps) remByMetals = inv.metals / mPer;
  if (wPer > eps) remByMach = inv.machinery / wPer;
  const double remByStock = std::min(remByMetals, remByMach);

  if (remByStock <= eps) {
    r.outOfSupplies = true;
    r.done = true;
  }

  return r;
}

} // namespace stellar::sim
