#include "stellar/sim/LambertPlanner.h"

#include "stellar/sim/Units.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::sim {
namespace {

inline bool vecNearlyZero(const math::Vec3d& v, double eps = 1e-12) {
  return v.lengthSq() <= eps * eps;
}

inline double lerp(double a, double b, double t) {
  return a + (b - a) * t;
}

struct SanitizedPorkchopParams {
  int depSteps{1};
  int tofSteps{1};
  double depMin{0.0};
  double depMax{0.0};
  double tofMin{1.0};
  double tofMax{1.0};
};

SanitizedPorkchopParams sanitize(const LambertPorkchopParams& params) {
  SanitizedPorkchopParams s;
  s.depSteps = std::max(1, params.departSteps);
  s.tofSteps = std::max(1, params.tofSteps);

  s.depMin = params.departMinSec;
  s.depMax = params.departMaxSec;
  if (!std::isfinite(s.depMin)) s.depMin = 0.0;
  if (!std::isfinite(s.depMax)) s.depMax = s.depMin;
  if (s.depMax < s.depMin) std::swap(s.depMin, s.depMax);

  s.tofMin = params.tofMinSec;
  s.tofMax = params.tofMaxSec;
  if (!std::isfinite(s.tofMin)) s.tofMin = 1.0;
  if (!std::isfinite(s.tofMax)) s.tofMax = s.tofMin;
  if (s.tofMax < s.tofMin) std::swap(s.tofMin, s.tofMax);
  s.tofMin = std::max(1e-3, s.tofMin);
  s.tofMax = std::max(s.tofMin, s.tofMax);

  return s;
}

LambertOptions bakeOptions(const LambertOptions& in,
                           const math::Vec3d& r1Km,
                           const math::Vec3d& v0KmS,
                           const math::Vec3d& r2Km) {
  LambertOptions o = in;
  if (vecNearlyZero(o.refNormal)) {
    math::Vec3d n = math::cross(r1Km, v0KmS);
    if (vecNearlyZero(n)) n = math::cross(r1Km, r2Km);
    if (vecNearlyZero(n)) n = {0, 0, 1};
    o.refNormal = n;
  }
  return o;
}

double computeScore(LambertScoreMode mode,
                    double dvDepartKmS,
                    double arriveRelSpeedKmS,
                    double arrivalWeight) {
  dvDepartKmS = std::max(0.0, dvDepartKmS);
  arriveRelSpeedKmS = std::max(0.0, arriveRelSpeedKmS);

  switch (mode) {
    case LambertScoreMode::MinDepartDv:
      return dvDepartKmS;
    case LambertScoreMode::MinTotalDv:
      return dvDepartKmS + arriveRelSpeedKmS;
    case LambertScoreMode::MinArrivalRelSpeed:
      return arriveRelSpeedKmS;
    case LambertScoreMode::Weighted:
    default:
      arrivalWeight = std::max(0.0, arrivalWeight);
      return dvDepartKmS + arrivalWeight * arriveRelSpeedKmS;
  }
}

void insertBest(std::vector<LambertTransferMetrics>& best,
                const LambertTransferMetrics& cand,
                int topK) {
  if (topK <= 0) return;

  auto cmp = [](const LambertTransferMetrics& a, const LambertTransferMetrics& b) {
    if (a.score != b.score) return a.score < b.score;
    if (a.departAfterSec != b.departAfterSec) return a.departAfterSec < b.departAfterSec;
    return a.tofSec < b.tofSec;
  };

  const auto it = std::lower_bound(best.begin(), best.end(), cand, cmp);
  best.insert(it, cand);

  if ((int)best.size() > topK) {
    best.resize((std::size_t)topK);
  }
}

} // namespace

void LambertPorkchopStepper::start(double baseTimeDays,
                                   EphemerisFn departEphem,
                                   EphemerisFn arriveEphem,
                                   double muKm3S2,
                                   const LambertPorkchopParams& params) {
  reset();

  started_ = true;
  done_ = false;

  baseTimeDays_ = baseTimeDays;
  departEphem_ = std::move(departEphem);
  arriveEphem_ = std::move(arriveEphem);
  muKm3S2_ = muKm3S2;
  params_ = params;

  const SanitizedPorkchopParams sp = sanitize(params_);
  depSteps_ = sp.depSteps;
  tofSteps_ = sp.tofSteps;
  depMin_ = sp.depMin;
  depMax_ = sp.depMax;
  tofMin_ = sp.tofMin;
  tofMax_ = sp.tofMax;

  res_ = LambertPorkchopResult{};
  res_.departSteps = depSteps_;
  res_.tofSteps = tofSteps_;

  const int topK = std::max(0, params_.topK);
  if (topK > 0) {
    res_.best.reserve((std::size_t)topK);
  }
  if (params_.storeGrid) {
    res_.grid.reserve((std::size_t)depSteps_ * (std::size_t)tofSteps_);
  }

  iDep_ = 0;
  iTof_ = 0;
  processed_ = 0;
  total_ = depSteps_ * tofSteps_;

  // If anything is clearly invalid, mark done immediately.
  if (!std::isfinite(baseTimeDays_) || !(muKm3S2_ > 0.0) || !std::isfinite(muKm3S2_) ||
      !departEphem_ || !arriveEphem_ || total_ <= 0) {
    done_ = true;
  }
}

void LambertPorkchopStepper::reset() {
  started_ = false;
  done_ = false;
  baseTimeDays_ = 0.0;
  departEphem_ = {};
  arriveEphem_ = {};
  muKm3S2_ = 0.0;
  params_ = LambertPorkchopParams{};
  depSteps_ = 0;
  tofSteps_ = 0;
  depMin_ = 0.0;
  depMax_ = 0.0;
  tofMin_ = 0.0;
  tofMax_ = 0.0;
  iDep_ = 0;
  iTof_ = 0;
  processed_ = 0;
  total_ = 0;
  res_ = LambertPorkchopResult{};
}

double LambertPorkchopStepper::progress01() const {
  if (!started_) return 0.0;
  if (done_) return 1.0;
  if (total_ <= 0) return 0.0;
  const double p = (double)processed_ / (double)total_;
  if (!std::isfinite(p)) return 0.0;
  return std::clamp(p, 0.0, 1.0);
}

bool LambertPorkchopStepper::step(int maxCells) {
  if (!started_ || done_) return done_;
  if (maxCells <= 0) maxCells = 1;

  const int topK = std::max(0, params_.topK);

  for (int k = 0; k < maxCells && !done_; ++k) {
    if (iDep_ >= depSteps_) {
      done_ = true;
      break;
    }

    const double tDep = (depSteps_ == 1) ? depMin_
                                         : lerp(depMin_, depMax_, (double)iDep_ / (double)(depSteps_ - 1));
    const double tof = (tofSteps_ == 1) ? tofMin_
                                        : lerp(tofMin_, tofMax_, (double)iTof_ / (double)(tofSteps_ - 1));

    LambertTransferMetrics cand = evaluateLambertTransfer(baseTimeDays_,
                                                          tDep,
                                                          tof,
                                                          departEphem_,
                                                          arriveEphem_,
                                                          muKm3S2_,
                                                          params_.lambertOpt,
                                                          params_.scoreMode,
                                                          params_.arrivalWeight,
                                                          params_.maxRevolutions);

    if (params_.storeGrid) {
      LambertPorkchopCell cell{};
      cell.ok = cand.ok;
      cell.departAfterSec = tDep;
      cell.tofSec = tof;
      cell.dvDepartMagKmS = cand.dvDepartMagKmS;
      cell.arriveRelSpeedKmS = cand.arriveRelSpeedKmS;
      cell.totalDvKmS = cand.totalDvKmS;
      cell.score = cand.score;
      cell.revolutions = cand.revolutions;
      res_.grid.push_back(cell);
    }

    if (cand.ok) {
      insertBest(res_.best, cand, topK);
    }

    ++processed_;

    // Advance
    ++iTof_;
    if (iTof_ >= tofSteps_) {
      iTof_ = 0;
      ++iDep_;
      if (iDep_ >= depSteps_) {
        done_ = true;
        break;
      }
    }
  }

  return done_;
}

LambertTransferMetrics evaluateLambertTransfer(double baseTimeDays,
                                               double departAfterSec,
                                               double tofSec,
                                               const EphemerisFn& departEphem,
                                               const EphemerisFn& arriveEphem,
                                               double muKm3S2,
                                               const LambertOptions& opt,
                                               LambertScoreMode scoreMode,
                                               double arrivalWeight,
                                               int maxRevolutions) {
  LambertTransferMetrics out{};
  out.departAfterSec = departAfterSec;
  out.tofSec = tofSec;

  if (!departEphem || !arriveEphem) {
    return out;
  }

  if (!(muKm3S2 > 0.0) || !std::isfinite(muKm3S2)) {
    return out;
  }
  if (!std::isfinite(baseTimeDays) || !std::isfinite(departAfterSec) || !std::isfinite(tofSec)) {
    return out;
  }

  const double safeTofSec = std::max(1e-3, tofSec);
  const double departDays = baseTimeDays + departAfterSec / kSecondsPerDay;
  const double arriveDays = departDays + safeTofSec / kSecondsPerDay;

  math::Vec3d r1Km{0, 0, 0};
  math::Vec3d v0KmS{0, 0, 0};
  math::Vec3d r2Km{0, 0, 0};
  math::Vec3d vTgtKmS{0, 0, 0};

  departEphem(departDays, r1Km, v0KmS);
  arriveEphem(arriveDays, r2Km, vTgtKmS);

  LambertOptions baked = bakeOptions(opt, r1Km, v0KmS, r2Km);

  // Optionally consider multi-revolution solutions and keep whichever yields the best score.
  if (maxRevolutions > 0) {
    const LambertMultiRevResult multi = solveLambertUniversalMultiRev(r1Km, r2Km, safeTofSec, muKm3S2, maxRevolutions, baked);

    bool bestOk = false;
    double bestScore = std::numeric_limits<double>::infinity();
    LambertResult bestLambert{};
    math::Vec3d bestDvDepart{0, 0, 0};
    math::Vec3d bestArrRel{0, 0, 0};
    double bestDvDepartMag = 0.0;
    double bestArrRelSpeed = 0.0;

    for (const auto& sol : multi.solutions) {
      if (!sol.ok) continue;

      const math::Vec3d dvDepart = sol.v1KmS - v0KmS;
      const double dvDepartMag = dvDepart.length();
      const math::Vec3d arrRelVel = sol.v2KmS - vTgtKmS;
      const double arrRelSpeed = arrRelVel.length();
      const double score = computeScore(scoreMode, dvDepartMag, arrRelSpeed, arrivalWeight);

      if (!std::isfinite(score)) continue;

      if (!bestOk || score < bestScore) {
        bestOk = true;
        bestScore = score;
        bestLambert = sol;
        bestDvDepart = dvDepart;
        bestArrRel = arrRelVel;
        bestDvDepartMag = dvDepartMag;
        bestArrRelSpeed = arrRelSpeed;
      }
    }

    if (!bestOk) {
      return out;
    }

    out.lambert = bestLambert;
    out.ok = true;
    out.revolutions = bestLambert.revolutions;
    out.dvDepartKmS = bestDvDepart;
    out.dvDepartMagKmS = bestDvDepartMag;
    out.arriveRelVelKmS = bestArrRel;
    out.arriveRelSpeedKmS = bestArrRelSpeed;
  } else {
    out.lambert = solveLambertUniversal(r1Km, r2Km, safeTofSec, muKm3S2, baked);
    out.ok = out.lambert.ok;
    if (!out.ok) {
      return out;
    }

    out.revolutions = out.lambert.revolutions;
    out.dvDepartKmS = out.lambert.v1KmS - v0KmS;
    out.dvDepartMagKmS = out.dvDepartKmS.length();

    out.arriveRelVelKmS = out.lambert.v2KmS - vTgtKmS;
    out.arriveRelSpeedKmS = out.arriveRelVelKmS.length();
  }

  out.totalDvKmS = out.dvDepartMagKmS + out.arriveRelSpeedKmS;
  out.score = computeScore(scoreMode, out.dvDepartMagKmS, out.arriveRelSpeedKmS, arrivalWeight);
  if (!std::isfinite(out.score)) {
    out.ok = false;
  }

  return out;
}

LambertPorkchopResult searchLambertPorkchop(double baseTimeDays,
                                            const EphemerisFn& departEphem,
                                            const EphemerisFn& arriveEphem,
                                            double muKm3S2,
                                            const LambertPorkchopParams& params) {
  LambertPorkchopStepper stepper;
  stepper.start(baseTimeDays, departEphem, arriveEphem, muKm3S2, params);
  // Large chunk size; callers that need frame-friendly stepping should use the
  // stepper directly.
  while (!stepper.done()) {
    stepper.step(4096);
  }
  return stepper.result();
}

} // namespace stellar::sim
