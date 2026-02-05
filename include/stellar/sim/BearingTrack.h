#pragma once

#include "stellar/math/Mat3.h"
#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// BearingTrack — bearing-only 3D triangulation (headless, deterministic)
// -----------------------------------------------------------------------------
//
// Motivation:
//   A lot of interesting "space sim" sensing situations are bearing-only:
//     - passive EM/radiation bearings (no range)
//     - intermittent contacts (range drops out under jamming/occlusion)
//     - mission clues where the player must "cross-fix" a signal
//
//   Given multiple observer positions {o_i} and corresponding unit bearings
//   {d_i} pointing at the same target, we can estimate the 3D point p that
//   best fits all lines in a least-squares sense.
//
// Math:
//   We minimize sum_i || (I - d_i d_i^T) (p - o_i) ||^2.
//   This yields a symmetric linear system:
//     A p = b
//   where:
//     A = sum_i w_i (I - d_i d_i^T)
//     b = sum_i w_i (I - d_i d_i^T) o_i
//
// Design constraints:
//   - header-only
//   - deterministic
//   - no dynamic allocations
//
// Notes:
//   - With a single bearing, A is rank-2 (unsolvable for unique p). The track
//     remains invalid until enough geometric diversity exists.
//   - We support exponential forgetting (half-life) so the estimate adapts when
//     the target or observer moves.

struct BearingTrackParams {
  // Exponential forgetting for the normal-equation accumulators.
  // <= 0 disables decay (infinite memory).
  double observationHalfLifeSec{6.0};

  // Small diagonal regularization added when solving. This helps avoid numeric
  // blowups when bearings are near-parallel, but *does not* override the
  // determinant gate (we still require the unregularized system to be solvable).
  double solveRegularization{1.0e-6};

  // Relative determinant threshold used to gate solvability.
  // We scale this by (maxAbsElement^3) so it behaves sensibly under different
  // weight magnitudes.
  double determinantEps{1.0e-12};

  // Require at least this total (decayed) weight before attempting a solve.
  // This helps avoid "solving" from a single overweight measurement.
  double minEffectiveWeight{1.8};

  // Velocity estimation smoothing (half-life). <= 0 snaps to measured velocity.
  double velHalfLifeSec{0.65};

  // Clamp the velocity magnitude (safety against pathological solves).
  double maxSpeedKmS{25000.0};

  // Residual-derived 1-sigma clamp.
  double sigmaMinKm{0.0};
  double sigmaMaxKm{800000.0};
};

struct BearingTrack3d {
  // Normal equation accumulators.
  math::Mat3d A{};
  math::Vec3d b{0, 0, 0};
  double weight{0.0};

  // Residual accumulator (fit quality): sigma^2 ~= E[dist(line, p)^2].
  double residualWeight{0.0};
  double residualSqSum{0.0};

  // Current estimate (km / km/s).
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
  double sigmaKm{0.0};

  // Seconds since last used measurement.
  double ageSinceMeasSec{0.0};

  bool initialized{false};
};

struct BearingTrackResult {
  bool valid{false};
  bool hasMeasurement{false};

  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
  double sigmaKm{0.0};
  double ageSinceMeasSec{0.0};
};

inline double halfLifeDecay(double halfLifeSec, double dtSec) {
  if (!std::isfinite(dtSec) || dtSec <= 0.0) return 1.0;
  if (!std::isfinite(halfLifeSec) || halfLifeSec <= 0.0) return 1.0;
  return std::pow(0.5, dtSec / halfLifeSec);
}

inline double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

inline double lerp(double a, double b, double t) {
  return a + (b - a) * t;
}

inline math::Vec3d lerp(const math::Vec3d& a, const math::Vec3d& b, double t) {
  return a * (1.0 - t) + b * t;
}

inline double maxAbsElement(const math::Mat3d& M) {
  double m = 0.0;
  for (double v : M.m) m = std::max(m, std::abs(v));
  return m;
}

// Update the bearing-only track.
//
// Inputs:
//  - dtSec: time step (seconds)
//  - hasMeasurement: if false, the track simply coasts and ages/decays
//  - observerPosKm: observer position in world space (km)
//  - bearingDirWorld: (not necessarily normalized) direction in world space
//  - measWeight: relative confidence (>= 0). Typical range: 0..1.
inline BearingTrackResult updateBearingTrack(BearingTrack3d& track,
                                             double dtSec,
                                             bool hasMeasurement,
                                             const math::Vec3d& observerPosKm,
                                             const math::Vec3d& bearingDirWorld,
                                             double measWeight = 1.0,
                                             const BearingTrackParams& params = {}) {
  if (!std::isfinite(dtSec) || dtSec < 0.0) dtSec = 0.0;

  const math::Vec3d prevPos = track.posKm;

  // Predict/coast.
  if (track.initialized && dtSec > 0.0) {
    track.posKm += track.velKmS * dtSec;
  }

  // Age always advances.
  if (dtSec > 0.0) track.ageSinceMeasSec += dtSec;

  // Exponential decay of accumulators.
  const double decay = halfLifeDecay(params.observationHalfLifeSec, dtSec);
  if (decay != 1.0) {
    track.A *= decay;
    track.b *= decay;
    track.weight *= decay;
    track.residualWeight *= decay;
    track.residualSqSum *= decay;
  }

  bool usedMeas = false;
  math::Vec3d dUnit{0, 0, 0};
  double wUsed = 0.0;

  if (hasMeasurement && math::isFinite(observerPosKm) && math::isFinite(bearingDirWorld)) {
    dUnit = math::safeNormalized(bearingDirWorld, {0, 0, 0}, 1e-18);

    if (dUnit.lengthSq() > 1e-18) {
      wUsed = std::max(0.0, measWeight);
      if (std::isfinite(wUsed) && wUsed > 0.0) {
        usedMeas = true;

        // P = I - d d^T (projector onto the plane perpendicular to d).
        const math::Mat3d ddT = math::Mat3d::outerProduct(dUnit, dUnit);
        math::Mat3d P = math::Mat3d::identity() - ddT;

        track.A += P * wUsed;
        track.b += (P * observerPosKm) * wUsed;
        track.weight += wUsed;

        track.ageSinceMeasSec = 0.0;
      }
    }
  }

  // Determine if the system is geometrically solvable.
  bool canSolve = false;
  if (track.weight >= params.minEffectiveWeight) {
    const double scale = std::max(1.0, maxAbsElement(track.A));
    double detEps = params.determinantEps;
    if (!std::isfinite(detEps) || detEps < 0.0) detEps = 0.0;

    const double detThresh = detEps * scale * scale * scale;
    const double detAbs = std::abs(track.A.determinant());
    if (detAbs > detThresh) {
      canSolve = true;
    }
  }

  bool solved = false;
  math::Vec3d solvedPos = track.posKm;

  if (canSolve) {
    math::Mat3d Areg = track.A;

    // Add gentle diagonal regularization for numeric stability.
    double lambda = params.solveRegularization;
    if (!std::isfinite(lambda) || lambda < 0.0) lambda = 0.0;

    const double tr = std::abs(Areg.trace());
    lambda *= std::max(1.0, tr / 3.0);
    if (lambda > 0.0) Areg.addToDiagonal(lambda);

    // Use an absolute epsilon derived from the matrix scale.
    const double scale = std::max(1.0, maxAbsElement(Areg));
    double detEps = params.determinantEps;
    if (!std::isfinite(detEps) || detEps < 0.0) detEps = 0.0;
    const double detEpsAbs = std::max(1.0e-18, detEps * scale * scale * scale);

    math::Vec3d p{0, 0, 0};
    if (Areg.solve(p, track.b, detEpsAbs) && math::isFinite(p)) {
      solved = true;
      solvedPos = p;
    }
  }

  // Update state when solvable.
  if (solved) {
    // Velocity estimate from successive position solutions.
    if (track.initialized && dtSec > 1.0e-6) {
      const math::Vec3d measVel = (solvedPos - prevPos) / dtSec;

      const double velAlpha = 1.0 - halfLifeDecay(params.velHalfLifeSec, dtSec);
      const double a = clamp01(velAlpha);
      track.velKmS = lerp(track.velKmS, measVel, a);
    } else if (!track.initialized) {
      track.velKmS = {0, 0, 0};
    }

    if (params.maxSpeedKmS > 0.0 && std::isfinite(params.maxSpeedKmS)) {
      track.velKmS = math::clampMagnitude(track.velKmS, params.maxSpeedKmS);
    }

    track.posKm = solvedPos;
    track.initialized = true;
  }

  // Update residual fit quality whenever we have an estimate and a measurement.
  if (track.initialized && usedMeas) {
    const math::Vec3d rel = (track.posKm - observerPosKm);
    // Distance from point to line: |(p - o) x d| (d must be unit).
    const double dist = math::cross(rel, dUnit).length();
    if (std::isfinite(dist)) {
      track.residualSqSum += wUsed * dist * dist;
      track.residualWeight += wUsed;
    }
  }

  if (track.residualWeight > 1e-12) {
    const double sigma2 = std::max(0.0, track.residualSqSum / track.residualWeight);
    track.sigmaKm = std::sqrt(sigma2);
  }

  if (!std::isfinite(track.sigmaKm)) track.sigmaKm = params.sigmaMaxKm;
  track.sigmaKm = std::clamp(track.sigmaKm, params.sigmaMinKm, params.sigmaMaxKm);

  BearingTrackResult out{};
  out.valid = track.initialized;
  out.hasMeasurement = usedMeas;
  out.posKm = track.posKm;
  out.velKmS = track.velKmS;
  out.sigmaKm = track.sigmaKm;
  out.ageSinceMeasSec = track.ageSinceMeasSec;
  return out;
}

} // namespace stellar::sim
