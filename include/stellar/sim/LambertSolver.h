#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Lambert solver (universal variables, 0-rev)
// -----------------------------------------------------------------------------
//
// Solves the classic two-body Lambert boundary value problem:
//   given r1, r2, dt, mu -> find v1, v2 such that the Kepler arc from
//   (r1, v1) reaches r2 after dt.
//
// Notes:
//  - This is a *single-revolution* (0-rev) solver intended for gameplay tools.
//  - Uses a robust bisection on the universal variable z.
//  - Units:
//      r: km
//      v: km/s
//      mu: km^3/s^2
//      dt: seconds

struct LambertOptions {
  // When true, pick the solution with transfer angle > pi.
  // (Often called the "indirect" or "long-way" solution.)
  bool longWay{false};

  // If refNormal is non-zero, prograde chooses the transfer that preserves
  // the sign of dot(r1 x r2, refNormal). Retrograde flips it.
  bool prograde{true};

  // Reference plane normal used to disambiguate the direction of motion.
  // If left as {0,0,0}, the solver will use (r1 x r2) as the reference.
  math::Vec3d refNormal{0, 0, 0};

  int maxIterations{64};

  // Time-of-flight tolerance (seconds) for convergence.
  double tolSec{1e-4};
};

struct LambertResult {
  bool ok{false};

  math::Vec3d v1KmS{0, 0, 0};
  math::Vec3d v2KmS{0, 0, 0};

  // Transfer angle used (radians). In [0, 2pi].
  double transferAngleRad{0.0};

  // Diagnostics
  double z{0.0};
  int iterations{0};

  // Number of complete revolutions (M) of the transfer orbit.
  // 0 corresponds to the classic single-revolution solution.
  int revolutions{0};
};

LambertResult solveLambertUniversal(const math::Vec3d& r1Km,
                                   const math::Vec3d& r2Km,
                                   double dtSec,
                                   double muKm3S2,
                                   const LambertOptions& opt = {});


// Find *all* feasible Lambert solutions up to maxRevolutions.
//
// For a given geometry and time-of-flight there can be multiple feasible
// multi-revolution solutions (typically up to two per M>0 branch).
// The returned solutions are sorted by (revolutions, z).
struct LambertMultiRevResult {
  bool ok{false};
  std::vector<LambertResult> solutions;
};

LambertMultiRevResult solveLambertUniversalMultiRev(const math::Vec3d& r1Km,
                                                    const math::Vec3d& r2Km,
                                                    double dtSec,
                                                    double muKm3S2,
                                                    int maxRevolutions,
                                                    const LambertOptions& opt = {});

} // namespace stellar::sim
