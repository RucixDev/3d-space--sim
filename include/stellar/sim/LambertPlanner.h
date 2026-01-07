#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/LambertSolver.h"

#include <functional>
#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// LambertPlanner — gameplay-friendly transfer window search ("porkchop")
// -----------------------------------------------------------------------------
//
// This module builds on the existing universal-variable Lambert solver and adds
// a deterministic way to search a *window* of departure times and times-of-flight
// to find low-Δv transfers between two moving objects (stations/planets).
//
// In astrodynamics, visualizing Δv over (departure, TOF) forms a "porkchop plot".
// Here we keep it lightweight and gameplay-oriented:
//  - fixed-size sampling grid (bounded CPU)
//  - deterministic best-K extraction
//  - optional full grid capture for tooling/plotting
//
// Units:
//   - timeDays: simulation time in days (for ephemerides)
//   - departAfterSec / tofSec: seconds
//   - positions: km
//   - velocities: km/s

enum class LambertScoreMode : core::u8 {
  // Minimize only the departure burn magnitude (good for fly-bys).
  MinDepartDv = 0,

  // Minimize departure + arrival relative speed (good for rendezvous).
  MinTotalDv = 1,

  // Minimize only arrival relative speed (good for "soft arrival" / capture).
  MinArrivalRelSpeed = 2,

  // Minimize departure + arrivalWeight * arrivalRelSpeed.
  Weighted = 3,
};

// Ephemeris callback: fill outPosKm/outVelKmS for an object at timeDays.
using EphemerisFn = std::function<void(double timeDays, math::Vec3d& outPosKm, math::Vec3d& outVelKmS)>;

// Single evaluated transfer candidate.
struct LambertTransferMetrics {
  bool ok{false};

  // Inputs
  double departAfterSec{0.0};
  double tofSec{0.0};

  // Lambert diagnostics/solution
  LambertResult lambert{};

  // Selected number of complete revolutions (M) for this solution.
  int revolutions{0};

  // Derived burns/metrics
  math::Vec3d dvDepartKmS{0, 0, 0};
  double dvDepartMagKmS{0.0};

  math::Vec3d arriveRelVelKmS{0, 0, 0};
  double arriveRelSpeedKmS{0.0};

  double totalDvKmS{0.0};

  // The scalar score used by the search.
  double score{0.0};
};

// Light-weight grid cell for plotting/debug output.
struct LambertPorkchopCell {
  bool ok{false};
  double departAfterSec{0.0};
  double tofSec{0.0};
  double dvDepartMagKmS{0.0};
  double arriveRelSpeedKmS{0.0};
  double totalDvKmS{0.0};
  double score{0.0};

  // Number of complete revolutions (M) for the selected solution.
  int revolutions{0};
};

struct LambertPorkchopParams {
  // Departure window (seconds relative to baseTimeDays).
  double departMinSec{0.0};
  double departMaxSec{0.0};
  int departSteps{1};

  // Time-of-flight window (seconds).
  double tofMinSec{0.25 * 3600.0};
  double tofMaxSec{48.0 * 3600.0};
  int tofSteps{64};

  LambertScoreMode scoreMode{LambertScoreMode::MinTotalDv};

  // Used only when scoreMode == Weighted.
  double arrivalWeight{0.25};

  // Keep the best K candidates (sorted by score). 0 disables collection.
  int topK{5};

  // If true, populate LambertPorkchopResult::grid.
  bool storeGrid{false};

  // Maximum number of complete revolutions (M) to consider.
  //
  // When >0, the planner will evaluate all feasible solutions for M in [0..maxRevolutions]
  // and keep whichever yields the best score. This can be significantly more expensive
  // than a classic 0-rev solve, so keep grid sizes modest.
  int maxRevolutions{0};

  // Lambert options used for each solve.
  // If refNormal is left as {0,0,0}, the planner will choose a stable reference
  // normal based on the departure state.
  LambertOptions lambertOpt{};
};

struct LambertPorkchopResult {
  // Best candidates sorted by (score, departAfterSec, tofSec).
  std::vector<LambertTransferMetrics> best;

  // Optional grid (row-major: depart outer, tof inner).
  std::vector<LambertPorkchopCell> grid;

  int departSteps{0};
  int tofSteps{0};
};

// -----------------------------------------------------------------------------
// LambertPorkchopStepper — incremental search (frame-friendly)
// -----------------------------------------------------------------------------
//
// `searchLambertPorkchop(...)` is great for headless tooling, but in interactive
// UIs you typically want to spread work across frames.
//
// This stepper evaluates the same fixed-size (departure, TOF) grid in a
// deterministic order and incrementally maintains:
//  - best-K candidates (sorted by score)
//  - optional full-grid capture (if params.storeGrid)
//
// Usage:
//   LambertPorkchopStepper s;
//   s.start(baseTimeDays, depEphem, arrEphem, mu, params);
//   while (!s.done()) s.step(256); // process 256 cells at a time
//   auto res = s.result();
class LambertPorkchopStepper {
public:
  LambertPorkchopStepper() = default;

  // Start a new search, resetting any previous state.
  void start(double baseTimeDays,
             EphemerisFn departEphem,
             EphemerisFn arriveEphem,
             double muKm3S2,
             const LambertPorkchopParams& params);

  // Cancel / reset to an empty state.
  void reset();

  bool started() const { return started_; }
  bool done() const { return started_ && done_; }

  // 0..1 progress through the grid (best-effort; returns 1 when done).
  double progress01() const;

  // Process up to maxCells additional grid cells. Returns true if the search
  // is now complete.
  bool step(int maxCells = 256);

  // Current partial/final result.
  const LambertPorkchopResult& result() const { return res_; }

private:
  // Sanitized inputs
  bool started_{false};
  bool done_{false};

  double baseTimeDays_{0.0};
  EphemerisFn departEphem_{};
  EphemerisFn arriveEphem_{};
  double muKm3S2_{0.0};
  LambertPorkchopParams params_{};

  int depSteps_{0};
  int tofSteps_{0};
  double depMin_{0.0};
  double depMax_{0.0};
  double tofMin_{0.0};
  double tofMax_{0.0};

  int iDep_{0};
  int iTof_{0};
  int processed_{0};
  int total_{0};

  LambertPorkchopResult res_{};
};

// Evaluate a single candidate transfer.
LambertTransferMetrics evaluateLambertTransfer(double baseTimeDays,
                                               double departAfterSec,
                                               double tofSec,
                                               const EphemerisFn& departEphem,
                                               const EphemerisFn& arriveEphem,
                                               double muKm3S2,
                                               const LambertOptions& opt,
                                               LambertScoreMode scoreMode,
                                               double arrivalWeight = 0.25,
                                               int maxRevolutions = 0);

// Search a (departure, TOF) window for the best transfer(s).
LambertPorkchopResult searchLambertPorkchop(double baseTimeDays,
                                            const EphemerisFn& departEphem,
                                            const EphemerisFn& arriveEphem,
                                            double muKm3S2,
                                            const LambertPorkchopParams& params);

} // namespace stellar::sim
