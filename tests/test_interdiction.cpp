#include "stellar/sim/Interdiction.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

static bool approxVec(const math::Vec3d& a, const math::Vec3d& b, double eps = 1e-9) {
  return approx(a.x, b.x, eps) && approx(a.y, b.y, eps) && approx(a.z, b.z, eps);
}

int test_interdiction() {
  int fails = 0;

  // Deterministic begin: same inputs => same escape vector.
  {
    sim::InterdictionParams params;
    const core::u64 seed = 0xBADC0FFEEull;
    const core::u64 pirateId = 12345;
    const double timeDays = 42.5;
    const math::Vec3d fwd{0.0, 0.0, 1.0};
    const math::Vec3d up{0.0, 1.0, 0.0};

    const auto a = sim::beginInterdiction(seed, pirateId, timeDays, fwd, up, 0.6, 1.2, params);
    const auto b = sim::beginInterdiction(seed, pirateId, timeDays, fwd, up, 0.6, 1.2, params);

    if (a.phase != sim::InterdictionPhase::Warning || b.phase != sim::InterdictionPhase::Warning) {
      std::cerr << "[test_interdiction] expected beginInterdiction() to start in Warning phase.\n";
      ++fails;
    }

    if (!approxVec(a.escapeDir, b.escapeDir, 1e-12)) {
      std::cerr << "[test_interdiction] beginInterdiction() not deterministic for same inputs.\n";
      ++fails;
    }

    const double dot = math::dot(fwd.normalized(), a.escapeDir.normalized());
    if (dot < 0.10 || dot > 0.999) {
      std::cerr << "[test_interdiction] escapeDir dot forward out of expected range: " << dot << "\n";
      ++fails;
    }
  }

  // Evade path: with perfect alignment and weak pull, meter should reach 1.0.
  {
    sim::InterdictionParams params;
    params.escapeDriftRadPerSec = 0.0;
    params.escapeJitterRadPerSec = 0.0;

    sim::InterdictionState s = sim::beginInterdiction(/*seed=*/777, /*pirateId=*/9,
                                                      /*timeDays=*/100.25,
                                                      /*fwd=*/{0,0,1}, /*up=*/{0,1,0},
                                                      /*close=*/0.15, /*strength=*/0.7,
                                                      params);

    // Advance through warning into active.
    sim::stepInterdiction(s, params.warningSec + 0.01, {0,0,1}, 0.15, 0.7, false, params);
    if (s.phase != sim::InterdictionPhase::Active) {
      std::cerr << "[test_interdiction] expected to enter Active phase after warning.\n";
      ++fails;
    }

    bool evaded = false;
    for (int i = 0; i < 500 && sim::interdictionInProgress(s); ++i) {
      const math::Vec3d fwd = s.escapeDir; // perfect alignment
      const auto out = sim::stepInterdiction(s, 0.05, fwd, 0.15, 0.7, false, params);
      if (out.failed) {
        std::cerr << "[test_interdiction] unexpected failure during evade path.\n";
        ++fails;
        break;
      }
      if (out.evaded) {
        evaded = true;
        break;
      }
    }
    if (!evaded) {
      std::cerr << "[test_interdiction] expected to evade with strong alignment.\n";
      ++fails;
    }
  }

  // Fail path: with poor alignment and strong pull, meter should hit 0.
  {
    sim::InterdictionParams params;
    params.escapeDriftRadPerSec = 0.0;
    params.escapeJitterRadPerSec = 0.0;

    sim::InterdictionState s = sim::beginInterdiction(/*seed=*/888, /*pirateId=*/77,
                                                      /*timeDays=*/3.0,
                                                      /*fwd=*/{0,0,1}, /*up=*/{0,1,0},
                                                      /*close=*/1.0, /*strength=*/1.6,
                                                      params);
    sim::stepInterdiction(s, params.warningSec + 0.01, {0,0,1}, 1.0, 1.6, false, params);

    bool failed = false;
    for (int i = 0; i < 800 && sim::interdictionInProgress(s); ++i) {
      const math::Vec3d anti{0, 0, -1};
      const auto out = sim::stepInterdiction(s, 0.05, anti, 1.0, 1.6, false, params);
      if (out.evaded) {
        std::cerr << "[test_interdiction] unexpected evade during fail path.\n";
        ++fails;
        break;
      }
      if (out.failed) {
        failed = true;
        break;
      }
    }

    if (!failed) {
      std::cerr << "[test_interdiction] expected to fail with poor alignment + strong pull.\n";
      ++fails;
    }
  }

  // Submit should end immediately.
  {
    sim::InterdictionParams params;
    sim::InterdictionState s = sim::beginInterdiction(/*seed=*/999, /*pirateId=*/1,
                                                      /*timeDays=*/1.0,
                                                      /*fwd=*/{0,0,1}, /*up=*/{0,1,0},
                                                      /*close=*/0.5, /*strength=*/1.0,
                                                      params);

    const auto out = sim::stepInterdiction(s, 0.016, {0,0,1}, 0.5, 1.0, /*submit=*/true, params);
    if (!out.submitted || !out.endedThisFrame || sim::interdictionInProgress(s)) {
      std::cerr << "[test_interdiction] submit should end interdiction immediately.\n";
      ++fails;
    }
  }

  // Trigger chance should respond to inputs.
  {
    sim::InterdictionTriggerParams tp;
    const double dt = 0.2;
    const double lo = sim::interdictionTriggerChance(dt, 0.1, 0.0, 1.0, tp);
    const double hi = sim::interdictionTriggerChance(dt, 0.9, 12'000.0, 1.5, tp);
    if (!(hi > lo)) {
      std::cerr << "[test_interdiction] triggerChance should increase with closeness/cargo/strength.\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_interdiction] PASS\n";
  }
  return fails;
}
