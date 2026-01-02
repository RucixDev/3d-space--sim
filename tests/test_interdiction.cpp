#include "stellar/sim/Interdiction.h"
#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

int test_interdiction() {
  using stellar::core::u64;
  namespace sim = stellar::sim;

  int fails = 0;

  const u64 seed = 123;
  const u64 pirateId = 999;
  const double timeDays = 10.0;
  const stellar::math::Vec3d fwd{0.0, 0.0, 1.0};
  const stellar::math::Vec3d up{0.0, 1.0, 0.0};

  // 1) Warning phase should transition into Active once the timer elapses.
  {
    sim::InterdictionParams p{};
    p.warningSec = 0.5;
    p.activeSec = 10.0;
    p.startMeter = 0.5;
    // Freeze the meter so this test only validates phase transition logic.
    p.playerGainPerSec = 0.0;
    p.pullBasePerSec = 0.0;
    p.pullExtraPerSec = 0.0;
    p.escapeDriftRadPerSec = 0.0;
    p.escapeJitterRadPerSec = 0.0;

    auto st = sim::beginInterdiction(seed, pirateId, timeDays, fwd, up, /*closeness=*/0.5, /*strength=*/1.0, p);
    if (st.phase != sim::InterdictionPhase::Warning) {
      std::fprintf(stderr, "[test_interdiction] expected to start in Warning phase.\n");
      ++fails;
    }

    bool beganActive = false;
    for (int i = 0; i < 40 && sim::interdictionInProgress(st); ++i) {
      auto out = sim::stepInterdiction(st, /*dt=*/0.05, fwd, /*closeness=*/0.5, /*strength=*/1.0, /*submit=*/false, p);
      beganActive |= out.beganActive;
    }

    if (!beganActive || st.phase != sim::InterdictionPhase::Active) {
      std::fprintf(stderr, "[test_interdiction] expected to transition to Active after warning.\n");
      ++fails;
    }
  }

  // 2) Strong alignment with zero pull should eventually evade.
  {
    sim::InterdictionParams p{};
    p.warningSec = 0.0;
    p.activeSec = 10.0;
    p.startMeter = 0.2;
    p.playerGainPerSec = 0.8;
    p.pullBasePerSec = 0.0;
    p.pullExtraPerSec = 0.0;
    p.escapeDriftRadPerSec = 0.0;
    p.escapeJitterRadPerSec = 0.0;

    auto st = sim::beginInterdiction(seed, pirateId, timeDays, fwd, up, /*closeness=*/0.5, /*strength=*/1.0, p);
    // Force a deterministic, perfectly alignable escape vector for this test.
    st.escapeDir = fwd;
    st.phase = sim::InterdictionPhase::Active;
    st.activeRemainingSec = p.activeSec;
    st.meter = p.startMeter;

    bool evaded = false;
    for (int i = 0; i < 200 && sim::interdictionInProgress(st); ++i) {
      auto out = sim::stepInterdiction(st, /*dt=*/0.1, fwd, /*closeness=*/0.5, /*strength=*/1.0, /*submit=*/false, p);
      if (out.evaded) {
        evaded = true;
        break;
      }
    }

    if (!evaded) {
      std::fprintf(stderr, "[test_interdiction] expected to evade with strong alignment.\n");
      ++fails;
    }
  }

  // 3) Poor alignment with strong pull should fail.
  {
    sim::InterdictionParams p{};
    p.warningSec = 0.0;
    p.activeSec = 10.0;
    p.startMeter = 0.8;
    p.playerGainPerSec = 0.0;
    p.pullBasePerSec = 1.0;
    p.pullExtraPerSec = 0.0;
    p.escapeDriftRadPerSec = 0.0;
    p.escapeJitterRadPerSec = 0.0;

    auto st = sim::beginInterdiction(seed, pirateId, timeDays, fwd, up, /*closeness=*/1.0, /*strength=*/1.0, p);
    st.escapeDir = fwd;
    st.phase = sim::InterdictionPhase::Active;
    st.activeRemainingSec = p.activeSec;
    st.meter = p.startMeter;

    const stellar::math::Vec3d badFwd{0.0, 0.0, -1.0};
    bool failed = false;
    for (int i = 0; i < 200 && sim::interdictionInProgress(st); ++i) {
      auto out = sim::stepInterdiction(st, /*dt=*/0.1, badFwd, /*closeness=*/1.0, /*strength=*/1.0, /*submit=*/false, p);
      if (out.failed) {
        failed = true;
        break;
      }
    }

    if (!failed) {
      std::fprintf(stderr, "[test_interdiction] expected to fail with poor alignment + strong pull.\n");
      ++fails;
    }
  }

  return fails;
}
