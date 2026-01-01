#include "stellar/sim/Distress.h"

#include <cmath>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

int test_distress() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  const core::u64 seed = 424242ull;
  const SystemId sys = 12345ull;
  const core::u64 sig = 999000111ull;
  const core::u32 faction = 7;
  const double t = 1000.25;

  // Determinism.
  {
    const DistressPlan a = planDistressEncounter(seed, sys, sig, t, faction);
    const DistressPlan b = planDistressEncounter(seed, sys, sig, t, faction);

    if (a.scenario != b.scenario || a.hasVictim != b.hasVictim || a.ambush != b.ambush || a.pirateCount != b.pirateCount) {
      std::cerr << "[test_distress] plan not deterministic (flags)\n";
      ++fails;
    }
    if (a.needCommodity != b.needCommodity) {
      std::cerr << "[test_distress] plan not deterministic (commodity)\n";
      ++fails;
    }
    if (!nearly(a.needUnits, b.needUnits) || !nearly(a.rewardCr, b.rewardCr) || !nearly(a.repReward, b.repReward) || !nearly(a.risk, b.risk)) {
      std::cerr << "[test_distress] plan not deterministic (numeric)\n";
      ++fails;
    }
    if (a.payerFactionId != b.payerFactionId) {
      std::cerr << "[test_distress] plan not deterministic (payer)\n";
      ++fails;
    }
  }

  // Basic bounds.
  {
    const DistressPlan p = planDistressEncounter(seed, sys, sig, t, faction);

    if (p.risk < 0.0 || p.risk > 1.0) {
      std::cerr << "[test_distress] risk out of bounds\n";
      ++fails;
    }
    if (p.repReward < 0.0 || p.repReward > 5.0) {
      std::cerr << "[test_distress] repReward out of bounds\n";
      ++fails;
    }
    if (p.rewardCr < 0.0 || p.rewardCr > 10000.0) {
      std::cerr << "[test_distress] rewardCr out of bounds\n";
      ++fails;
    }
    if (p.ambush && p.pirateCount <= 0) {
      std::cerr << "[test_distress] ambush requires pirates\n";
      ++fails;
    }
    if (!p.ambush && p.pirateCount != 0) {
      std::cerr << "[test_distress] non-ambush should not have pirates\n";
      ++fails;
    }
    if (!p.hasVictim) {
      if (!(p.scenario == DistressScenario::Ambush || p.needUnits <= 1e-9)) {
        std::cerr << "[test_distress] no-victim plan should not require cargo\n";
        ++fails;
      }
    } else {
      if (p.needUnits <= 1e-6) {
        std::cerr << "[test_distress] victim plan requires non-zero needUnits\n";
        ++fails;
      }
      if (p.rewardCr <= 1e-6) {
        std::cerr << "[test_distress] victim plan requires non-zero rewardCr\n";
        ++fails;
      }
    }
  }

  // Different seeds should generally change the plan.
  {
    const DistressPlan a = planDistressEncounter(seed, sys, sig, t, faction);
    const DistressPlan b = planDistressEncounter(seed + 1, sys, sig, t, faction);

    const bool same = (a.scenario == b.scenario && a.hasVictim == b.hasVictim && a.ambush == b.ambush && a.pirateCount == b.pirateCount &&
                       a.needCommodity == b.needCommodity && nearly(a.needUnits, b.needUnits) && nearly(a.rewardCr, b.rewardCr));
    if (same) {
      std::cerr << "[test_distress] changing universe seed should usually change output\n";
      ++fails;
    }
  }

  return fails;
}
