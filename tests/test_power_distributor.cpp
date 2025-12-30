#include "stellar/sim/PowerDistributor.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_power_distributor() {
  int fails = 0;

  // --- Pips normalization invariants ---
  {
    sim::Pips p{10, -2, 0};
    sim::normalizePips(p);

    const int sum = p.eng + p.wep + p.sys;
    if (sum != sim::kPipTotal) {
      std::cerr << "[test_power_distributor] normalizePips sum mismatch. got "
                << sum << " expected " << sim::kPipTotal << "\n";
      ++fails;
    }
    if (p.eng < 0 || p.eng > sim::kPipMax || p.wep < 0 || p.wep > sim::kPipMax || p.sys < 0 ||
        p.sys > sim::kPipMax) {
      std::cerr << "[test_power_distributor] normalizePips clamp failed. got eng=" << p.eng
                << " wep=" << p.wep << " sys=" << p.sys << "\n";
      ++fails;
    }
  }

  // --- Sum=0 resets to default 2/2/2 ---
  {
    sim::Pips p{0, 0, 0};
    sim::normalizePips(p);
    if (p.eng != 2 || p.wep != 2 || p.sys != 2) {
      std::cerr << "[test_power_distributor] normalizePips default failed. got eng=" << p.eng
                << " wep=" << p.wep << " sys=" << p.sys << "\n";
      ++fails;
    }
  }

  // --- Recharge redistributes from full channels ---
  {
    sim::DistributorConfig cfg{};
    cfg.capEng = 10.0;
    cfg.capWep = 10.0;
    cfg.capSys = 10.0;
    cfg.rechargePerSimSec = 30.0;

    sim::DistributorState st{};
    st.eng = 10.0;
    st.wep = 0.0;
    st.sys = 0.0;

    sim::Pips p{2, 2, 2};
    sim::stepDistributor(st, cfg, p, 1.0);

    if (!approx(st.eng, 10.0) || !approx(st.wep, 10.0) || !approx(st.sys, 10.0)) {
      std::cerr << "[test_power_distributor] stepDistributor redistribution failed. got eng="
                << st.eng << " wep=" << st.wep << " sys=" << st.sys << "\n";
      ++fails;
    }
  }

  // --- Partial boost consumption ---
  {
    sim::DistributorConfig cfg{};
    cfg.capEng = 100.0;
    cfg.boostCostPerSimSec = 10.0;

    sim::DistributorState st{};
    st.eng = 5.0;

    const double frac = sim::consumeBoostFraction(st, cfg, 1.0);
    if (!approx(frac, 0.5, 1e-9) || !approx(st.eng, 0.0, 1e-9)) {
      std::cerr << "[test_power_distributor] consumeBoostFraction failed. frac=" << frac
                << " eng=" << st.eng << "\n";
      ++fails;
    }
  }

  // --- Weapon capacitor cost heuristic (beam laser) ---
  {
    const auto& w = sim::weaponDef(sim::WeaponType::BeamLaser);
    const double cost = sim::weaponCapacitorCost(w);
    const double expected = 5.0 * w.cooldownSimSec + 0.10 * w.dmg;
    if (!approx(cost, expected, 1e-12)) {
      std::cerr << "[test_power_distributor] weaponCapacitorCost mismatch. got=" << cost
                << " expected=" << expected << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_power_distributor] PASS\n";
  }
  return fails;
}
