#include "stellar/sim/PowerDistributor.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

namespace {

struct BaseTuning {
  double capEng;
  double capWep;
  double capSys;
  double rechargePerSimSec;
  double boostCostPerSimSec;
  double shieldRegenCostPerPoint;
};

// Baseline tuning per hull. Distributor Mk scaling (kDistributors[mk].mult)
// applies to capacity + recharge, and improves efficiency (boost + shield regen
// cost) at higher Mk.
static constexpr BaseTuning kBaseByHull[] = {
  // ENG cap, WEP cap, SYS cap, recharge, boost cost, shield cost
  {  90.0,   90.0,  100.0,   24.0,     15.0,      0.60}, // Scout
  {  95.0,   85.0,  110.0,   24.0,     16.0,      0.65}, // Hauler
  {  90.0,  105.0,   95.0,   24.0,     14.0,      0.60}, // Fighter
};

int hullIndex(ShipHullClass hull) {
  const int idx = static_cast<int>(hull);
  return std::clamp(idx, 0, 2);
}

double weightForPips(int p) {
  const int pp = std::clamp(p, 0, kPipMax);
  // Non-linear weighting: concentrating pips gives a strong bias.
  // p=0 -> 1.0, p=4 -> 1 + 0.7*16 = 12.2
  return 1.0 + 0.7 * double(pp * pp);
}

} // namespace

void normalizePips(Pips& p) {
  auto clampP = [](int v) { return std::clamp(v, 0, kPipMax); };
  p.eng = clampP(p.eng);
  p.wep = clampP(p.wep);
  p.sys = clampP(p.sys);

  int sum = p.eng + p.wep + p.sys;
  if (sum == 0) {
    p.eng = 2;
    p.wep = 2;
    p.sys = 2;
    return;
  }

  // Reduce the largest channel(s) until we hit the total.
  while (sum > kPipTotal) {
    if (p.eng >= p.wep && p.eng >= p.sys && p.eng > 0) {
      --p.eng;
    } else if (p.wep >= p.eng && p.wep >= p.sys && p.wep > 0) {
      --p.wep;
    } else if (p.sys > 0) {
      --p.sys;
    } else {
      break;
    }
    --sum;
  }

  // Increase the smallest channel(s) until we hit the total.
  while (sum < kPipTotal) {
    if (p.eng <= p.wep && p.eng <= p.sys && p.eng < kPipMax) {
      ++p.eng;
    } else if (p.wep <= p.eng && p.wep <= p.sys && p.wep < kPipMax) {
      ++p.wep;
    } else if (p.sys < kPipMax) {
      ++p.sys;
    } else {
      break;
    }
    ++sum;
  }
}

DistributorConfig distributorConfig(ShipHullClass hull, int distributorMk) {
  const int h = hullIndex(hull);
  const int mk = std::clamp(distributorMk, 1, 3);
  const double mult = kDistributors[mk].mult;

  const BaseTuning base = kBaseByHull[h];

  DistributorConfig cfg{};
  cfg.capEng = base.capEng * mult;
  cfg.capWep = base.capWep * mult;
  cfg.capSys = base.capSys * mult;

  cfg.rechargePerSimSec = base.rechargePerSimSec * mult;

  // Efficiency gains at higher Mk: same output for less capacitor drain.
  cfg.boostCostPerSimSec = base.boostCostPerSimSec / mult;
  cfg.shieldRegenCostPerPoint = base.shieldRegenCostPerPoint / mult;

  return cfg;
}

DistributorState makeFull(const DistributorConfig& cfg) {
  DistributorState st{};
  st.eng = std::max(0.0, cfg.capEng);
  st.wep = std::max(0.0, cfg.capWep);
  st.sys = std::max(0.0, cfg.capSys);
  return st;
}

double consumeBoostFraction(DistributorState& st, const DistributorConfig& cfg, double dtSim) {
  if (dtSim <= 0.0) return 0.0;

  const double costPerSec = std::max(0.0, cfg.boostCostPerSimSec);
  const double need = costPerSec * dtSim;
  if (need <= 1e-12) return 0.0;

  const double cap = std::max(0.0, cfg.capEng);
  const double have = std::clamp(st.eng, 0.0, cap);

  const double frac = std::clamp(have / need, 0.0, 1.0);

  const double used = need * frac;
  st.eng = std::max(0.0, have - used);
  return frac;
}

void stepDistributor(DistributorState& st,
                     const DistributorConfig& cfg,
                     const Pips& pipsIn,
                     double dtSim) {
  if (dtSim <= 0.0) return;

  Pips p = pipsIn;
  normalizePips(p);

  const double capEng = std::max(0.0, cfg.capEng);
  const double capWep = std::max(0.0, cfg.capWep);
  const double capSys = std::max(0.0, cfg.capSys);

  st.eng = std::clamp(st.eng, 0.0, capEng);
  st.wep = std::clamp(st.wep, 0.0, capWep);
  st.sys = std::clamp(st.sys, 0.0, capSys);

  double remaining = std::max(0.0, cfg.rechargePerSimSec) * dtSim;
  if (remaining <= 1e-12) return;

  // Iteratively distribute recharge, redistributing shares that would overflow full caps.
  for (int iter = 0; iter < 6 && remaining > 1e-10; ++iter) {
    const bool openEng = st.eng + 1e-12 < capEng;
    const bool openWep = st.wep + 1e-12 < capWep;
    const bool openSys = st.sys + 1e-12 < capSys;

    const double wEng = openEng ? weightForPips(p.eng) : 0.0;
    const double wWep = openWep ? weightForPips(p.wep) : 0.0;
    const double wSys = openSys ? weightForPips(p.sys) : 0.0;

    const double wSum = wEng + wWep + wSys;
    if (wSum <= 0.0) break;

    double used = 0.0;

    auto add = [&](double& v, double cap, double w) {
      if (w <= 0.0) return;
      const double share = remaining * (w / wSum);
      const double space = cap - v;
      const double delta = std::min(space, share);
      v += delta;
      used += delta;
    };

    add(st.eng, capEng, wEng);
    add(st.wep, capWep, wWep);
    add(st.sys, capSys, wSys);

    remaining -= used;
    if (used <= 1e-12) break;
  }

  st.eng = std::clamp(st.eng, 0.0, capEng);
  st.wep = std::clamp(st.wep, 0.0, capWep);
  st.sys = std::clamp(st.sys, 0.0, capSys);
}

double shieldRegenMultiplierFromPips(int sysPips) {
  const int p = std::clamp(sysPips, 0, kPipMax);
  return 0.6 + 0.2 * double(p);
}

double weaponCapacitorCost(const WeaponDef& w) {
  // A compact heuristic that plays well with the default Mk1 recharge tuning:
  //  - rapid weapons get a low per-shot cost
  //  - heavy weapons pay via high per-shot damage
  const double baseDrain = 5.0;   // energy/sec baseline
  const double dmgCoeff = 0.10;   // energy per damage

  const double cd = std::max(0.0, w.cooldownSimSec);
  const double dmg = std::max(0.0, w.dmg);
  return std::max(0.0, baseDrain * cd + dmgCoeff * dmg);
}

} // namespace stellar::sim
