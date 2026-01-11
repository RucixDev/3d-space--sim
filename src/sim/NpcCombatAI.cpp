#include "stellar/sim/NpcCombatAI.h"

#include "stellar/core/Clamp.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

// Greedy pip allocator that respects kPipMax and always returns exactly kPipTotal.
static Pips allocatePips(double wEng, double wWep, double wSys, core::u64 seed) {
  // Safety: keep weights non-negative.
  wEng = std::max(0.0, wEng);
  wWep = std::max(0.0, wWep);
  wSys = std::max(0.0, wSys);

  // If everything is zero, fall back.
  if (wEng + wWep + wSys < 1e-9) {
    Pips p{2, 2, 2};
    normalizePips(p);
    return p;
  }

  Pips p{0, 0, 0};

  core::SplitMix64 rng(seed);

  auto score = [&](double w, int cur) {
    if (cur >= kPipMax) return -1e30;
    // Diminishing returns: prefer spreading pips unless weight is significantly larger.
    const double denom = 1.0 + 0.85 * (double)cur;
    // Tiny deterministic jitter to break ties.
    const double jitter = 1.0 + (rng.nextDouble() - 0.5) * 0.00025;
    return (w / denom) * jitter;
  };

  for (int i = 0; i < kPipTotal; ++i) {
    const double sE = score(wEng, p.eng);
    const double sW = score(wWep, p.wep);
    const double sS = score(wSys, p.sys);

    if (sE >= sW && sE >= sS) {
      ++p.eng;
    } else if (sW >= sS) {
      ++p.wep;
    } else {
      ++p.sys;
    }
  }

  normalizePips(p);
  return p;
}

NpcPipDecision decideNpcPips(const NpcPipContext& in, core::u64 seed) {
  NpcPipContext ctx = in;

  ctx.aiSkill = clamp01(ctx.aiSkill);
  ctx.shieldFrac = clamp01(ctx.shieldFrac);
  ctx.hullFrac = clamp01(ctx.hullFrac);
  ctx.engCapFrac = clamp01(ctx.engCapFrac);
  ctx.wepCapFrac = clamp01(ctx.wepCapFrac);
  ctx.sysCapFrac = clamp01(ctx.sysCapFrac);

  const double s = ctx.aiSkill;

  // Baselines by role.
  double wEng = 1.0;
  double wWep = 1.0;
  double wSys = 1.0;

  switch (ctx.role) {
    case NpcCombatRole::Pirate:
      // Aggressive, but not suicidal.
      wEng = 0.95;
      wWep = 1.18;
      wSys = 0.85;
      break;
    case NpcCombatRole::Police:
      // More durable and "steady".
      wEng = 0.90;
      wWep = 1.00;
      wSys = 1.18;
      break;
    case NpcCombatRole::Trader:
      // Wants to run and keep shields up.
      wEng = 1.30;
      wWep = 0.55;
      wSys = 1.05;
      break;
  }

  // Context adjustments.
  if (ctx.inCombat) {
    // Skilled pilots bias a bit more to WEP/ENG in combat.
    wWep += 0.30 + 0.35 * s;
    wEng += 0.10 + 0.15 * s;
  }

  if (ctx.wantsToFire) {
    wWep += 0.75 + 0.45 * s;
  }

  // If weapon capacitor is low, strongly bias to WEP to recover.
  if (ctx.inCombat) {
    const double low = 0.35;
    if (ctx.wepCapFrac < low) {
      const double t = (low - ctx.wepCapFrac) / std::max(1e-6, low);
      wWep += t * (0.90 + 0.55 * s);
    }
  }

  // Defensive pressure.
  if (ctx.underFire) {
    wSys += 0.65 + 0.45 * (1.0 - ctx.shieldFrac);
    wEng += 0.15; // keep maneuver authority
  }

  if (ctx.missileThreat) {
    // Missiles demand maneuvering + some shield sustain.
    wEng += 0.95;
    wSys += 0.65;
  }

  if (ctx.shieldFrac < 0.40) {
    const double t = (0.40 - ctx.shieldFrac) / 0.40;
    wSys += t * (1.05 + 0.25 * s);
  }

  if (ctx.hullFrac < 0.50) {
    // When hull is getting scary, prioritize escape over DPS.
    const double t = (0.50 - ctx.hullFrac) / 0.50;
    wEng += t * (0.55 + 0.35 * s);
    wWep *= (1.0 - 0.35 * t);
  }

  if (ctx.wantsToBoost) {
    wEng += 0.85 + 0.35 * s;
  }

  // If ENG capacitor is low, try to refill it.
  {
    const double low = 0.30;
    if (ctx.engCapFrac < low) {
      const double t = (low - ctx.engCapFrac) / std::max(1e-6, low);
      wEng += t * (0.85 + 0.25 * s);
    }
  }

  // If SYS capacitor is low while shields aren't full, bias to SYS.
  if ((ctx.shieldFrac < 0.95 || ctx.underFire) && ctx.sysCapFrac < 0.25) {
    const double t = (0.25 - ctx.sysCapFrac) / 0.25;
    wSys += t * 0.55;
  }

  // Fleeing overrides most other considerations.
  if (ctx.fleeing) {
    wEng += 3.0;
    wSys += 0.85;
    wWep *= 0.10;
  }

  NpcPipDecision out{};
  out.pips = allocatePips(wEng, wWep, wSys, seed);

  // Boost policy: allow in situations where it looks intentional, not constant.
  out.allowBoost = ctx.fleeing || ctx.missileThreat || (ctx.inCombat && (ctx.wantsToBoost || s > 0.60) && ctx.engCapFrac > 0.15);

  // Suggested orbit/jink intensity.
  if (ctx.fleeing) {
    double base = 0.10 + 0.18 * s;
    if (ctx.missileThreat) base += 0.10;
    out.orbitStrafe01 = std::clamp(base, 0.0, 0.60);
  } else if (ctx.inCombat) {
    double base = 0.08;
    if (ctx.role == NpcCombatRole::Pirate) base = 0.10;
    if (ctx.role == NpcCombatRole::Police) base = 0.09;

    base += 0.22 * s;

    if (ctx.underFire) base += 0.13;
    if (ctx.shieldFrac < 0.35) base += 0.10;
    if (ctx.missileThreat) base += 0.06;

    out.orbitStrafe01 = std::clamp(base, 0.0, 0.75);
  } else {
    out.orbitStrafe01 = 0.0;
  }

  return out;
}

} // namespace stellar::sim
