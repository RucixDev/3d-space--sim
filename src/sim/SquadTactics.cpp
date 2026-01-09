#include "stellar/sim/SquadTactics.h"

#include "stellar/core/Hash.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {
namespace {

static double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

// Convert a 64-bit hash to a stable double in [0, 1).
static double unitFromHash(core::u64 h) {
  // Use the top 53 bits as a mantissa to avoid precision issues.
  const core::u64 mant = (h >> 11) & ((core::u64(1) << 53) - 1);
  return (double)mant / (double)(core::u64(1) << 53);
}

// Small deterministic jitter to avoid every squad reacting on exactly the same frame.
static double moraleJitter(core::u64 universeSeed, core::u64 groupId, core::u64 salt) {
  const core::u64 h = core::hashCombine(core::hashCombine(universeSeed, groupId), salt);
  const double u = unitFromHash(h);
  // [-1, +1] scaled.
  return (u * 2.0 - 1.0) * 0.04;
}

} // namespace

PirateSquadDecision evaluatePirateSquad(const SquadSnapshot& snap,
                                        const PirateSquadContext& ctx,
                                        core::u64 universeSeed,
                                        core::u64 groupId) {
  PirateSquadDecision out{};

  if (snap.totalCount <= 0 || snap.aliveCount <= 0) {
    out.morale01 = 0.0;
    out.retreat = false;
    out.retreatDurationSec = 0.0;
    return out;
  }

  const double alive = (double)snap.aliveCount;
  const double total = (double)std::max(1, snap.totalCount);
  const double casualtyFrac = clamp01(1.0 - (alive / total));
  const double underFireFrac = clamp01((double)snap.underFireCount / std::max(1.0, alive));

  // Base morale: pirates are opportunistic and vary with system danger.
  double m = 0.52;

  // Bigger packs are bolder.
  if (snap.aliveCount >= 3) m += 0.18;
  else if (snap.aliveCount == 2) m += 0.10;

  // System knobs.
  m += 0.22 * (clamp01(ctx.piracy01) - 0.5);
  m -= 0.24 * (clamp01(ctx.security01) - 0.5);

  // More valuable target => more motivated.
  {
    const double cargoN = clamp01(std::log10(1.0 + std::max(0.0, ctx.playerCargoValueCr) / 450.0));
    m += 0.10 * cargoN;
  }

  if (ctx.committedRaid) m += 0.08;

  // Damage/casualties.
  m -= 0.55 * casualtyFrac;
  m -= 0.35 * (1.0 - clamp01(snap.avgHullFrac));
  m -= 0.12 * (1.0 - clamp01(snap.avgShieldFrac));

  // Losing the leader hurts morale.
  if (!snap.leaderAlive) m -= 0.28;

  // Sustained incoming fire makes opportunists bail faster.
  m -= 0.10 * underFireFrac;

  // Small deterministic noise.
  m += moraleJitter(universeSeed, groupId, 0xA11CE5u);

  m = clamp01(m);
  out.morale01 = m;

  // Retreat rules: pirates break off when morale is low, when leader is down,
  // or after heavy casualties.
  const bool leaderDown = !snap.leaderAlive;
  const bool heavyLosses = casualtyFrac > 0.60;

  out.retreat = (m < 0.18) || (leaderDown && m < 0.34) || (heavyLosses && m < 0.48);

  if (out.retreat) {
    // How long they stay away scales with how badly morale collapsed.
    const double sec = 70.0 + 180.0 * (1.0 - m) + 45.0 * casualtyFrac;
    out.retreatDurationSec = std::clamp(sec, 60.0, 300.0);
  }

  return out;
}

PoliceSquadDecision evaluatePoliceSquad(const SquadSnapshot& snap,
                                        const PoliceSquadContext& ctx,
                                        core::u64 universeSeed,
                                        core::u64 groupId) {
  PoliceSquadDecision out{};

  if (snap.totalCount <= 0 || snap.aliveCount <= 0) {
    out.morale01 = 0.0;
    out.retreat = false;
    out.retreatDurationSec = 0.0;
    out.callReinforcements = false;
    out.reinforcementUrgency01 = 0.0;
    return out;
  }

  const double alive = (double)snap.aliveCount;
  const double total = (double)std::max(1, snap.totalCount);
  const double casualtyFrac = clamp01(1.0 - (alive / total));
  const double underFireFrac = clamp01((double)snap.underFireCount / std::max(1.0, alive));
  const double heatN = clamp01(ctx.policeHeat / 6.0);
  const double bountyN = clamp01(std::max(0.0, ctx.playerBountyCr) / 4500.0);

  // Base morale: trained/security forces are more disciplined.
  double m = 0.78;

  // System security makes them bolder.
  m += 0.18 * (clamp01(ctx.security01) - 0.5);

  // Being on an active pursuit increases commitment.
  if (ctx.playerWanted) m += 0.06;

  // High heat + high bounty => more willingness to stick in.
  m += 0.06 * heatN;
  m += 0.08 * bountyN;

  // Damage/casualties.
  m -= 0.48 * casualtyFrac;
  m -= 0.28 * (1.0 - clamp01(snap.avgHullFrac));
  m -= 0.10 * (1.0 - clamp01(snap.avgShieldFrac));

  // Leader loss matters, but doesn't instantly rout them.
  if (!snap.leaderAlive) m -= 0.10;

  // Under sustained fire, they rely more on reinforcements.
  m -= 0.06 * underFireFrac;

  // Small deterministic noise.
  m += moraleJitter(universeSeed, groupId, 0xBADC0FFEu);

  m = clamp01(m);
  out.morale01 = m;

  // Police almost never "retreat" in the prototype; instead they should call for backup.
  out.retreat = (m < 0.08) && (casualtyFrac > 0.55) && ctx.fightingPlayer;
  if (out.retreat) {
    const double sec = 80.0 + 140.0 * (1.0 - m);
    out.retreatDurationSec = std::clamp(sec, 60.0, 240.0);
  }

  // Reinforcement logic.
  out.callReinforcements = false;
  out.reinforcementUrgency01 = 0.0;

  if (ctx.fightingPlayer && ctx.playerWanted) {
    // Urgency rises if morale is slipping, if casualties are high, if they're under fire,
    // and if the player has a meaningful bounty.
    double urg = 0.0;
    urg = std::max(urg, clamp01((0.55 - m) / 0.55));
    urg = std::max(urg, 0.75 * casualtyFrac);
    urg = std::max(urg, 0.55 * underFireFrac);
    urg = std::max(urg, 0.35 * bountyN);

    // System security + existing heat reduce urgency slightly (help avoid infinite escalation).
    urg *= (1.0 - 0.10 * clamp01(ctx.security01));
    urg *= (1.0 - 0.18 * heatN);

    urg = clamp01(urg);
    out.reinforcementUrgency01 = urg;
    out.callReinforcements = (urg > 0.32);
  }

  return out;
}

} // namespace stellar::sim
