#include "stellar/sim/SystemEvents.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double clamp01(double v) { return std::clamp(v, 0.0, 1.0); }

const char* systemEventKindName(SystemEventKind k) {
  switch (k) {
    case SystemEventKind::TradeBoom: return "Trade Boom";
    case SystemEventKind::TradeBust: return "Trade Bust";
    case SystemEventKind::PirateRaid: return "Pirate Raid";
    case SystemEventKind::SecurityCrackdown: return "Security Crackdown";
    case SystemEventKind::CivilUnrest: return "Civil Unrest";
    case SystemEventKind::ResearchBreakthrough: return "Research Breakthrough";
    default: return "None";
  }
}

static double safeCycleDays(double d) {
  if (!std::isfinite(d) || d < 0.25) return 6.0;
  return d;
}

static double safeChance(double p) {
  if (!std::isfinite(p)) return 0.0;
  return std::clamp(p, 0.0, 1.0);
}

static double safeRange(double v, double fallback) {
  if (!std::isfinite(v)) return fallback;
  return v;
}

static double bias01(double x01) {
  return std::clamp(x01, 0.0, 1.0) - 0.5; // [-0.5, +0.5]
}

static double clampWeight(double w) {
  if (!std::isfinite(w)) return 0.0;
  return std::max(0.0, w);
}

SystemEvent generateSystemEvent(core::u64 universeSeed,
                                SystemId systemId,
                                double timeDays,
                                const SystemSecurityProfile& baseProfile,
                                const SystemEventParams& params) {
  SystemEvent ev{};
  ev.systemId = systemId;

  if (systemId == 0) return ev;

  const double cycleDays = safeCycleDays(params.cycleDays);
  const double tDays = std::max(0.0, safeRange(timeDays, 0.0));
  const core::i64 cycleIndex = (core::i64)std::floor(tDays / cycleDays);

  ev.startDay = (double)cycleIndex * cycleDays;
  ev.endDay = ev.startDay + cycleDays;

  core::u64 seed = universeSeed;
  seed = core::hashCombine(seed, (core::u64)systemId);
  seed = core::hashCombine(seed, (core::u64)cycleIndex);
  seed = core::hashCombine(seed, core::fnv1a64("system_events"));

  core::SplitMix64 rng(seed);

  // First gate: chance that this cycle has any event at all.
  if (!rng.chance(safeChance(params.eventChance))) {
    return ev;
  }

  // Base weights (sum doesn't matter; roulette selection normalizes).
  double wNone = 0.14;
  double wBoom = 0.16;
  double wBust = 0.14;
  double wRaid = 0.16;
  double wCrack = 0.14;
  double wUnrest = 0.14;
  double wResearch = 0.12;

  const double bSec = bias01(baseProfile.security01);
  const double bPir = bias01(baseProfile.piracy01);
  const double bTrf = bias01(baseProfile.traffic01);
  const double bCon = bias01(baseProfile.contest01);
  const auto& tr = baseProfile.traits;

  // Bias weights by the system's current conditions.
  wRaid *= (1.0 + 2.1 * bPir - 1.0 * bSec + 0.8 * std::max(0.0, bCon));
  wCrack *= (1.0 - 1.6 * bPir + 0.7 * bSec + 1.0 * std::max(0.0, tr.authority - 0.5));
  wBoom *= (1.0 + 1.8 * bTrf + 1.0 * std::max(0.0, tr.wealth - 0.5) + 0.35 * (tr.tech - 0.5));
  wBust *= (1.0 - 1.8 * bTrf + 1.0 * std::max(0.0, 0.5 - tr.wealth) + 0.40 * (0.5 - bSec));
  wUnrest *= (1.0 + 1.2 * std::max(0.0, 0.5 - tr.stability)
              + 0.9 * std::max(0.0, 0.5 - baseProfile.security01)
              + 0.7 * std::max(0.0, bCon)
              + 0.6 * std::max(0.0, tr.corruption - 0.5));
  wResearch *= (1.0 + 1.6 * std::max(0.0, tr.tech - 0.5)
                + 0.7 * std::max(0.0, tr.wealth - 0.5)
                - 0.9 * std::max(0.0, bPir));

  // If we don't have a controlling faction, de-emphasize economic/social states.
  if (baseProfile.controllingFactionId == 0) {
    wBoom *= 0.35;
    wBust *= 0.35;
    wCrack *= 0.55;
    wUnrest *= 0.55;
    wResearch *= 0.35;
  }

  // Sanitization.
  wNone = clampWeight(wNone);
  wBoom = clampWeight(wBoom);
  wBust = clampWeight(wBust);
  wRaid = clampWeight(wRaid);
  wCrack = clampWeight(wCrack);
  wUnrest = clampWeight(wUnrest);
  wResearch = clampWeight(wResearch);

  const double sum = wNone + wBoom + wBust + wRaid + wCrack + wUnrest + wResearch;
  if (sum <= 1e-12) {
    return ev;
  }

  const double r = rng.nextDouble() * sum;
  double t = 0.0;

  auto choose = [&](double w, SystemEventKind k) -> bool {
    t += w;
    if (r <= t) {
      ev.kind = k;
      return true;
    }
    return false;
  };

  if (!choose(wNone, SystemEventKind::None) &&
      !choose(wBoom, SystemEventKind::TradeBoom) &&
      !choose(wBust, SystemEventKind::TradeBust) &&
      !choose(wRaid, SystemEventKind::PirateRaid) &&
      !choose(wCrack, SystemEventKind::SecurityCrackdown) &&
      !choose(wUnrest, SystemEventKind::CivilUnrest) &&
      !choose(wResearch, SystemEventKind::ResearchBreakthrough)) {
    ev.kind = SystemEventKind::None;
  }

  if (ev.kind == SystemEventKind::None) {
    // Event slot exists this cycle, but it rolled "None".
    return ev;
  }

  const double minS = std::clamp(safeRange(params.minSeverity, 0.35), 0.0, 1.0);
  const double maxS = std::clamp(safeRange(params.maxSeverity, 1.00), 0.0, 1.0);
  const double sev = (maxS >= minS) ? rng.range(minS, maxS) : rng.range(maxS, minS);

  ev.active = true;
  ev.severity01 = sev;

  // Channel deltas (scaled by severity). Keep these fairly small so the layer
  // feels like "weather" on top of baseline + player-driven state.
  switch (ev.kind) {
    case SystemEventKind::TradeBoom:
      ev.securityDelta = +0.05 * sev;
      ev.piracyDelta = -0.06 * sev;
      ev.trafficDelta = +0.16 * sev;
      break;

    case SystemEventKind::TradeBust:
      ev.securityDelta = -0.05 * sev;
      ev.piracyDelta = +0.06 * sev;
      ev.trafficDelta = -0.16 * sev;
      break;

    case SystemEventKind::PirateRaid:
      ev.securityDelta = -0.08 * sev;
      ev.piracyDelta = +0.19 * sev;
      ev.trafficDelta = -0.11 * sev;
      break;

    case SystemEventKind::SecurityCrackdown:
      ev.securityDelta = +0.19 * sev;
      ev.piracyDelta = -0.17 * sev;
      ev.trafficDelta = -0.05 * sev;
      break;

    case SystemEventKind::CivilUnrest:
      ev.securityDelta = -0.16 * sev;
      ev.piracyDelta = +0.11 * sev;
      ev.trafficDelta = -0.09 * sev;
      break;

    case SystemEventKind::ResearchBreakthrough:
      ev.securityDelta = +0.05 * sev;
      ev.piracyDelta = -0.03 * sev;
      ev.trafficDelta = +0.07 * sev;
      break;

    default:
      break;
  }

  const double maxAbs = std::max(0.0, safeRange(params.maxAbsDelta, 0.22));
  ev.securityDelta = std::clamp(ev.securityDelta, -maxAbs, maxAbs);
  ev.piracyDelta = std::clamp(ev.piracyDelta, -maxAbs, maxAbs);
  ev.trafficDelta = std::clamp(ev.trafficDelta, -maxAbs, maxAbs);

  return ev;
}

SystemSecurityProfile applySystemEventToProfile(const SystemSecurityProfile& base,
                                                const SystemEvent& ev) {
  if (!ev.active) return base;

  SystemSecurityProfile out = base;
  out.security01 = clamp01(out.security01 + ev.securityDelta);
  out.piracy01 = clamp01(out.piracy01 + ev.piracyDelta);
  out.traffic01 = clamp01(out.traffic01 + ev.trafficDelta);
  return out;
}

} // namespace stellar::sim
