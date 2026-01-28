#include "stellar/sim/GalNet.h"

#include "stellar/math/Math.h"
#include "stellar/sim/SystemEventEconomy.h"
#include "stellar/sim/Universe.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

namespace stellar::sim {

GalNetAutoBroadcastDecision galNetMaybeAutoBroadcast(GalNetAnnounceState& ioState,
                                                     const SystemEvent& ev,
                                                     double minSeverity01,
                                                     bool autoEnabled,
                                                     bool broadcastEventEnds) {
  GalNetAutoBroadcastDecision d{};

  // Treat NaNs as disabled.
  if (!std::isfinite(minSeverity01)) minSeverity01 = 1.0;
  minSeverity01 = std::clamp(minSeverity01, 0.0, 1.0);

  // Only evaluate once per cycle boundary.
  // SystemEvents guarantees startDay is stable within a cycle.
  if (ioState.lastCycleStartDay == ev.startDay) {
    return d;
  }

  const bool severeEnough = ev.active && (std::clamp(ev.severity01, 0.0, 1.0) + 1e-9 >= minSeverity01);

  if (autoEnabled && severeEnough) {
    d.shouldPublish = true;
    d.allowWhenNoEvent = false;
  } else if (autoEnabled && broadcastEventEnds && ioState.lastCycleHadActiveEvent && !ev.active) {
    d.shouldPublish = true;
    d.allowWhenNoEvent = true;
  }

  // Update tracking even if autoEnabled is off, to avoid retroactive spam when enabling.
  ioState.lastCycleStartDay = ev.startDay;
  ioState.lastCycleHadActiveEvent = ev.active;
  return d;
}

const char* galNetSystemEventKindLabel(SystemEventKind kind) {
  // Keep the labels stable and explicit, even if systemEventKindName() changes.
  switch (kind) {
    case SystemEventKind::TradeBoom:            return "Trade Boom";
    case SystemEventKind::TradeBust:            return "Trade Bust";
    case SystemEventKind::PirateRaid:           return "Pirate Raid";
    case SystemEventKind::CivilUnrest:          return "Civil Unrest";
    case SystemEventKind::SecurityCrackdown:    return "Security Crackdown";
    case SystemEventKind::ResearchBreakthrough: return "Research Breakthrough";
    case SystemEventKind::None:
    default:                                    return "None";
  }
}

const char* galNetSystemEventKindTip(SystemEventKind kind) {
  switch (kind) {
    case SystemEventKind::TradeBoom:
      return "Higher traffic and more lawful activity; good time for honest trade routes.";
    case SystemEventKind::TradeBust:
      return "Thin margins and jittery markets; expect fewer traders and rougher deals.";
    case SystemEventKind::PirateRaid:
      return "Pirates are active; keep your shields up and consider escorts.";
    case SystemEventKind::CivilUnrest:
      return "Unpredictable conditions; opportunists and agitators may be active.";
    case SystemEventKind::SecurityCrackdown:
      return "Increased patrol presence; illegal cargo and wanted pilots face more scrutiny.";
    case SystemEventKind::ResearchBreakthrough:
      return "Experimental traffic and unusual activity; keep an eye out for odd signals.";
    case SystemEventKind::None:
    default:
      return "";
  }
}

std::string factionNameOrId(core::u32 factionId, const std::vector<Faction>& factions) {
  if (factionId == 0) return "Independent";
  for (const auto& f : factions) {
    if (f.id == factionId) return f.name;
  }
  return "Faction " + std::to_string(factionId);
}

static int pct01(double x) {
  if (!std::isfinite(x)) return 0;
  return (int)std::lround(std::clamp(x, 0.0, 1.0) * 100.0);
}

static std::string fmtDeltaPct(int d) {
  if (d == 0) return "0";
  return (d > 0 ? "+" : "") + std::to_string(d);
}

static std::string fmtEndsIn(double remainingDays) {
  if (!std::isfinite(remainingDays)) return "n/a";
  remainingDays = std::max(0.0, remainingDays);
  const int remainingHours = (int)std::ceil(remainingDays * 24.0);
  if (remainingHours < 48) {
    return std::to_string(remainingHours) + "h";
  }
  return std::to_string((int)std::ceil(remainingDays)) + "d";
}

static CommsChannel channelForEvent(const SystemEvent& ev) {
  if (!ev.active) return CommsChannel::System;
  switch (ev.kind) {
    case SystemEventKind::TradeBoom:
    case SystemEventKind::TradeBust:
      return CommsChannel::Trade;
    case SystemEventKind::PirateRaid:
      return CommsChannel::Pirate;
    case SystemEventKind::SecurityCrackdown:
      return CommsChannel::Security;
    case SystemEventKind::CivilUnrest:
    case SystemEventKind::ResearchBreakthrough:
    case SystemEventKind::None:
    default:
      return CommsChannel::System;
  }
}

GalNetBulletinResult makeGalNetBulletin(const StarSystem& sys,
                                        const SystemConditionsSnapshot& snap,
                                        double timeDays,
                                        std::string_view controllingFactionName,
                                        std::string_view contextTag,
                                        bool allowWhenNoEvent) {
  GalNetBulletinResult r{};

  if (sys.stub.id == 0) {
    r.ok = false;
    r.reason = "Invalid system id.";
    return r;
  }

  if (!snap.event.active && !allowWhenNoEvent) {
    r.ok = false;
    r.reason = "No active event.";
    return r;
  }

  const int baseSec = pct01(snap.base.security01);
  const int effSec  = pct01(snap.effective.security01);
  const int basePir = pct01(snap.base.piracy01);
  const int effPir  = pct01(snap.effective.piracy01);
  const int baseTrf = pct01(snap.base.traffic01);
  const int effTrf  = pct01(snap.effective.traffic01);

  const int dSec = effSec - baseSec;
  const int dPir = effPir - basePir;
  const int dTrf = effTrf - baseTrf;

  const double remainingDays = (std::isfinite(snap.event.endDay) ? std::max(0.0, snap.event.endDay - timeDays)
                                                                 : 0.0);
  const std::string endsIn = fmtEndsIn(remainingDays);

  const char* eventLabel = galNetSystemEventKindLabel(snap.event.kind);
  const char* tip = galNetSystemEventKindTip(snap.event.kind);

  std::string subject;
  if (snap.event.active) {
    subject = std::string("GalNet: ") + eventLabel + " — " + sys.stub.name;
  } else {
    subject = std::string("GalNet: System status update — ") + sys.stub.name;
  }

  std::ostringstream body;
  if (!contextTag.empty()) {
    body << contextTag << "\n\n";
  }

  body << "System: " << sys.stub.name << "\n";

  if (snap.base.controllingFactionId != 0) {
    if (!controllingFactionName.empty()) {
      body << "Controlling faction: " << controllingFactionName << "\n";
    } else {
      body << "Controlling faction: Faction " << snap.base.controllingFactionId << "\n";
    }
  }

  if (snap.event.active) {
    const int sev = pct01(snap.event.severity01);
    body << "Event: " << eventLabel << "\n";
    body << "Severity: " << sev << "%\n";
    body << "Ends in: " << endsIn << "\n\n";
  } else {
    body << "Event: None\n";
    body << "Cycle ends in: " << endsIn << "\n\n";
  }

  body << "Conditions:\n";
  body << "  Security: " << effSec << "% (" << fmtDeltaPct(dSec) << "%)\n";
  body << "  Piracy: " << effPir << "% (" << fmtDeltaPct(dPir) << "%)\n";
  body << "  Traffic: " << effTrf << "% (" << fmtDeltaPct(dTrf) << "%)\n";

  {
    const std::string econLine = systemEventEconomySummary(snap.event);
    if (!econLine.empty()) {
      body << "\n" << econLine << "\n";
    }
  }

  if (snap.event.active && tip && tip[0] != 0) {
    body << "\nTip: " << tip << "\n";
  }

  GalNetBulletin b{};
  b.hasActiveEvent = snap.event.active;
  b.importance01 = b.hasActiveEvent ? std::clamp(snap.event.severity01, 0.0, 1.0) : 0.25;

  b.msg.timeDays = timeDays;
  b.msg.channel = channelForEvent(snap.event);
  b.msg.from = "GalNet";
  b.msg.subject = std::move(subject);
  b.msg.body = body.str();
  b.msg.systemId = sys.stub.id;
  b.msg.factionId = snap.base.controllingFactionId;

  r.ok = true;
  r.bulletin = std::move(b);
  return r;
}

GalNetBulletinResult makeGalNetBulletin(const Universe& universe,
                                        SystemId systemId,
                                        double timeDays,
                                        std::string_view contextTag,
                                        bool allowWhenNoEvent) {
  GalNetBulletinResult r{};
  if (systemId == 0) {
    r.ok = false;
    r.reason = "Invalid system id.";
    return r;
  }

  const StarSystem& sys = universe.getSystem(systemId);

  const SystemSecurityDeltaState* delta = nullptr;
  if (const auto* map = universe.systemSecurityDeltaMap()) {
    auto it = map->find(systemId);
    if (it != map->end()) delta = &it->second;
  }

  const SystemConditionsSnapshot snap = snapshotSystemConditions(
      universe.seed(),
      sys,
      timeDays,
      delta,
      universe.systemSecurityDynamicsParams(),
      universe.systemEventParams());

  const std::string factionName = factionNameOrId(snap.base.controllingFactionId, universe.factions());

  return makeGalNetBulletin(sys,
                            snap,
                            timeDays,
                            /*controllingFactionName=*/factionName,
                            contextTag,
                            allowWhenNoEvent);
}

} // namespace stellar::sim
