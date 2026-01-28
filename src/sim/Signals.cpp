#include "stellar/sim/Signals.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/SecurityModel.h"
#include "stellar/sim/SystemConditions.h"
#include "stellar/sim/TrafficConvoyLayer.h"
#include "stellar/sim/Units.h"
#include "stellar/sim/WorldIds.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

const char* signalKindName(SignalKind k) {
  switch (k) {
    case SignalKind::ResourceField: return "ResourceField";
    case SignalKind::Derelict: return "Derelict";
    case SignalKind::Distress: return "Distress";
    case SignalKind::MissionSalvage: return "MissionSalvage";
    case SignalKind::TrafficConvoy: return "TrafficConvoy";
    default: return "Unknown";
  }
}

bool isSignalResolved(const std::vector<core::u64>& resolvedSignalIds, core::u64 signalId) {
  for (const auto id : resolvedSignalIds) {
    if (id == signalId) return true;
  }
  return false;
}

static math::Vec3d randUnitVec(core::SplitMix64& rng) {
  // Rejection sampling inside unit sphere.
  for (int i = 0; i < 12; ++i) {
    math::Vec3d v{rng.range(-1.0, 1.0), rng.range(-1.0, 1.0), rng.range(-1.0, 1.0)};
    const double d2 = v.lengthSq();
    if (d2 > 1e-12 && d2 <= 1.0) return v * (1.0 / std::sqrt(d2));
  }
  return {1,0,0};
}

static const Station* pickAnchorStation(const StarSystem& system) {
  // Prefer stations that would plausibly support mining gameplay.
  for (const auto& st : system.stations) {
    if (st.type == econ::StationType::Mining || st.type == econ::StationType::Refinery) return &st;
  }
  if (!system.stations.empty()) return &system.stations.front();
  return nullptr;
}

SystemSignalPlan generateSystemSignals(core::u64 universeSeed,
                                      const StarSystem& system,
                                      double timeDays,
                                      const std::vector<Mission>& activeMissions,
                                      const std::vector<core::u64>& resolvedSignalIds,
                                      const SignalGenParams& params,
                                      const TrafficLedger* trafficLedger,
                                      const SystemConditionsSnapshot* conditions) {
  SystemSignalPlan out{};

  const Station* anchor = pickAnchorStation(system);
  if (!anchor) return out;

  // Integer day stamp (so non-traffic signals don't jitter on small dt).
  const core::u64 dayStamp = (core::u64)std::max(0.0, std::floor(timeDays));

  // Use a stable "anchor time" inside the current day so static sites don't
  // shift slightly when queried at different times (e.g. system entry vs. UI
  // refresh).
  const double anchorTimeDays = (double)dayStamp + 0.5;

  const math::Vec3d anchorPosKm = stationPosKm(*anchor, anchorTimeDays);
  const double anchorCommsKm = std::max(0.0, anchor->commsRangeKm);

  // Prefer aligning procedural belts/sheets to the anchor station's orbital plane.
  // (Cross(position, velocity) gives an angular-momentum direction.)
  math::Vec3d preferredPlaneN = math::cross(anchorPosKm, stationVelKmS(*anchor, anchorTimeDays));
  if (preferredPlaneN.lengthSq() > 1e-12) {
    preferredPlaneN = preferredPlaneN.normalized();
  } else {
    preferredPlaneN = {0.0, 0.0, 0.0};
  }

  // System-stable key used by the renderer/game prototype as a base for procedural ids.
  // (We keep this scheme here so it can be shared when the game integrates this module.)
  const core::u64 sysKey = core::hashCombine(universeSeed, static_cast<core::u64>(system.stub.id));

  // Derive compact system-level security knobs used to shape derelict content.
  const auto sec = systemSecurityProfile(universeSeed, system);

  // Optional: caller-provided conditions snapshot for event-reactive signal generation.
  //
  // When provided, we can bias signal density (distress / derelicts / convoys) based on the
  // deterministic system event layer and pass the *effective* security knobs into encounter
  // planners. This makes system events feel "real" by manifesting as on-the-ground signals.
  const SystemConditionsSnapshot* cond = nullptr;
  if (conditions && conditions->systemId == system.stub.id) {
    cond = conditions;
  }

  const SystemEvent* ev = nullptr;
  double sev01 = 0.0;

  // Default to deterministic baseline security for stability. When conditions are provided,
  // use the effective values (baseline + dynamics + event) to shape encounter plans.
  double piracy01 = sec.piracy01;
  double security01 = sec.security01;
  double contest01 = sec.contest01;
  if (cond) {
    piracy01 = cond->effective.piracy01;
    security01 = cond->effective.security01;
    contest01 = cond->effective.contest01;
    ev = &cond->event;
    if (ev->active) {
      sev01 = std::clamp(ev->severity01, 0.0, 1.0);
    }
  }

  // Apply event-reactive modifiers to the caller's generation knobs.
  SignalGenParams p = params;
  int extraDerelictsPerDay = 0;
  if (ev && ev->active) {
    switch (ev->kind) {
      case SystemEventKind::PirateRaid: {
        // Pirate raids generate lots of distress chatter and wreckage.
        if (p.includeDistress) {
          const int bonus = 1 + (int)std::llround(sev01 * 3.0); // +1..+4
          p.distressPerDay = std::clamp(p.distressPerDay + bonus, 0, 8);
        }
        if (p.includeDailyDerelict) {
          extraDerelictsPerDay = std::clamp((int)std::ceil(sev01 * 2.0), 0, 3);
        }
      } break;
      case SystemEventKind::CivilUnrest: {
        // Civil unrest also produces distress calls, but usually less than a raid.
        if (p.includeDistress) {
          const int bonus = (int)std::llround(sev01 * 2.0); // +0..+2
          p.distressPerDay = std::clamp(p.distressPerDay + bonus, 0, 8);
        }
        if (p.includeDailyDerelict) {
          extraDerelictsPerDay = std::clamp((int)std::ceil(sev01 * 1.0), 0, 2);
        }
      } break;
      case SystemEventKind::ResearchBreakthrough: {
        // Breakthroughs can leave behind experimental debris / unmanned probes.
        if (p.includeDailyDerelict) {
          extraDerelictsPerDay = std::clamp((int)std::ceil(sev01 * 1.2), 0, 2);
        }
      } break;
      case SystemEventKind::TradeBoom: {
        if (p.includeTrafficConvoys) {
          const int bonus = (int)std::llround(sev01 * 3.0);
          p.trafficLaneParams.convoysPerDayBase = std::max(0, p.trafficLaneParams.convoysPerDayBase + bonus);
          p.trafficLaneParams.maxConvoysPerDay = std::max(p.trafficLaneParams.maxConvoysPerDay, p.trafficLaneParams.convoysPerDayBase * 4);
        }
      } break;
      case SystemEventKind::TradeBust: {
        if (p.includeTrafficConvoys) {
          const int cut = (int)std::llround(sev01 * 2.0);
          p.trafficLaneParams.convoysPerDayBase = std::max(0, p.trafficLaneParams.convoysPerDayBase - cut);
        }
      } break;
      case SystemEventKind::SecurityCrackdown: {
        // Crackdowns tend to reduce opportunistic distress spam and convoy piracy.
        if (p.includeDistress) {
          const int cut = (int)std::llround(sev01 * 1.0);
          p.distressPerDay = std::max(0, p.distressPerDay - cut);
        }
      } break;
      default: break;
    }
  }

  // --- Persistent resource fields ---
  if (p.resourceFieldCount > 0) {
    out.resourceFields = generateResourceFields(universeSeed, system.stub.id, anchorPosKm, anchorCommsKm, p.resourceFieldCount, preferredPlaneN);
    out.sites.reserve(out.sites.size() + out.resourceFields.fields.size());
    for (const auto& f : out.resourceFields.fields) {
      SignalSite s{};
      s.id = f.id;
      s.kind = SignalKind::ResourceField;
      s.posKm = f.posKm;
      s.expireDay = 0.0;
      s.resolved = false;
      s.fieldKind = f.kind;
      out.sites.push_back(s);
    }
  }

  // --- Daily derelict salvage site ---
  if (p.includeDailyDerelict) {
    // typeCode=2 in the current game prototype for "daily derelict".
    const core::u64 id = makeDeterministicWorldId(core::hashCombine(sysKey, 2ull), dayStamp);

    core::SplitMix64 srng(core::hashCombine(id, 0xD311E1C7ull));
    const math::Vec3d dir = randUnitVec(srng);

    SignalSite s{};
    s.id = id;
    s.kind = SignalKind::Derelict;
    s.posKm = anchorPosKm + dir * (anchorCommsKm * 1.6 + 190000.0);
    s.expireDay = (double)dayStamp + 1.0;
    s.resolved = isSignalResolved(resolvedSignalIds, id);

    s.hasDerelictPlan = true;
    s.derelict = planDerelictEncounter(universeSeed,
                                       system.stub.id,
                                       id,
                                       anchorTimeDays,
                                       piracy01,
                                       security01,
                                       contest01,
                                       /*missionSite=*/false,
                                       /*includeDayStamp=*/true);
    out.sites.push_back(s);
  }


  // --- Event-reactive derelict "aftermath" sites ---
  //
  // These are extra short-lived salvage signals that manifest during certain deterministic
  // system events (e.g. pirate raids). They make events feel tangible: you can *see*
  // the chaos in normal space, not just in UI numbers.
  if (extraDerelictsPerDay > 0 && ev && p.includeDailyDerelict) {
    const int count = std::clamp(extraDerelictsPerDay, 0, 6);

    // typeCode=6 reserved for event-reactive derelict sites. Mix in event kind so a new
    // event within the same day doesn't reuse ids (avoids weird "already resolved" states).
    const core::u64 baseKey = core::hashCombine(sysKey, 6ull);
    const core::u64 kindKey = core::hashCombine(baseKey, (core::u64)ev->kind);

    for (int i = 0; i < count; ++i) {
      const core::u64 salt = dayStamp * 16ull + (core::u64)i;
      const core::u64 id = makeDeterministicWorldId(kindKey, salt);

      core::SplitMix64 srng(core::hashCombine(id, 0xE6E1D11Cull));
      const math::Vec3d dir = randUnitVec(srng);
      const double distKm = anchorCommsKm * 1.55 + 210000.0 + srng.range(0.0, 140000.0);

      SignalSite s{};
      s.id = id;
      s.kind = SignalKind::Derelict;
      s.posKm = anchorPosKm + dir * distKm;
      s.expireDay = (double)dayStamp + 1.0;
      s.resolved = isSignalResolved(resolvedSignalIds, id);

      s.hasDerelictPlan = true;
      s.derelict = planDerelictEncounter(universeSeed,
                                         system.stub.id,
                                         id,
                                         anchorTimeDays,
                                         piracy01,
                                         security01,
                                         contest01,
                                         /*missionSite=*/false,
                                         /*includeDayStamp=*/true);

      out.sites.push_back(s);
    }
  }
  // --- Distress calls ---
  if (p.includeDistress && p.distressPerDay > 0) {
    const int count = std::clamp(p.distressPerDay, 0, 8);
    const double ttl = (p.distressTtlDays > 1e-6) ? p.distressTtlDays : 1.0;

    // Use the station faction if available; fall back to system faction.
    const core::u32 localFactionId = (anchor->factionId != 0) ? anchor->factionId : system.stub.factionId;

    for (int i = 0; i < count; ++i) {
      // typeCode=3 reserved for distress calls.
      const core::u64 salt = dayStamp * 16ull + (core::u64)i;
      const core::u64 id = makeDeterministicWorldId(core::hashCombine(sysKey, 3ull), salt);

      core::SplitMix64 srng(core::hashCombine(id, 0xD15A2E55ull));
      const math::Vec3d dir = randUnitVec(srng);
      const double distKm = anchorCommsKm * 1.35 + 160000.0 + srng.range(0.0, 110000.0);

      SignalSite s{};
      s.id = id;
      s.kind = SignalKind::Distress;
      s.posKm = anchorPosKm + dir * distKm;
      // Anchor TTL to the integer day stamp so repeated queries within a day
      // can't extend the lifetime of a deterministic distress call.
      s.expireDay = (double)dayStamp + ttl;
      s.resolved = isSignalResolved(resolvedSignalIds, id);

      s.hasDistressPlan = true;
      // Use the stable anchor time inside the day for determinism.
      s.distress = planDistressEncounter(universeSeed, system.stub.id, id, anchorTimeDays, localFactionId);

      out.sites.push_back(s);
    }
  }

  // --- Mission salvage sites (from active missions) ---
  // Mirrors the current renderer's behavior:
  //   - Use mission.targetNpcId as the stable signal id.
  //   - Place the site near the destination station.
  for (const auto& m : activeMissions) {
    if (m.type != MissionType::Salvage) continue;
    if (m.completed || m.failed) continue;
    if (m.toSystem != system.stub.id) continue;
    if (m.scanned) continue;
    if (m.targetNpcId == 0) continue;

    // Anchor near the destination station when available.
    const Station* baseSt = anchor;
    for (const auto& st : system.stations) {
      if (st.id == m.toStation) { baseSt = &st; break; }
    }

    const math::Vec3d basePos = stationPosKm(*baseSt, anchorTimeDays);
    const double commsKm = std::max(0.0, baseSt->commsRangeKm);

    core::SplitMix64 srng(core::hashCombine((core::u64)system.stub.id, (core::u64)m.targetNpcId));
    const math::Vec3d dir = randUnitVec(srng);
    const double distKm = commsKm * 1.8 + srng.range(150000.0, 260000.0);

    SignalSite s{};
    s.id = m.targetNpcId;
    s.kind = SignalKind::MissionSalvage;
    s.posKm = basePos + dir * distKm;

    // TTL is heuristic; keep it long enough for the player to travel.
    const double ttl = std::clamp(m.deadlineDay - timeDays + 1.0, 0.5, 6.0);
    s.expireDay = timeDays + ttl;
    s.resolved = false;
    s.missionId = m.id;

    // Mission salvage sites should keep their content stable across time (don't mix the day).
    s.hasDerelictPlan = true;
    s.derelict = planDerelictEncounter(universeSeed,
                                       system.stub.id,
                                       s.id,
                                       anchorTimeDays,
                                       piracy01,
                                       security01,
                                       sec.contest01,
                                       /*missionSite=*/true,
                                       /*includeDayStamp=*/false);

    out.sites.push_back(s);
  }

  // --- Traffic convoy signals (moving lane traffic) ---
  if (p.includeTrafficConvoys) {
    std::vector<TrafficConvoyView> views;

    // Integration hook: if the caller provides a TrafficLedger recorded from
    // simulateNpcTradeTraffic(...), prefer turning those *actual shipments* into
    // convoy sites. If no shipments are available, fall back to the deterministic
    // lane prototype so the UI doesn't go empty.
    if (trafficLedger) {
      views = generateTrafficConvoysFromLedger(*trafficLedger,
                                              system,
                                              timeDays,
                                              p.trafficLaneParams.genWindowDays,
                                              p.trafficLaneParams.includeInactive,
                                              p.trafficLaneParams);
    }
    if (views.empty()) {
      views = generateTrafficConvoys(universeSeed, system, timeDays, p.trafficLaneParams);
    }

    out.sites.reserve(out.sites.size() + views.size());
    for (const auto& v : views) {
      // Skip already-arrived convoys even when includeInactive=true.
      if (v.convoy.arriveDay > 0.0 && timeDays > v.convoy.arriveDay) continue;

      // When includeInactive=false, generateTrafficConvoysFromLedger() already
      // filters, but keep this guard for safety (and for the fallback path).
      if (!p.trafficLaneParams.includeInactive && !v.state.active) continue;

      SignalSite s{};
      s.id = v.convoy.id;
      s.kind = SignalKind::TrafficConvoy;
      s.posKm = v.state.posKm;
      s.expireDay = v.convoy.arriveDay;
      s.resolved = false;

      s.hasTrafficConvoy = true;
      s.trafficConvoy = v.convoy;
      s.trafficState = v.state;

      out.sites.push_back(s);
    }
  }

  // Stable ordering makes consumer code simpler.
  std::sort(out.sites.begin(), out.sites.end(), [](const SignalSite& a, const SignalSite& b) {
    if (a.kind != b.kind) return (core::u8)a.kind < (core::u8)b.kind;
    if (a.kind == SignalKind::TrafficConvoy && a.hasTrafficConvoy && b.hasTrafficConvoy) {
      if (a.trafficConvoy.departDay != b.trafficConvoy.departDay) {
        return a.trafficConvoy.departDay < b.trafficConvoy.departDay;
      }
    }
    return a.id < b.id;
  });

  return out;
}

} // namespace stellar::sim
