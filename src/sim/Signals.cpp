#include "stellar/sim/Signals.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"
#include "stellar/sim/SaveGame.h"
#include "stellar/sim/SecurityModel.h"
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
                                      const TrafficLedger* trafficLedger) {
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

  // System-stable key used by the renderer/game prototype as a base for procedural ids.
  // (We keep this scheme here so it can be shared when the game integrates this module.)
  const core::u64 sysKey = core::hashCombine(universeSeed, static_cast<core::u64>(system.stub.id));

  // Derive compact system-level security knobs used to shape derelict content.
  const auto sec = systemSecurityProfile(universeSeed, system);

  // --- Persistent resource fields ---
  if (params.resourceFieldCount > 0) {
    out.resourceFields = generateResourceFields(universeSeed, system.stub.id, anchorPosKm, anchorCommsKm, params.resourceFieldCount);
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
  if (params.includeDailyDerelict) {
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
                                       sec.piracy01,
                                       sec.security01,
                                       sec.contest01,
                                       /*missionSite=*/false,
                                       /*includeDayStamp=*/true);
    out.sites.push_back(s);
  }

  // --- Distress calls ---
  if (params.includeDistress && params.distressPerDay > 0) {
    const int count = std::clamp(params.distressPerDay, 0, 8);
    const double ttl = (params.distressTtlDays > 1e-6) ? params.distressTtlDays : 1.0;

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
                                       sec.piracy01,
                                       sec.security01,
                                       sec.contest01,
                                       /*missionSite=*/true,
                                       /*includeDayStamp=*/false);

    out.sites.push_back(s);
  }

  // --- Traffic convoy signals (moving lane traffic) ---
  if (params.includeTrafficConvoys) {
    std::vector<TrafficConvoyView> views;

    // Integration hook: if the caller provides a TrafficLedger recorded from
    // simulateNpcTradeTraffic(...), prefer turning those *actual shipments* into
    // convoy sites. If no shipments are available, fall back to the deterministic
    // lane prototype so the UI doesn't go empty.
    if (trafficLedger) {
      views = generateTrafficConvoysFromLedger(*trafficLedger,
                                              system,
                                              timeDays,
                                              params.trafficLaneParams.genWindowDays,
                                              params.trafficLaneParams.includeInactive,
                                              params.trafficLaneParams);
    }
    if (views.empty()) {
      views = generateTrafficConvoys(universeSeed, system, timeDays, params.trafficLaneParams);
    }

    out.sites.reserve(out.sites.size() + views.size());
    for (const auto& v : views) {
      // Skip already-arrived convoys even when includeInactive=true.
      if (v.convoy.arriveDay > 0.0 && timeDays > v.convoy.arriveDay) continue;

      // When includeInactive=false, generateTrafficConvoysFromLedger() already
      // filters, but keep this guard for safety (and for the fallback path).
      if (!params.trafficLaneParams.includeInactive && !v.state.active) continue;

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
