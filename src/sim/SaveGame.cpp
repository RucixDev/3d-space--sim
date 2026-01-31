#include "stellar/sim/SaveGame.h"

#include "stellar/core/Clamp.h"
#include "stellar/core/Base64.h"
#include "stellar/core/AtomicWriteFile.h"
#include "stellar/core/Log.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>

namespace stellar::sim {

bool saveToFile(const SaveGame& s, const std::string& path) {
  const std::filesystem::path dstPath(path);

  // Best-effort safety net: keep the previous successful save around as "*.bak".
  {
    std::error_code ec;
    if (std::filesystem::exists(dstPath, ec) && !ec) {
      auto bakPath = dstPath;
      bakPath += ".bak";
      std::filesystem::copy_file(dstPath, bakPath, std::filesystem::copy_options::overwrite_existing, ec);
      if (ec) {
        stellar::core::log(stellar::core::LogLevel::Warn,
                           "SaveGame: failed to write backup '" + bakPath.string() + "': " + ec.message());
      }
    }
  }

  // Serialize to memory first, then atomically replace the destination file.
  // This avoids corrupting the user's save on crashes or partial writes.
  std::ostringstream f;

  f.setf(std::ios::fixed);
  f.precision(8);

  f << "StellarForgeSave " << s.version << "\n";
  f << "seed " << s.seed << "\n";
  f << "timeDays " << s.timeDays << "\n";
  f << "currentSystem " << s.currentSystem << "\n";
  f << "dockedStation " << s.dockedStation << "\n";
  f << "lastDockedSystem " << s.lastDockedSystem << "\n";
  f << "lastDockedStation " << s.lastDockedStation << "\n";

  f << "shipPosKm " << s.shipPosKm.x << " " << s.shipPosKm.y << " " << s.shipPosKm.z << "\n";
  f << "shipVelKmS " << s.shipVelKmS.x << " " << s.shipVelKmS.y << " " << s.shipVelKmS.z << "\n";
  f << "shipOrient " << s.shipOrient.w << " " << s.shipOrient.x << " " << s.shipOrient.y << " " << s.shipOrient.z << "\n";
  f << "shipAngVel " << s.shipAngVelRadS.x << " " << s.shipAngVelRadS.y << " " << s.shipAngVelRadS.z << "\n";

  f << "credits " << s.credits << "\n";
  f << "insuranceDebtCr " << s.insuranceDebtCr << "\n";

  // Persistent progression unlock state (license-style gating).
  f << "creditsHighWaterCr " << s.creditsHighWaterCr << "\n";
  f << "repHighWater " << s.repHighWater << "\n";
  f << "unlockedHullMask " << s.unlockedHullMask << "\n";
  f << "unlockedCoreMk " << (int)s.unlockedCoreMk << "\n";
  f << "unlockedSmuggleHoldMk " << (int)s.unlockedSmuggleHoldMk << "\n";
  f << "unlockedFuelScoopMk " << (int)s.unlockedFuelScoopMk << "\n";
  f << "unlockedWeaponMask " << s.unlockedWeaponMask << "\n";


  // Ship meta / progression
  f << "fuel " << s.fuel << "\n";
  f << "fuelMax " << s.fuelMax << "\n";
  f << "fsdRangeLy " << s.fsdRangeLy << "\n";
  f << "hull " << s.hull << "\n";
  f << "shield " << s.shield << "\n";
  f << "heat " << s.heat << "\n";

  // Countermeasures (ammo / consumables)
  f << "cmFlares " << (int)s.cmFlares << "\n";
  f << "cmChaff " << (int)s.cmChaff << "\n";
  f << "cmHeatSinks " << (int)s.cmHeatSinks << "\n";

  f << "pipsEng " << s.pipsEng << "\n";
  f << "pipsWep " << s.pipsWep << "\n";
  f << "pipsSys " << s.pipsSys << "\n";
  f << "capEngFrac " << s.capEngFrac << "\n";
  f << "capWepFrac " << s.capWepFrac << "\n";
  f << "capSysFrac " << s.capSysFrac << "\n";

  f << "cargoCapacityKg " << s.cargoCapacityKg << "\n";
  f << "passengerSeats " << s.passengerSeats << "\n";
  // Search & Rescue: rescued escape pods carried onboard.
  f << "rescued_pods " << s.rescuedPods.size() << "\n";
  for (const auto& p : s.rescuedPods) {
    f << "rescued_pod "
      << p.id << " "
      << p.recoveredSystem << " "
      << p.registryFactionId << " "
      << p.recoveredDay << " "
      << p.lifeSupportEndDay << " "
      << (p.fromPlayerKill ? 1 : 0)
      << "\n";
  }

  f << "fsdReadyDay " << s.fsdReadyDay << "\n";

  // Navigation UI state (quality-of-life).
  f << "navAutoRun " << (s.navAutoRun ? 1 : 0) << "\n";
  f << "navRouteHop " << s.navRouteHop << "\n";
  f << "navRouteMode " << (int)s.navRouteMode << "\n";
  f << "navConstrainToCurrentFuelRange " << (s.navConstrainToCurrentFuelRange ? 1 : 0) << "\n";
  f << "navHazardWeight " << s.navHazardWeight << "\n";
  f << "pendingArrivalStation " << s.pendingArrivalStation << "\n";

  f << "navRoute " << s.navRoute.size() << "\n";
  for (auto sysId : s.navRoute) {
    f << "nav " << sysId << "\n";
  }

  // Time Trials (leaderboard + optional best-run ghost recordings).
  f << "timeTrialBestRecords " << s.timeTrialBest.size() << "\n";
  for (const auto& r : s.timeTrialBest) {
    f << "timeTrialBest "
      << r.courseKey << " "
      << r.bestTimeSec << " "
      << r.ghostCodec << " "
      << (r.ghostB64.empty() ? "~" : r.ghostB64)
      << "\n";
  }

  // Loadout / progression
  f << "shipHull " << (int)s.shipHull << "\n";
  f << "thrusterMk " << (int)s.thrusterMk << "\n";
  f << "shieldMk " << (int)s.shieldMk << "\n";
  f << "distributorMk " << (int)s.distributorMk << "\n";
  f << "weaponPrimary " << (int)s.weaponPrimary << "\n";
  f << "weaponSecondary " << (int)s.weaponSecondary << "\n";
  f << "weaponAmmoPrimary " << (int)s.weaponAmmoPrimary << "\n";
  f << "weaponAmmoSecondary " << (int)s.weaponAmmoSecondary << "\n";
  f << "smuggleHoldMk " << (int)s.smuggleHoldMk << "\n";
  f << "fuelScoopMk " << (int)s.fuelScoopMk << "\n";

  f << "cargo";
  for (double u : s.cargo) f << " " << u;
  f << "\n";

  // Exploration
  f << "explorationDataCr " << s.explorationDataCr << "\n";
  f << "scannedKeys " << s.scannedKeys.size() << "\n";
  for (core::u64 k : s.scannedKeys) {
    f << "scan " << k << "\n";
  }

  f << "logbook " << s.logbook.size() << "\n";
  for (const auto& e : s.logbook) {
    f << "log "
      << e.key << " "
      << static_cast<int>(e.kind) << " "
      << static_cast<int>(e.subKind) << " "
      << e.systemId << " "
      << e.stationId << " "
      << e.objectId << " "
      << static_cast<int>(e.commodity) << " "
      << e.units << " "
      << e.discoveredDay << " "
      << e.valueCr << " "
      << (e.sold ? 1 : 0)
      << "\n";
  }

  auto b64tok = [](const std::string& s) -> std::string {
    return s.empty() ? std::string("~") : core::base64Encode(s);
  };

  // Comms (diegetic inbox)
  f << "comms " << s.comms.size() << "\n";
  for (const auto& m : s.comms) {
    f << "comms_msg "
      << m.id << " "
      << m.timeDays << " "
      << (int)m.channel << " "
      << m.factionId << " "
      << m.systemId << " "
      << m.stationId << " "
      << m.relatedId << " "
      << (m.unread ? 1 : 0) << " "
      << (m.pinned ? 1 : 0) << " "
      << b64tok(m.from) << " "
      << b64tok(m.subject) << " "
      << b64tok(m.body)
      << "\n";
  }

  // GalNet watchlist (system event bulletins)
  f << "galnetWatchSystems " << s.galNetWatchSystems.size() << "\n";
  for (auto sysId : s.galNetWatchSystems) {
    f << "watch " << sysId << "\n";
  }


  // Integration Hub automations (stellar_game)
  //
  // Persisted so Integration Hub (IFTTT-style) rules survive quicksave/quickload.
  f << "hubAutomationsEnabled " << (s.hubAutomationsEnabled ? 1 : 0) << "\n";
  f << "hubRules " << s.hubAutomationRules.size() << "\n";
  for (const auto& r : s.hubAutomationRules) {
    f << "hub_rule "
      << (r.enabled ? 1 : 0) << " "
      << (int)r.eventKind << " "
      << (int)r.tagMatch << " "
      << r.cooldownSec << " "
      << b64tok(r.name) << " "
      << b64tok(r.tagText) << " "
      << r.actions.size()
      << "\n";

    for (const auto& a : r.actions) {
      f << "hub_action "
        << (int)a.kind << " "
        << (int)a.u64aSource << " "
        << a.u64aConst << " "
        << (int)a.u64bSource << " "
        << a.u64bConst << " "
        << a.i32Const << " "
        << (a.bConst ? 1 : 0) << " "
        << a.delaySec << " "
        << b64tok(a.msgTemplate)
        << "\n";
    }
  }
  // World state (signals / mining depletion)
  f << "resolved_signals " << s.resolvedSignalIds.size() << "\n";
  for (core::u64 id : s.resolvedSignalIds) {
    f << "signal_resolved " << id << "\n";
  }

  f << "asteroid_states " << s.asteroidStates.size() << "\n";
  for (const auto& a : s.asteroidStates) {
    f << "asteroid " << a.asteroidId << " " << a.remainingUnits << "\n";
  }

  // Missions
  f << "nextMissionId " << s.nextMissionId << "\n";
  f << "missions " << s.missions.size() << "\n";
  for (const auto& m : s.missions) {
    f << "mission "
      << m.id << " "
      << static_cast<int>(m.type) << " "
      << m.fromSystem << " " << m.fromStation << " "
      << m.toSystem << " " << m.toStation << " "
      << static_cast<int>(m.commodity) << " "
      << m.units << " "
      << m.targetNpcId << " "
      << m.reward << " "
      << m.deadlineDay << " "
      << (m.completed ? 1 : 0) << " "
      << (m.failed ? 1 : 0) << " "
      << (m.cargoProvided ? 1 : 0) << " "
      // Optional / newer fields (kept at end for backward compatibility)
      << m.factionId << " "
      << m.viaSystem << " " << m.viaStation << " "
      << static_cast<int>(m.leg) << " "
      << (m.scanned ? 1 : 0) << " "
      // Mission provenance (Round 17+) -- appended for backward compatibility
      << static_cast<int>(m.sourceKind) << " "
      << static_cast<int>(m.sourceEventKind) << " "
      << m.sourceEventSeverity01 << " "
      << m.sourceEventEndDay << " "
      << (m.priority ? 1 : 0)
      << "\n";
  }

  // Escort convoys (mission-critical NPC persistence).
  // Stored as lightweight kinematic + health state so escort missions survive save/load.
  f << "escort_convoys " << s.escortConvoys.size() << "\n";
  for (const auto& c : s.escortConvoys) {
    f << "convoy "
      << c.convoyId << " "
      << c.missionId << " "
      << c.systemId << " "
      << c.fromStation << " "
      << c.toStation << " "
      << c.posKm.x << " " << c.posKm.y << " " << c.posKm.z << " "
      << c.velKmS.x << " " << c.velKmS.y << " " << c.velKmS.z << " "
      << c.orient.w << " " << c.orient.x << " " << c.orient.y << " " << c.orient.z << " "
      << c.angVelRadS.x << " " << c.angVelRadS.y << " " << c.angVelRadS.z << " "
      << c.hullFrac << " "
      << c.shieldFrac << " "
      << c.cargoValueCr << " "
      << c.tooFarSec << " "
      << (c.ambushSpawned ? 1 : 0) << " "
      << c.nextAmbushDays
      << "\n";
  }


  // Bounty targets (mission-critical NPC persistence).
  // Stored as lightweight kinematic + health state so bounty missions don't reset on save/load.
  f << "bounty_targets " << s.bountyTargets.size() << "\n";
  for (const auto& b : s.bountyTargets) {
    f << "bounty_target "
      << b.targetId << " "
      << b.missionId << " "
      << b.systemId << " "
      << b.hideoutStation << " "
      << b.posKm.x << " " << b.posKm.y << " " << b.posKm.z << " "
      << b.velKmS.x << " " << b.velKmS.y << " " << b.velKmS.z << " "
      << b.orient.w << " " << b.orient.x << " " << b.orient.y << " " << b.orient.z << " "
      << b.angVelRadS.x << " " << b.angVelRadS.y << " " << b.angVelRadS.z << " "
      << b.hullFrac << " "
      << b.shieldFrac
      << "\n";
  }

  // Industry orders (station processing queues)
  f << "nextIndustryOrderId " << s.nextIndustryOrderId << "\n";
  f << "industry_orders " << s.industryOrders.size() << "\n";
  for (const auto& o : s.industryOrders) {
    f << "industry "
      << o.id << " "
      << static_cast<int>(o.recipe) << " "
      << o.stationId << " "
      << static_cast<int>(o.inputA) << " " << o.inputAUnits << " "
      << static_cast<int>(o.inputB) << " " << o.inputBUnits << " "
      << static_cast<int>(o.output) << " " << o.outputUnits << " "
      << o.submittedDay << " "
      << o.readyDay << " "
      << (o.claimed ? 1 : 0)
      << "\n";
  }

  f << "trackedMissionId " << s.trackedMissionId << "\n";

  // Mission board (cached offers)
  f << "mission_offers_station " << s.missionOffersStationId << "\n";
  f << "mission_offers_day " << s.missionOffersDayStamp << "\n";
  f << "mission_offers " << s.missionOffers.size() << "\n";
  for (const auto& m : s.missionOffers) {
    f << "offer "
      << m.id << " "
      << static_cast<int>(m.type) << " "
      << m.fromSystem << " " << m.fromStation << " "
      << m.toSystem << " " << m.toStation << " "
      << static_cast<int>(m.commodity) << " "
      << m.units << " "
      << m.targetNpcId << " "
      << m.reward << " "
      << m.deadlineDay << " "
      << (m.completed ? 1 : 0) << " "
      << (m.failed ? 1 : 0) << " "
      << (m.cargoProvided ? 1 : 0) << " "
      // Optional / newer fields (kept at end for backward compatibility)
      << m.factionId << " "
      << m.viaSystem << " " << m.viaStation << " "
      << static_cast<int>(m.leg) << " "
      << (m.scanned ? 1 : 0) << " "
      // Mission provenance (Round 17+) -- appended for backward compatibility
      << static_cast<int>(m.sourceKind) << " "
      << static_cast<int>(m.sourceEventKind) << " "
      << m.sourceEventSeverity01 << " "
      << m.sourceEventEndDay << " "
      << (m.priority ? 1 : 0)
      << "\n";
  }

  // Reputation
  f << "reputation " << s.reputation.size() << "\n";
  for (const auto& r : s.reputation) {
    f << "rep " << r.factionId << " " << r.rep << "\n";
  }

  // Bounties
  f << "bounties " << s.bounties.size() << "\n";
  for (const auto& b : s.bounties) {
    f << "bounty " << b.factionId << " " << b.bountyCr << "\n";
  }

  // Fines (minor offenses, payable later)
  f << "fines " << s.fines.size() << "\n";
  for (const auto& fi : s.fines) {
    f << "fine " << fi.factionId << " " << fi.fineCr << " " << fi.dueDay << "\n";
  }

  // Bounty vouchers (earned from destroying criminals; redeemed later)
  f << "bounty_vouchers " << s.bountyVouchers.size() << "\n";
  for (const auto& v : s.bountyVouchers) {
    f << "voucher " << v.factionId << " " << v.bountyCr << "\n";
  }

  // Background NPC traffic day stamps (per visited system)
  f << "trafficStamps " << s.trafficStamps.size() << "\n";
  for (const auto& t : s.trafficStamps) {
    f << "traffic " << t.systemId << " " << t.dayStamp << "\n";
  }

  // Recent NPC trade shipments (TrafficLedger replay).
  //
  // Stored as a flat list with systemId on each record to keep the format simple
  // and robust across version changes.
  f << "trafficShipments " << s.trafficShipments.size() << "\n";
  for (const auto& sh : s.trafficShipments) {
    f << "shipment "
      << sh.id << " "
      << sh.systemId << " "
      << sh.dayStamp << " "
      << sh.fromStation << " "
      << sh.toStation << " "
      << sh.factionId << " "
      << static_cast<int>(sh.commodity) << " "
      << sh.units << " "
      << sh.departDay << " "
      << sh.arriveDay << " "
      << sh.distKm << " "
      << sh.speedKmS
      << "\n";
  }

  // Recently disrupted traffic convoys (anti-farm + economy impact persistence).
  f << "trafficInterdictions " << s.trafficInterdictions.size() << "\n";
  for (const auto& d : s.trafficInterdictions) {
    f << "interdict "
      << d.convoyId << " "
      << d.systemId << " "
      << d.fromStation << " "
      << d.toStation << " "
      << static_cast<int>(d.commodity) << " "
      << d.units << " "
      << d.expireDay
      << "\n";
  }

  // Ambient traffic escort contract (accepted on-grid from traffic convoys).
  // Persisted so quicksave/quickload doesn't insta-fail the contract.
  {
    const auto& c = s.trafficEscort;
    f << "trafficEscortContract "
      << (c.active ? 1 : 0) << " "
      << c.convoyId << " "
      << c.payerFactionId << " "
      << c.toStationId << " "
      << c.startDays << " "
      << c.untilDays << " "
      << c.maxRangeKm << " "
      << c.tooFarSec << " "
      << c.rewardCr << " "
      << c.bonusPerPirateCr << " "
      << c.repReward << " "
      << c.piratesKilled << " "
      << (c.piratesPresentAtStart ? 1 : 0)
      << "\n";
  }

  // Anti-farm record for completed traffic escorts.
  f << "trafficEscortSettlements " << s.trafficEscortSettlements.size() << "\n";
  for (const auto& e : s.trafficEscortSettlements) {
    f << "escortSettled " << e.convoyId << " " << e.settledDay << "\n";
  }

  // Per-system dynamic security deltas.
  // These are applied by stellar_game on top of the deterministic security model.
  f << "systemSecurityDeltas " << s.systemSecurityDeltas.size() << "\n";
  for (const auto& d : s.systemSecurityDeltas) {
    f << "secDelta "
      << d.systemId << " "
      << d.securityDelta << " "
      << d.piracyDelta << " "
      << d.trafficDelta << " "
      << d.lastUpdateDay
      << "\n";
  }


  // Station warehouse / cargo storage (player-owned, per station).
  f << "station_storage " << s.stationStorage.size() << "\n";
  for (const auto& st : s.stationStorage) {
    f << "storage "
      << st.stationId << " "
      << static_cast<int>(st.stationType) << " "
      << st.factionId << " "
      << st.feeRate << " "
      << st.lastFeeDay << " "
      << st.feesDueCr;
    for (double u : st.cargo) f << " " << u;
    f << "\n";
  }

  f << "station_overrides " << s.stationOverrides.size() << "\n";
  for (const auto& ov : s.stationOverrides) {
    f << "station " << ov.stationId << "\n";
    f << "lastUpdateDay " << ov.state.lastUpdateDay << "\n";
    f << "lastSampleDay " << ov.state.lastSampleDay << "\n";

    f << "inventory";
    for (double v : ov.state.inventory) f << " " << v;
    f << "\n";

    // history per commodity
    for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
      const auto& hist = ov.state.history[i];
      f << "history " << i << " " << hist.size();
      for (const auto& p : hist) f << " " << p.day << " " << p.price;
      f << "\n";
    }

    f << "endstation\n";
  }

  if (!f.good()) {
    stellar::core::log(stellar::core::LogLevel::Error, "SaveGame: failed to serialize save: " + path);
    return false;
  }

  const std::string payload = f.str();
  std::string err;
  const bool ok = stellar::core::atomicWriteFile(
      dstPath,
      [&](std::ostream& out, std::string* outErr) -> bool {
        out.write(payload.data(), static_cast<std::streamsize>(payload.size()));
        if (!out.good()) {
          if (outErr) *outErr = "ostream write failed";
          return false;
        }
        return true;
      },
      &err);

  if (!ok) {
    stellar::core::log(stellar::core::LogLevel::Error,
                       "SaveGame: failed to write '" + dstPath.string() + "': " +
                           (err.empty() ? std::string("unknown error") : err));
    return false;
  }

  return true;
}

static bool expectToken(std::istream& in, const char* tok) {
  std::string s;
  if (!(in >> s)) return false;
  return s == tok;
}

bool loadFromFile(const std::string& path, SaveGame& out) {
  std::ifstream f(path);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Warn, "SaveGame: file not found: " + path);
    return false;
  }

  std::string header;
  if (!(f >> header)) return false;
  if (header != "StellarForgeSave") {
    stellar::core::log(stellar::core::LogLevel::Error, "SaveGame: bad header");
    return false;
  }

  int version = 0;
  if (!(f >> version)) return false;
  out = SaveGame{};
  out.version = version;

  std::string key;
  while (f >> key) {
    if (key == "seed") {
      f >> out.seed;
    } else if (key == "timeDays") {
      f >> out.timeDays;
    } else if (key == "currentSystem") {
      f >> out.currentSystem;
    } else if (key == "dockedStation") {
      f >> out.dockedStation;
    } else if (key == "lastDockedSystem") {
      f >> out.lastDockedSystem;
    } else if (key == "lastDockedStation") {
      f >> out.lastDockedStation;
    } else if (key == "shipPosKm") {
      f >> out.shipPosKm.x >> out.shipPosKm.y >> out.shipPosKm.z;
    } else if (key == "shipVelKmS") {
      f >> out.shipVelKmS.x >> out.shipVelKmS.y >> out.shipVelKmS.z;
    } else if (key == "shipOrient") {
      f >> out.shipOrient.w >> out.shipOrient.x >> out.shipOrient.y >> out.shipOrient.z;
    } else if (key == "shipAngVel") {
      f >> out.shipAngVelRadS.x >> out.shipAngVelRadS.y >> out.shipAngVelRadS.z;
    } else if (key == "credits") {
      f >> out.credits;
    } else if (key == "insuranceDebtCr") {
      f >> out.insuranceDebtCr;
    } else if (key == "creditsHighWaterCr") {
      f >> out.creditsHighWaterCr;
      out.creditsHighWaterCr = std::max(0.0, out.creditsHighWaterCr);
    } else if (key == "repHighWater") {
      f >> out.repHighWater;
      out.repHighWater = std::clamp(out.repHighWater, 0.0, 100.0);
    } else if (key == "unlockedHullMask") {
      f >> out.unlockedHullMask;
    } else if (key == "unlockedCoreMk") {
      int v = 1;
      f >> v;
      out.unlockedCoreMk = core::clampCast<core::u8>(v, 1, 3);
    } else if (key == "unlockedSmuggleHoldMk") {
      int v = 0;
      f >> v;
      out.unlockedSmuggleHoldMk = core::clampCast<core::u8>(v, 0, 3);
    } else if (key == "unlockedFuelScoopMk") {
      int v = 0;
      f >> v;
      out.unlockedFuelScoopMk = core::clampCast<core::u8>(v, 0, 3);
    } else if (key == "unlockedWeaponMask") {
      f >> out.unlockedWeaponMask;
      // Ensure starter weapons remain available even if a file is edited by hand.
      out.unlockedWeaponMask |= ((1u << 0) | (1u << 2));
    } else if (key == "fuel") {
      f >> out.fuel;
    } else if (key == "fuelMax") {
      f >> out.fuelMax;
    } else if (key == "fsdRangeLy") {
      f >> out.fsdRangeLy;
    } else if (key == "hull") {
      f >> out.hull;
    } else if (key == "shield") {
      f >> out.shield;
    } else if (key == "heat") {
      f >> out.heat;
    } else if (key == "cmFlares") {
      int v = 0;
      f >> v;
      out.cmFlares = core::clampCast<core::u8>(v, 0, 255);
    } else if (key == "cmChaff") {
      int v = 0;
      f >> v;
      out.cmChaff = core::clampCast<core::u8>(v, 0, 255);
    } else if (key == "cmHeatSinks") {
      int v = 0;
      f >> v;
      out.cmHeatSinks = core::clampCast<core::u8>(v, 0, 255);
    } else if (key == "pipsEng") {
      f >> out.pipsEng;
      out.pipsEng = std::clamp(out.pipsEng, 0, 4);
    } else if (key == "pipsWep") {
      f >> out.pipsWep;
      out.pipsWep = std::clamp(out.pipsWep, 0, 4);
    } else if (key == "pipsSys") {
      f >> out.pipsSys;
      out.pipsSys = std::clamp(out.pipsSys, 0, 4);
    } else if (key == "capEngFrac") {
      f >> out.capEngFrac;
      out.capEngFrac = std::clamp(out.capEngFrac, 0.0, 1.0);
    } else if (key == "capWepFrac") {
      f >> out.capWepFrac;
      out.capWepFrac = std::clamp(out.capWepFrac, 0.0, 1.0);
    } else if (key == "capSysFrac") {
      f >> out.capSysFrac;
      out.capSysFrac = std::clamp(out.capSysFrac, 0.0, 1.0);

    } else if (key == "cargoCapacityKg") {
      f >> out.cargoCapacityKg;
    } else if (key == "passengerSeats") {
      f >> out.passengerSeats;
    } else if (key == "rescued_pods") {
      std::size_t n = 0;
      f >> n;
      out.rescuedPods.clear();
      out.rescuedPods.reserve(std::min<std::size_t>(n, 4096));

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "rescued_pod") {
          f.clear();
          f.seekg(pos);
          break;
        }

        RescuedPod p{};
        int fromKill = 0;
        if (!(f >> p.id >> p.recoveredSystem >> p.registryFactionId >> p.recoveredDay >> p.lifeSupportEndDay >> fromKill)) break;
        p.fromPlayerKill = (fromKill != 0);

        // Defensive clamps for partially-written / corrupt saves.
        if (!std::isfinite(p.recoveredDay)) p.recoveredDay = 0.0;
        if (!std::isfinite(p.lifeSupportEndDay)) p.lifeSupportEndDay = 0.0;
        if (p.recoveredDay < 0.0) p.recoveredDay = 0.0;

        out.rescuedPods.push_back(p);
      }
    } else if (key == "fsdReadyDay") {
      f >> out.fsdReadyDay;
    } else if (key == "navAutoRun") {
      int v = 0;
      f >> v;
      out.navAutoRun = (v != 0);
    } else if (key == "navRouteHop") {
      f >> out.navRouteHop;
    } else if (key == "navRouteMode") {
      int v = 0;
      f >> v;
      out.navRouteMode = core::clampCast<core::u8>(v, 0, 2);
    } else if (key == "navConstrainToCurrentFuelRange") {
      int v = 0;
      f >> v;
      out.navConstrainToCurrentFuelRange = (v != 0);
    } else if (key == "navHazardWeight") {
      f >> out.navHazardWeight;
    } else if (key == "pendingArrivalStation") {
      f >> out.pendingArrivalStation;
    } else if (key == "navRoute") {
      std::size_t n = 0;
      f >> n;
      out.navRoute.clear();
      out.navRoute.reserve(std::min<std::size_t>(n, 100000));
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "nav") {
          f.clear();
          f.seekg(pos);
          break;
        }
        SystemId id = 0;
        if (!(f >> id)) break;
        out.navRoute.push_back(id);
      }
    } else if (key == "timeTrialBestRecords") {
      std::size_t n = 0;
      f >> n;
      out.timeTrialBest.clear();
      out.timeTrialBest.reserve(std::min<std::size_t>(n, 10000));
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "timeTrialBest") {
          f.clear();
          f.seekg(pos);
          break;
        }

        TimeTrialBestRecord r{};
        int codec = 0;
        std::string b64;
        if (!(f >> r.courseKey >> r.bestTimeSec >> codec >> b64)) break;

        if (!std::isfinite(r.bestTimeSec) || r.bestTimeSec < 0.0) r.bestTimeSec = 0.0;
        r.ghostCodec = (core::u32)std::clamp(codec, 0, 1000000);
        if (b64 != "~") r.ghostB64 = b64;

        out.timeTrialBest.push_back(std::move(r));
      }
    } else if (key == "shipHull") {
      int v = 0;
      f >> v;
      out.shipHull = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "thrusterMk") {
      int v = 0;
      f >> v;
      out.thrusterMk = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "shieldMk") {
      int v = 0;
      f >> v;
      out.shieldMk = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "distributorMk") {
      int v = 0;
      f >> v;
      out.distributorMk = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "weaponPrimary") {
      int v = 0;
      f >> v;
      out.weaponPrimary = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "weaponSecondary") {
      int v = 0;
      f >> v;
      out.weaponSecondary = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "weaponAmmoPrimary") {
      int v = 255;
      f >> v;
      out.weaponAmmoPrimary = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "weaponAmmoSecondary") {
      int v = 255;
      f >> v;
      out.weaponAmmoSecondary = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "smuggleHoldMk") {
      int v = 0;
      f >> v;
      out.smuggleHoldMk = (core::u8)std::clamp(v, 0, 255);
    } else if (key == "fuelScoopMk") {
      int v = 0;
      f >> v;
      out.fuelScoopMk = (core::u8)std::clamp(v, 0, 3);
    } else if (key == "cargo") {
      // Cargo line can grow over time as new commodities are added.
      // Read the rest of the line and fill missing entries with 0 so older saves stay loadable.
      std::string line;
      std::getline(f, line);
      std::istringstream iss(line);
      for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
        double v = 0.0;
        if (iss >> v) out.cargo[i] = v;
        else out.cargo[i] = 0.0;
      }
    } else if (key == "explorationDataCr") {
      f >> out.explorationDataCr;
    } else if (key == "scannedKeys") {
      std::size_t n = 0;
      f >> n;
      out.scannedKeys.clear();
      out.scannedKeys.reserve(n);

      // Read scan entries robustly.
      // If the count is stale/corrupt, fall back to stopping when we no longer
      // see a "scan" tag, *without* consuming the next top-level key.
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "scan") {
          f.clear();
          f.seekg(pos);
          break;
        }
        core::u64 k = 0;
        if (!(f >> k)) break;
        out.scannedKeys.push_back(k);
      }
    } else if (key == "logbook") {
      std::size_t n = 0;
      f >> n;
      out.logbook.clear();
      out.logbook.reserve(std::min<std::size_t>(n, 50000));

      // Robust read strategy similar to scannedKeys: stop when tag no longer matches.
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "log") {
          f.clear();
          f.seekg(pos);
          break;
        }

        LogbookEntry e{};
        int kind = 0;
        int sub = 0;
        int commodity = 0;
        int sold = 0;

        if (!(f >> e.key >> kind >> sub >> e.systemId >> e.stationId >> e.objectId
              >> commodity >> e.units >> e.discoveredDay >> e.valueCr >> sold)) {
          break;
        }

        e.kind = static_cast<LogbookEntryKind>(core::clampCast<core::u8>(kind, 0, (int)LogbookEntryKind::Count - 1));
        e.subKind = core::clampCast<core::u8>(sub, 0, 255);
        e.commodity = static_cast<econ::CommodityId>(core::clampCast<core::u8>(commodity, 0, (int)econ::kCommodityCount - 1));
        e.units = std::max(0.0, e.units);
        e.valueCr = std::max(0.0, e.valueCr);
        e.sold = (sold != 0);
        out.logbook.push_back(e);
      }
    } else if (key == "comms") {
      std::size_t n = 0;
      f >> n;
      out.comms.clear();
      out.comms.reserve(std::min<std::size_t>(n, 4096));

      auto decodeTok = [&](const std::string& tok) -> std::string {
        if (tok == "~") return {};
        std::string s;
        if (!core::base64Decode(tok, &s)) return {};
        return s;
      };

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "comms_msg") {
          f.clear();
          f.seekg(pos);
          break;
        }

        CommsMessage m{};
        int chan = 0;
        int unread = 1;
        int pinned = 0;
        std::string fromB64;
        std::string subjectB64;
        std::string bodyB64;

        if (!(f >> m.id >> m.timeDays >> chan >> m.factionId >> m.systemId >> m.stationId >> m.relatedId >> unread >> pinned >> fromB64 >> subjectB64 >> bodyB64)) {
          break;
        }

        m.channel = static_cast<CommsChannel>(core::clampCast<core::u8>(chan, 0, (int)CommsChannel::Custom));
        m.unread = (unread != 0);
        m.pinned = (pinned != 0);
        m.from = decodeTok(fromB64);
        m.subject = decodeTok(subjectB64);
        m.body = decodeTok(bodyB64);

        out.comms.push_back(std::move(m));
      }
    } else if (key == "galnetWatchSystems") {
      std::size_t n = 0;
      f >> n;
      out.galNetWatchSystems.clear();
      out.galNetWatchSystems.reserve(std::min<std::size_t>(n, 4096));

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "watch") {
          f.clear();
          f.seekg(pos);
          break;
        }
        SystemId id = 0;
        if (!(f >> id)) break;
        if (id != 0) out.galNetWatchSystems.push_back(id);
      }

      std::sort(out.galNetWatchSystems.begin(), out.galNetWatchSystems.end());
      out.galNetWatchSystems.erase(
          std::unique(out.galNetWatchSystems.begin(), out.galNetWatchSystems.end()),
          out.galNetWatchSystems.end());
    } else if (key == "hubAutomationsEnabled") {
      int v = 1;
      f >> v;
      out.hubAutomationsEnabled = (v != 0);
    } else if (key == "hubRules") {
      std::size_t n = 0;
      f >> n;
      out.hubAutomationRules.clear();
      out.hubAutomationRules.reserve(std::min<std::size_t>(n, 256));

      auto decodeTok = [&](const std::string& tok) -> std::string {
        if (tok == "~") return {};
        std::string s;
        if (!core::base64Decode(tok, &s)) return {};
        return s;
      };

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "hub_rule") {
          f.clear();
          f.seekg(pos);
          break;
        }

        IntegrationHubAutomationRuleState r{};
        int enabled = 0;
        int eventKind = 0;
        int tagMatch = 0;
        double cooldownSec = 0.0;
        std::string nameB64;
        std::string tagTextB64;
        std::size_t actionCount = 0;

        if (!(f >> enabled >> eventKind >> tagMatch >> cooldownSec >> nameB64 >> tagTextB64 >> actionCount)) {
          break;
        }

        r.enabled = (enabled != 0);
        r.eventKind = core::clampCast<core::u8>(eventKind, 0, 255);
        r.tagMatch = core::clampCast<core::u8>(tagMatch, 0, 255);
        r.cooldownSec = std::max(0.0, cooldownSec);
        r.name = decodeTok(nameB64);
        r.tagText = decodeTok(tagTextB64);

        r.actions.clear();
        r.actions.reserve(std::min<std::size_t>(actionCount, 64));
        for (std::size_t j = 0; j < actionCount; ++j) {
          const std::streampos pos2 = f.tellg();
          std::string atag;
          if (!(f >> atag)) break;
          if (atag != "hub_action") {
            f.clear();
            f.seekg(pos2);
            break;
          }

          IntegrationHubAutomationActionState a{};
          int kind = 0;
          int u64aSource = 0;
          int u64bSource = 0;
          int bConst = 0;
          std::string msgB64;

          if (!(f >> kind >> u64aSource >> a.u64aConst >> u64bSource >> a.u64bConst >> a.i32Const >> bConst >> a.delaySec >> msgB64)) {
            break;
          }

          a.kind = core::clampCast<core::u8>(kind, 0, 255);
          a.u64aSource = core::clampCast<core::u8>(u64aSource, 0, 255);
          a.u64bSource = core::clampCast<core::u8>(u64bSource, 0, 255);
          a.bConst = (bConst != 0);
          a.delaySec = std::max(0.0, a.delaySec);
          a.msgTemplate = decodeTok(msgB64);

          r.actions.push_back(std::move(a));
        }

        out.hubAutomationRules.push_back(std::move(r));
      }
    } else if (key == "resolved_signals") {
      std::size_t n = 0;
      f >> n;
      out.resolvedSignalIds.clear();
      out.resolvedSignalIds.reserve(std::min<std::size_t>(n, 8192));

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "signal_resolved") {
          f.clear();
          f.seekg(pos);
          break;
        }
        core::u64 id = 0;
        if (!(f >> id)) break;
        out.resolvedSignalIds.push_back(id);
      }
    } else if (key == "asteroid_states") {
      std::size_t n = 0;
      f >> n;
      out.asteroidStates.clear();
      out.asteroidStates.reserve(std::min<std::size_t>(n, 16384));

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "asteroid") {
          f.clear();
          f.seekg(pos);
          break;
        }
        AsteroidState st{};
        if (!(f >> st.asteroidId >> st.remainingUnits)) break;
        out.asteroidStates.push_back(st);
      }
    } else if (key == "trackedMissionId") {
      f >> out.trackedMissionId;
    } else if (key == "nextMissionId") {
      f >> out.nextMissionId;
    } else if (key == "nextIndustryOrderId") {
      f >> out.nextIndustryOrderId;
    } else if (key == "missions") {
      std::size_t n = 0;
      f >> n;
      out.missions.clear();
      out.missions.reserve(std::min<std::size_t>(n, 4096));
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tok;
        if (!(f >> tok)) break;
        if (tok != "mission") {
          f.clear();
          f.seekg(pos);
          break;
        }

        // Parse the rest of the mission line so we can be backward compatible
        // with older saves that did not include newer fields.
        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        Mission m{};
        int type = 0;
        int commodity = 0;
        int completed = 0;
        int failed = 0;
        int cargoProvided = 0;

        iss >> m.id
            >> type
            >> m.fromSystem >> m.fromStation
            >> m.toSystem >> m.toStation
            >> commodity
            >> m.units
            >> m.targetNpcId
            >> m.reward
            >> m.deadlineDay
            >> completed
            >> failed
            >> cargoProvided;

        m.type = static_cast<MissionType>(type);
        m.commodity = static_cast<econ::CommodityId>(commodity);
        m.completed = (completed != 0);
        m.failed = (failed != 0);
        m.cargoProvided = (cargoProvided != 0);

        // Optional fields (save version >= 3)
        int leg = 0;
        int scanned = 0;
        if (iss >> m.factionId >> m.viaSystem >> m.viaStation >> leg >> scanned) {
          m.leg = static_cast<core::u8>(std::clamp(leg, 0, 255));
          m.scanned = (scanned != 0);

          // Optional: mission provenance / priority (Round 17+).
          int sourceKind = 0;
          int sourceEventKind = 0;
          double sev = 0.0;
          double endDay = 0.0;
          int priority = 0;

          if (iss >> sourceKind >> sourceEventKind >> sev >> endDay >> priority) {
            m.sourceKind = static_cast<MissionSourceKind>(std::clamp(sourceKind, 0, 255));
            m.sourceEventKind = static_cast<SystemEventKind>(std::clamp(sourceEventKind, 0, 255));
            m.sourceEventSeverity01 = std::clamp(sev, 0.0, 1.0);
            m.sourceEventEndDay = std::isfinite(endDay) ? endDay : 0.0;
            m.priority = (priority != 0);
          }
        }

        out.missions.push_back(std::move(m));
      }
    } else if (key == "escort_convoys") {
      std::size_t n = 0;
      f >> n;
      out.escortConvoys.clear();
      out.escortConvoys.reserve(std::min<std::size_t>(n, 1024));
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tok;
        if (!(f >> tok)) break;
        if (tok != "convoy") {
          f.clear();
          f.seekg(pos);
          break;
        }

        // Parse the rest of the convoy line so we can extend fields over time.
        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        EscortConvoyState c{};
        int ambushSpawned = 0;
        iss >> c.convoyId
            >> c.missionId
            >> c.systemId
            >> c.fromStation
            >> c.toStation
            >> c.posKm.x >> c.posKm.y >> c.posKm.z
            >> c.velKmS.x >> c.velKmS.y >> c.velKmS.z
            >> c.orient.w >> c.orient.x >> c.orient.y >> c.orient.z
            >> c.angVelRadS.x >> c.angVelRadS.y >> c.angVelRadS.z
            >> c.hullFrac
            >> c.shieldFrac
            >> c.cargoValueCr
            >> c.tooFarSec
            >> ambushSpawned
            >> c.nextAmbushDays;

        c.hullFrac = std::clamp(c.hullFrac, 0.0, 1.0);
        c.shieldFrac = std::clamp(c.shieldFrac, 0.0, 1.0);
        c.tooFarSec = std::max(0.0, c.tooFarSec);
        c.ambushSpawned = (ambushSpawned != 0);

        if (c.convoyId != 0) {
          out.escortConvoys.push_back(std::move(c));
        }
      }
    } else if (key == "bounty_targets") {
      std::size_t n = 0;
      f >> n;
      out.bountyTargets.clear();
      out.bountyTargets.reserve(std::min<std::size_t>(n, 2048));

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tok;
        if (!(f >> tok)) break;
        if (tok != "bounty_target") {
          f.clear();
          f.seekg(pos);
          break;
        }

        // Parse the rest of the target line so we can extend fields over time.
        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        BountyTargetState b{};
        iss >> b.targetId
            >> b.missionId
            >> b.systemId
            >> b.hideoutStation
            >> b.posKm.x >> b.posKm.y >> b.posKm.z
            >> b.velKmS.x >> b.velKmS.y >> b.velKmS.z
            >> b.orient.w >> b.orient.x >> b.orient.y >> b.orient.z
            >> b.angVelRadS.x >> b.angVelRadS.y >> b.angVelRadS.z
            >> b.hullFrac
            >> b.shieldFrac;

        b.hullFrac = std::clamp(b.hullFrac, 0.0, 1.0);
        b.shieldFrac = std::clamp(b.shieldFrac, 0.0, 1.0);

        if (b.targetId != 0) {
          out.bountyTargets.push_back(std::move(b));
        }
      }
    } else if (key == "industry_orders") {
      std::size_t n = 0;
      f >> n;
      out.industryOrders.clear();
      out.industryOrders.reserve(std::min<std::size_t>(n, 4096));
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tok;
        if (!(f >> tok)) break;
        if (tok != "industry") {
          f.clear();
          f.seekg(pos);
          break;
        }

        // Parse the rest of the order line so we can grow fields over time.
        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        IndustryOrder o{};
        int recipe = 0;
        int inA = 0;
        int inB = 0;
        int outC = 0;
        int claimed = 0;

        iss >> o.id
            >> recipe
            >> o.stationId
            >> inA >> o.inputAUnits
            >> inB >> o.inputBUnits
            >> outC >> o.outputUnits
            >> o.submittedDay
            >> o.readyDay
            >> claimed;

        o.recipe = static_cast<IndustryRecipeId>(recipe);
        o.inputA = static_cast<econ::CommodityId>(inA);
        o.inputB = static_cast<econ::CommodityId>(inB);
        o.output = static_cast<econ::CommodityId>(outC);
        o.claimed = (claimed != 0);

        out.industryOrders.push_back(o);
      }
    } else if (key == "mission_offers_station") {
      f >> out.missionOffersStationId;
    } else if (key == "mission_offers_day") {
      f >> out.missionOffersDayStamp;
    } else if (key == "mission_offers") {
      std::size_t n = 0;
      f >> n;
      out.missionOffers.clear();
      out.missionOffers.reserve(std::min<std::size_t>(n, 4096));
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tok;
        if (!(f >> tok)) break;
        if (tok != "offer") {
          f.clear();
          f.seekg(pos);
          break;
        }

        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        Mission m{};
        int type = 0;
        int commodity = 0;
        int completed = 0;
        int failed = 0;
        int cargoProvided = 0;

        iss >> m.id
            >> type
            >> m.fromSystem >> m.fromStation
            >> m.toSystem >> m.toStation
            >> commodity
            >> m.units
            >> m.targetNpcId
            >> m.reward
            >> m.deadlineDay
            >> completed
            >> failed
            >> cargoProvided;

        m.type = static_cast<MissionType>(type);
        m.commodity = static_cast<econ::CommodityId>(commodity);
        m.completed = (completed != 0);
        m.failed = (failed != 0);
        m.cargoProvided = (cargoProvided != 0);

        int leg = 0;
        int scanned = 0;
        if (iss >> m.factionId >> m.viaSystem >> m.viaStation >> leg >> scanned) {
          m.leg = static_cast<core::u8>(std::clamp(leg, 0, 255));
          m.scanned = (scanned != 0);

          // Optional: mission provenance / priority (Round 17+).
          int sourceKind = 0;
          int sourceEventKind = 0;
          double sev = 0.0;
          double endDay = 0.0;
          int priority = 0;

          if (iss >> sourceKind >> sourceEventKind >> sev >> endDay >> priority) {
            m.sourceKind = static_cast<MissionSourceKind>(std::clamp(sourceKind, 0, 255));
            m.sourceEventKind = static_cast<SystemEventKind>(std::clamp(sourceEventKind, 0, 255));
            m.sourceEventSeverity01 = std::clamp(sev, 0.0, 1.0);
            m.sourceEventEndDay = std::isfinite(endDay) ? endDay : 0.0;
            m.priority = (priority != 0);
          }
        }

        out.missionOffers.push_back(std::move(m));
      }
    } else if (key == "reputation") {
      std::size_t n = 0;
      f >> n;
      out.reputation.clear();
      out.reputation.reserve(n);
      for (std::size_t i = 0; i < n; ++i) {
        std::string tok;
        if (!(f >> tok) || tok != "rep") return false;
        FactionReputation r{};
        f >> r.factionId >> r.rep;
        out.reputation.push_back(std::move(r));
      }
    } else if (key == "bounties") {
      std::size_t n = 0;
      f >> n;
      out.bounties.clear();
      out.bounties.reserve(n);

      // Same robustness strategy as scannedKeys.
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "bounty") {
          f.clear();
          f.seekg(pos);
          break;
        }
        FactionBounty b{};
        if (!(f >> b.factionId >> b.bountyCr)) break;
        out.bounties.push_back(std::move(b));
      }
    } else if (key == "fines") {
      std::size_t n = 0;
      f >> n;
      out.fines.clear();
      out.fines.reserve(n);

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "fine") {
          f.clear();
          f.seekg(pos);
          break;
        }
        FactionFine fi{};
        if (!(f >> fi.factionId >> fi.fineCr)) break;
        // dueDay is optional for forward/backward compatibility.
        if (!(f >> fi.dueDay)) {
          f.clear();
          fi.dueDay = 0.0;
        }
        out.fines.push_back(std::move(fi));
      }
    } else if (key == "bounty_vouchers") {
      std::size_t n = 0;
      f >> n;
      out.bountyVouchers.clear();
      out.bountyVouchers.reserve(n);

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "voucher") {
          f.clear();
          f.seekg(pos);
          break;
        }
        FactionBounty v{};
        if (!(f >> v.factionId >> v.bountyCr)) break;
        out.bountyVouchers.push_back(std::move(v));
      }
    } else if (key == "trafficStamps") {
      std::size_t n = 0;
      f >> n;
      out.trafficStamps.clear();
      out.trafficStamps.reserve(n);

      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "traffic") {
          f.clear();
          f.seekg(pos);
          break;
        }
        SystemTrafficStamp t{};
        if (!(f >> t.systemId >> t.dayStamp)) break;
        out.trafficStamps.push_back(std::move(t));
      }
    } else if (key == "trafficShipments") {
      std::size_t n = 0;
      f >> n;
      out.trafficShipments.clear();
      out.trafficShipments.reserve(std::min<std::size_t>(n, 100000));

      // Defensive: allow corrupted or hand-edited save files to contain duplicate
      // shipment ids. These can otherwise surface as duplicated convoy signals.
      std::unordered_set<core::u64> seenIds;
      seenIds.reserve(std::min<std::size_t>(n, 4096) * 2 + 1);

      // Robustness strategy: if the count is stale/corrupt, stop when we no longer
      // see "shipment" and continue parsing later top-level keys.
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "shipment") {
          f.clear();
          f.seekg(pos);
          break;
        }

        // Parse the rest of the shipment line so we can extend fields over time.
        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        TrafficShipmentState sh{};
        int commodity = 0;
        // Mandatory fields.
        if (!(iss >> sh.id
                  >> sh.systemId
                  >> sh.dayStamp
                  >> sh.fromStation
                  >> sh.toStation
                  >> sh.factionId
                  >> commodity
                  >> sh.units)) {
          continue;
        }

        // Optional schedule metadata (older saves may omit).
        iss >> sh.departDay >> sh.arriveDay >> sh.distKm >> sh.speedKmS;

        sh.commodity = static_cast<econ::CommodityId>(std::clamp(commodity, 0, (int)econ::kCommodityCount - 1));

        if (sh.id != 0 && sh.systemId != 0 && sh.fromStation != 0 && sh.toStation != 0) {
          if (!seenIds.insert(sh.id).second) continue;
          out.trafficShipments.push_back(std::move(sh));
        }
      }
    
    } else if (key == "trafficInterdictions") {
      std::size_t n = 0;
      f >> n;
      out.trafficInterdictions.clear();
      out.trafficInterdictions.reserve(std::min<std::size_t>(n, 100000));

      // Defensive: allow corrupted or hand-edited save files to contain duplicate
      // convoy ids. These can otherwise allow re-spawning / farming.
      std::unordered_set<core::u64> seenIds;
      seenIds.reserve(std::min<std::size_t>(n, 4096) * 2 + 1);

      // Robustness strategy: if the count is stale/corrupt, stop when we no longer
      // see "interdict" and continue parsing later top-level keys.
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "interdict") {
          f.clear();
          f.seekg(pos);
          break;
        }

        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        TrafficInterdictionState d{};
        int commodity = 0;

        // Mandatory fields.
        if (!(iss >> d.convoyId
                  >> d.systemId
                  >> d.fromStation
                  >> d.toStation
                  >> commodity
                  >> d.units
                  >> d.expireDay)) {
          continue;
        }

        d.commodity = static_cast<econ::CommodityId>(std::clamp(commodity, 0, (int)econ::kCommodityCount - 1));

        if (d.convoyId != 0 && d.systemId != 0 && d.toStation != 0) {
          if (!seenIds.insert(d.convoyId).second) continue;
          out.trafficInterdictions.push_back(std::move(d));
        }
      }

    } else if (key == "trafficEscortContract") {
      // One-line record: contract state may be absent in older saves.
      int active = 0;
      TrafficEscortContractState c{};
      if (!(f >> active
              >> c.convoyId
              >> c.payerFactionId
              >> c.toStationId
              >> c.startDays
              >> c.untilDays
              >> c.maxRangeKm
              >> c.tooFarSec
              >> c.rewardCr
              >> c.bonusPerPirateCr
              >> c.repReward
              >> c.piratesKilled)) {
        continue;
      }
      int piratesPresent = 0;
      f >> piratesPresent;
      c.active = (active != 0);
      c.piratesPresentAtStart = (piratesPresent != 0);

      // Basic sanity: avoid "active but empty" state when reading corrupted saves.
      if (c.convoyId == 0 || c.maxRangeKm <= 0.0 || c.untilDays < c.startDays) {
        c.active = false;
        c.convoyId = 0;
        c.tooFarSec = 0.0;
      }

      out.trafficEscort = std::move(c);

    } else if (key == "trafficEscortSettlements") {
      std::size_t n = 0;
      f >> n;
      out.trafficEscortSettlements.clear();
      out.trafficEscortSettlements.reserve(std::min<std::size_t>(n, 200000));

      // Defensive: allow duplicates; keep the latest settledDay per convoy id.
      std::unordered_map<core::u64, double> latest;
      latest.reserve(std::min<std::size_t>(n, 4096) * 2 + 1);

      // Robustness: if the count is stale, stop when we no longer see "escortSettled".
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "escortSettled") {
          f.clear();
          f.seekg(pos);
          break;
        }

        TrafficEscortSettlementState e{};
        if (!(f >> e.convoyId >> e.settledDay)) break;
        if (e.convoyId == 0) continue;

        const auto it = latest.find(e.convoyId);
        if (it == latest.end() || e.settledDay > it->second) {
          latest[e.convoyId] = e.settledDay;
        }
      }

      out.trafficEscortSettlements.reserve(latest.size());
      for (const auto& kv : latest) {
        TrafficEscortSettlementState e{};
        e.convoyId = kv.first;
        e.settledDay = kv.second;
        out.trafficEscortSettlements.push_back(e);
      }
      std::sort(out.trafficEscortSettlements.begin(), out.trafficEscortSettlements.end(),
                [](const TrafficEscortSettlementState& a, const TrafficEscortSettlementState& b) {
                  if (a.settledDay != b.settledDay) return a.settledDay < b.settledDay;
                  return a.convoyId < b.convoyId;
                });

    } else if (key == "systemSecurityDeltas") {
      std::size_t n = 0;
      f >> n;
      out.systemSecurityDeltas.clear();
      out.systemSecurityDeltas.reserve(std::min<std::size_t>(n, 100000));

      // Defensive: allow corrupted / hand-edited saves to contain duplicate system ids.
      // Keep only the first occurrence so downstream systems don't double-apply impulses.
      std::unordered_set<SystemId> seenIds;
      seenIds.reserve(std::min<std::size_t>(n, 4096) * 2 + 1);

      // Robustness strategy: if the count is stale/corrupt, stop when we no longer see
      // "secDelta" and continue parsing later top-level keys.
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "secDelta") {
          f.clear();
          f.seekg(pos);
          break;
        }

        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        SystemSecurityDeltaState d{};
        if (!(iss >> d.systemId
                  >> d.securityDelta
                  >> d.piracyDelta
                  >> d.trafficDelta
                  >> d.lastUpdateDay)) {
          continue;
        }

        if (d.systemId != 0) {
          if (!seenIds.insert(d.systemId).second) continue;
          out.systemSecurityDeltas.push_back(std::move(d));
        }
      }
    } else if (key == "station_storage") {
      std::size_t n = 0;
      f >> n;

      out.stationStorage.clear();
      out.stationStorage.reserve(std::min<std::size_t>(n, 4096));

      // Robustness strategy: if the count is stale/corrupt, stop when we no longer see "storage"
      // and continue parsing subsequent top-level keys.
      for (std::size_t i = 0; i < n; ++i) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "storage") {
          f.clear();
          f.seekg(pos);
          break;
        }

        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);

        StationStorage st{};
        int stType = 0;
        iss >> st.stationId >> stType >> st.factionId >> st.feeRate >> st.lastFeeDay >> st.feesDueCr;
        st.stationType = static_cast<econ::StationType>(std::clamp(stType, 0, (int)econ::StationType::Count - 1));

        for (std::size_t ci = 0; ci < econ::kCommodityCount; ++ci) {
          double u = 0.0;
          if (iss >> u) st.cargo[ci] = u;
          else st.cargo[ci] = 0.0;
        }

        out.stationStorage.push_back(std::move(st));
      }
    } else if (key == "station_overrides") {
      std::size_t n = 0;
      f >> n;
      out.stationOverrides.clear();
      out.stationOverrides.reserve(std::min<std::size_t>(n, 2048));

      for (std::size_t si = 0; si < n; ++si) {
        const std::streampos pos = f.tellg();
        std::string tag;
        if (!(f >> tag)) break;
        if (tag != "station") {
          f.clear();
          f.seekg(pos);
          break;
        }

        StationEconomyOverride ov{};
        if (!(f >> ov.stationId)) break;

        // Parse station block until endstation
        std::string sk;
        while (f >> sk) {
          if (sk == "endstation") break;
          if (sk == "lastUpdateDay") {
            f >> ov.state.lastUpdateDay;
          } else if (sk == "lastSampleDay") {
            f >> ov.state.lastSampleDay;
          } else if (sk == "inventory") {
            // Like player cargo, station inventory lines can grow as commodities are added.
            // Parse the remainder of the line and default missing values to 0.
            std::string line;
            std::getline(f, line);
            std::istringstream iss(line);
            for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
              double v = 0.0;
              if (iss >> v) ov.state.inventory[i] = v;
              else ov.state.inventory[i] = 0.0;
            }
          } else if (sk == "history") {
            std::size_t cid = 0;
            std::size_t count = 0;
            f >> cid >> count;
            if (cid >= econ::kCommodityCount) return false;
            auto& hist = ov.state.history[cid];
            hist.clear();
            hist.reserve(count);
            for (std::size_t j = 0; j < count; ++j) {
              econ::PricePoint p{};
              f >> p.day >> p.price;
              hist.push_back(p);
            }
          } else {
            // Unknown token: attempt to skip line
            std::string line;
            std::getline(f, line);
          }
        }

        out.stationOverrides.push_back(std::move(ov));
      }
    } else {
      // Unknown key: skip rest of line
      std::string line;
      std::getline(f, line);
    }
  }

  return true;
}

} // namespace stellar::sim
