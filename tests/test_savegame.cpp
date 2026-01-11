#include "stellar/sim/SaveGame.h"

#include "stellar/sim/ShipLoadout.h"

#include "stellar/econ/Commodity.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

static std::string readAllText(const std::string& path) {
  std::ifstream f(path, std::ios::in);
  if (!f) return {};
  std::ostringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

static bool writeAllText(const std::string& path, const std::string& content) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) return false;
  f << content;
  return static_cast<bool>(f);
}

static std::string replaceCountLine(const std::string& content,
                                   const std::string& key,
                                   std::size_t newCount,
                                   bool& replaced) {
  replaced = false;
  std::istringstream iss(content);
  std::ostringstream oss;
  std::string line;
  const std::string prefix = key + " ";

  while (std::getline(iss, line)) {
    if (!replaced && line.rfind(prefix, 0) == 0) {
      oss << key << " " << newCount << "\n";
      replaced = true;
    } else {
      oss << line << "\n";
    }
  }
  return oss.str();
}

static std::string duplicateSecondLineWithPrefix(const std::string& content,
                                                 const std::string& prefix,
                                                 bool& duplicated) {
  duplicated = false;
  std::istringstream iss(content);
  std::ostringstream oss;
  std::string line;

  std::string firstMatch;
  int seen = 0;

  while (std::getline(iss, line)) {
    if (line.rfind(prefix, 0) == 0) {
      ++seen;
      if (seen == 1) {
        firstMatch = line;
        oss << line << "\n";
        continue;
      }
      if (seen == 2 && !firstMatch.empty()) {
        oss << firstMatch << "\n";
        duplicated = true;
        continue;
      }
    }
    oss << line << "\n";
  }
  return oss.str();
}

int test_savegame() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  SaveGame s{};
  s.version = 10;
  s.seed = 123456789ull;
  s.timeDays = 42.5;
  s.currentSystem = 999001;
  s.dockedStation = 1002003;
  s.credits = 9876.5;
  s.insuranceDebtCr = 321.0;
  s.shipPosKm = {1.0, 2.0, 3.0};
  s.shipVelKmS = {0.1, 0.0, -0.2};
  s.fuel = 12.3;
  s.fuelMax = 20.0;
  s.cargoCapacityKg = 240.0;
  s.passengerSeats = 6;
  // Navigation UI persistence (route plotting + auto-run).
  s.navAutoRun = true;
  s.navRouteMode = 2;
  s.navConstrainToCurrentFuelRange = false;
  s.navRouteHop = 1;
  s.pendingArrivalStation = 424242;
  s.navRoute = {1001, 1002, 1003, 1004};

  // Loadout (ammo-bearing weapons are new as of this patch).
  s.shipHull = (core::u8)ShipHullClass::Fighter;
  s.weaponPrimary = (core::u8)WeaponType::HomingMissile;
  s.weaponSecondary = (core::u8)WeaponType::BeamLaser;
  s.weaponAmmoPrimary = 4;
  s.weaponAmmoSecondary = 0;

  s.smuggleHoldMk = 2;
  s.hull = 0.75;
  s.shield = 0.15;
  s.cargo[static_cast<std::size_t>(econ::CommodityId::Metals)] = 5.5;
  s.cargo[static_cast<std::size_t>(econ::CommodityId::Food)] = 12.0;
  s.scannedKeys = {111, 222, 333};

  // Logbook (scan history + exploration data brokerage state).
  {
    LogbookEntry e{};
    e.key = 0xABC1;
    e.kind = LogbookEntryKind::StarScan;
    e.systemId = s.currentSystem;
    e.discoveredDay = 42.10;
    e.valueCr = 125.0;
    e.sold = false;
    s.logbook.push_back(e);
  }
  {
    LogbookEntry e{};
    e.key = 0xABC2;
    e.kind = LogbookEntryKind::AsteroidProspect;
    e.systemId = 10;
    e.objectId = 0x20000000000000AAull;
    e.commodity = econ::CommodityId::Metals;
    e.units = 42.25;
    e.discoveredDay = 42.20;
    e.valueCr = 0.0;
    e.sold = true;
    s.logbook.push_back(e);
  }
  s.resolvedSignalIds = {0x2000000000000001ull, 0x2000000000000002ull};
  s.asteroidStates.push_back({0x20000000000000AAull, 42.25});
  s.asteroidStates.push_back({0x20000000000000BBull, 0.0});
  s.trackedMissionId = 1;

  // One accepted mission.
  {
    Mission m{};
    m.id = 1;
    m.type = MissionType::Delivery;
    m.factionId = 7;
    m.fromSystem = 10;
    m.fromStation = 11;
    m.toSystem = 20;
    m.toStation = 21;
    m.viaSystem = 15;
    m.viaStation = 16;
    m.leg = 2;
    m.commodity = econ::CommodityId::Food;
    m.units = 8.0;
    m.reward = 1234.0;
    m.deadlineDay = 55.0;
    m.cargoProvided = true;
    m.scanned = true;
    s.missions.push_back(m);
  }

  // Persisted escort convoy state (mission-critical NPC).
  {
    EscortConvoyState c{};
    c.convoyId = 0x40000000000000CCull;
    c.missionId = 1;
    c.systemId = s.currentSystem;
    c.fromStation = 11;
    c.toStation = 21;
    c.posKm = {7.0, 8.0, 9.0};
    c.velKmS = {0.0, 0.1, 0.2};
    c.orient = {0.9238795325, 0.0, 0.3826834323, 0.0};
    c.angVelRadS = {0.01, 0.02, 0.03};
    c.hullFrac = 0.55;
    c.shieldFrac = 0.25;
    c.cargoValueCr = 1234.5;
    c.tooFarSec = 6.25;
    c.ambushSpawned = true;
    c.nextAmbushDays = 42.55;
    s.escortConvoys.push_back(c);
  }

  // Persisted bounty target state (mission-critical NPC).
  {
    BountyTargetState b{};
    b.targetId = 0x50000000000000DDull;
    b.missionId = 1;
    b.systemId = s.currentSystem;
    b.hideoutStation = 21;
    b.posKm = {11.0, 12.0, 13.0};
    b.velKmS = {0.2, 0.0, -0.1};
    b.orient = {0.9659258263, 0.0, 0.2588190451, 0.0};
    b.angVelRadS = {0.02, 0.01, 0.00};
    b.hullFrac = 0.45;
    b.shieldFrac = 0.80;
    s.bountyTargets.push_back(b);
  }

  // Mission offers (board persistence).
  s.missionOffersStationId = 11;
  s.missionOffersDayStamp = 42;
  {
    Mission m{};
    m.id = 2;
    m.type = MissionType::Courier;
    m.factionId = 7;
    m.fromSystem = 10;
    m.fromStation = 11;
    m.toSystem = 12;
    m.toStation = 13;
    m.reward = 250.0;
    m.deadlineDay = 45.0;
    s.missionOffers.push_back(m);
  }

  // Reputation + bounties.
  s.reputation.push_back({7, 12.5});
  s.bounties.push_back({7, 500.0});
  s.bountyVouchers.push_back({7, 200.0});

  // Background traffic day-stamps (per visited system).
  s.trafficStamps.push_back({10, 41});

  // Recent NPC-trade shipments (TrafficLedger replay).
  {
    TrafficShipmentState sh{};
    sh.id = 0x6000000000000011ull;
    sh.systemId = 10;
    sh.dayStamp = 41;
    sh.fromStation = 111;
    sh.toStation = 222;
    sh.factionId = 7;
    sh.commodity = econ::CommodityId::Food;
    sh.units = 12.5;
    sh.departDay = 41.25;
    sh.arriveDay = 41.32;
    sh.distKm = 120000.0;
    sh.speedKmS = 8000.0;
    s.trafficShipments.push_back(sh);
  }
  {
    TrafficShipmentState sh{};
    sh.id = 0x6000000000000012ull;
    sh.systemId = 10;
    sh.dayStamp = 42;
    sh.fromStation = 222;
    sh.toStation = 333;
    sh.factionId = 0;
    sh.commodity = econ::CommodityId::Metals;
    sh.units = 3.75;
    sh.departDay = 42.10;
    sh.arriveDay = 42.22;
    sh.distKm = 250000.0;
    sh.speedKmS = 6000.0;
    s.trafficShipments.push_back(sh);
  }

  // Recently interdicted traffic convoys (anti-farm + economy impact persistence).
  {
    TrafficInterdictionState d{};
    d.convoyId = 0x70000000000000A1ull;
    d.systemId = 10;
    d.fromStation = 111;
    d.toStation = 222;
    d.commodity = econ::CommodityId::Food;
    d.units = 12.5;
    d.expireDay = 41.32;
    s.trafficInterdictions.push_back(d);
  }

  // Ambient traffic escort contract + settlements (anti-farm).
  {
    s.trafficEscort.active = true;
    s.trafficEscort.convoyId = 0x6000000000000E11ull;
    s.trafficEscort.payerFactionId = 7;
    s.trafficEscort.toStationId = 222;
    s.trafficEscort.startDays = 42.40;
    s.trafficEscort.untilDays = 42.46;
    s.trafficEscort.maxRangeKm = 155000.0;
    s.trafficEscort.tooFarSec = 4.25;
    s.trafficEscort.rewardCr = 1200.0;
    s.trafficEscort.bonusPerPirateCr = 220.0;
    s.trafficEscort.repReward = 3.5;
    s.trafficEscort.piratesKilled = 2;
    s.trafficEscort.piratesPresentAtStart = true;
  }
  {
    s.trafficEscortSettlements.push_back(TrafficEscortSettlementState{0x6000000000000E12ull, 41.95});
    s.trafficEscortSettlements.push_back(TrafficEscortSettlementState{0x6000000000000E13ull, 42.10});
  }

  // Per-system security deltas (dynamic world state layered on top of SecurityModel).
  {
    SystemSecurityDeltaState d{};
    d.systemId = 10;
    d.securityDelta = 0.12;
    d.piracyDelta = -0.08;
    d.trafficDelta = 0.03;
    d.lastUpdateDay = 41.0;
    s.systemSecurityDeltas.push_back(d);
  }
  {
    SystemSecurityDeltaState d{};
    d.systemId = 999001;
    d.securityDelta = -0.04;
    d.piracyDelta = 0.11;
    d.trafficDelta = -0.02;
    d.lastUpdateDay = 42.25;
    s.systemSecurityDeltas.push_back(d);
  }
  {
    TrafficInterdictionState d{};
    d.convoyId = 0x70000000000000A2ull;
    d.systemId = 10;
    d.fromStation = 222;
    d.toStation = 333;
    d.commodity = econ::CommodityId::Metals;
    d.units = 3.75;
    d.expireDay = 42.22;
    s.trafficInterdictions.push_back(d);
  }

  // Station economy overrides (cached economies persisted).
  {
    StationEconomyOverride ov{};
    ov.stationId = 9999;
    ov.state.lastUpdateDay = 10.0;
    ov.state.lastSampleDay = 9.5;
    for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
      ov.state.inventory[i] = (double)i * 3.0;
    }
    ov.state.history[0].push_back(econ::PricePoint{1.0, 2.0});
    s.stationOverrides.push_back(ov);
  }

  const std::string path = "savegame_test.sav";

  if (!saveToFile(s, path)) {
    std::cerr << "[test_savegame] saveToFile failed\n";
    return 1;
  }

  SaveGame l{};
  if (!loadFromFile(path, l)) {
    std::cerr << "[test_savegame] loadFromFile failed\n";
    std::filesystem::remove(path);
    return 1;
  }

  // Basic scalar fields.
  if (l.seed != s.seed) {
    std::cerr << "[test_savegame] seed mismatch\n";
    ++fails;
  }
  if (!nearly(l.timeDays, s.timeDays)) {
    std::cerr << "[test_savegame] timeDays mismatch\n";
    ++fails;
  }
  if (l.currentSystem != s.currentSystem || l.dockedStation != s.dockedStation) {
    std::cerr << "[test_savegame] system/station mismatch\n";
    ++fails;
  }
  if (!nearly(l.credits, s.credits)) {
    std::cerr << "[test_savegame] credits mismatch\n";
    ++fails;
  }
  if (!nearly(l.insuranceDebtCr, s.insuranceDebtCr)) {
    std::cerr << "[test_savegame] insuranceDebtCr mismatch\n";
    ++fails;
  }

  if (l.passengerSeats != s.passengerSeats) {
    std::cerr << "[test_savegame] passengerSeats mismatch\n";
    ++fails;
  }

  // Navigation UI state.
  if (l.navAutoRun != s.navAutoRun) {
    std::cerr << "[test_savegame] navAutoRun mismatch\n";
    ++fails;
  }
  if (l.navRouteHop != s.navRouteHop) {
    std::cerr << "[test_savegame] navRouteHop mismatch\n";
    ++fails;
  }
  if (l.navRouteMode != s.navRouteMode) {
    std::cerr << "[test_savegame] navRouteMode mismatch\n";
    ++fails;
  }
  if (l.navConstrainToCurrentFuelRange != s.navConstrainToCurrentFuelRange) {
    std::cerr << "[test_savegame] navConstrainToCurrentFuelRange mismatch\n";
    ++fails;
  }
  if (l.pendingArrivalStation != s.pendingArrivalStation) {
    std::cerr << "[test_savegame] pendingArrivalStation mismatch\n";
    ++fails;
  }
  if (l.navRoute != s.navRoute) {
    std::cerr << "[test_savegame] navRoute mismatch\n";
    ++fails;
  }

  if (l.smuggleHoldMk != s.smuggleHoldMk) {
    std::cerr << "[test_savegame] smuggleHoldMk mismatch\n";
    ++fails;
  }

  // Loadout fields (incl. missile ammo).
  if (l.shipHull != s.shipHull) {
    std::cerr << "[test_savegame] shipHull mismatch\n";
    ++fails;
  }
  if (l.weaponPrimary != s.weaponPrimary || l.weaponSecondary != s.weaponSecondary) {
    std::cerr << "[test_savegame] weapons mismatch\n";
    ++fails;
  }
  if (l.weaponAmmoPrimary != s.weaponAmmoPrimary || l.weaponAmmoSecondary != s.weaponAmmoSecondary) {
    std::cerr << "[test_savegame] weaponAmmo mismatch\n";
    ++fails;
  }

  // Cargo array.
  if (!nearly(l.cargo[static_cast<std::size_t>(econ::CommodityId::Metals)],
              s.cargo[static_cast<std::size_t>(econ::CommodityId::Metals)])) {
    std::cerr << "[test_savegame] cargo[Metals] mismatch\n";
    ++fails;
  }

  // Missions.
  if (l.missions.size() != s.missions.size()) {
    std::cerr << "[test_savegame] missions size mismatch\n";
    ++fails;
  } else {
    const auto& a = l.missions.front();
    const auto& b = s.missions.front();
    if (a.id != b.id || a.type != b.type || a.factionId != b.factionId || a.toStation != b.toStation) {
      std::cerr << "[test_savegame] mission fields mismatch\n";
      ++fails;
    }
    if (a.leg != b.leg || a.scanned != b.scanned) {
      std::cerr << "[test_savegame] mission leg/scanned mismatch\n";
      ++fails;
    }
  }

  // Escort convoy persistence.
  if (l.escortConvoys.size() != s.escortConvoys.size()) {
    std::cerr << "[test_savegame] escortConvoys size mismatch\n";
    ++fails;
  } else if (!l.escortConvoys.empty()) {
    const auto& a = l.escortConvoys.front();
    const auto& b = s.escortConvoys.front();
    if (a.convoyId != b.convoyId || a.missionId != b.missionId || a.systemId != b.systemId) {
      std::cerr << "[test_savegame] escortConvoys id fields mismatch\n";
      ++fails;
    }
    if (!nearly(a.posKm.x, b.posKm.x) || !nearly(a.velKmS.y, b.velKmS.y)) {
      std::cerr << "[test_savegame] escortConvoys kinematics mismatch\n";
      ++fails;
    }
    if (!nearly(a.hullFrac, b.hullFrac) || !nearly(a.shieldFrac, b.shieldFrac)) {
      std::cerr << "[test_savegame] escortConvoys health mismatch\n";
      ++fails;
    }
    if (a.ambushSpawned != b.ambushSpawned || !nearly(a.nextAmbushDays, b.nextAmbushDays)) {
      std::cerr << "[test_savegame] escortConvoys ambush state mismatch\n";
      ++fails;
    }
  }
  // Bounty target persistence.
  if (l.bountyTargets.size() != s.bountyTargets.size()) {
    std::cerr << "[test_savegame] bountyTargets size mismatch\n";
    ++fails;
  } else if (!l.bountyTargets.empty()) {
    const auto& a = l.bountyTargets.front();
    const auto& b = s.bountyTargets.front();
    if (a.targetId != b.targetId || a.missionId != b.missionId || a.systemId != b.systemId || a.hideoutStation != b.hideoutStation) {
      std::cerr << "[test_savegame] bountyTargets id fields mismatch\n";
      ++fails;
    }
    if (!nearly(a.posKm.x, b.posKm.x) || !nearly(a.velKmS.z, b.velKmS.z)) {
      std::cerr << "[test_savegame] bountyTargets kinematics mismatch\n";
      ++fails;
    }
    if (!nearly(a.hullFrac, b.hullFrac) || !nearly(a.shieldFrac, b.shieldFrac)) {
      std::cerr << "[test_savegame] bountyTargets health mismatch\n";
      ++fails;
    }
  }


  if (l.trackedMissionId != s.trackedMissionId) {
    std::cerr << "[test_savegame] trackedMissionId mismatch\n";
    ++fails;
  }

  // World state persistence.
  if (l.resolvedSignalIds.size() != s.resolvedSignalIds.size()) {
    std::cerr << "[test_savegame] resolvedSignalIds size mismatch\n";
    ++fails;
  } else if (l.resolvedSignalIds != s.resolvedSignalIds) {
    std::cerr << "[test_savegame] resolvedSignalIds values mismatch\n";
    ++fails;
  }
  if (l.asteroidStates.size() != s.asteroidStates.size()) {
    std::cerr << "[test_savegame] asteroidStates size mismatch\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < l.asteroidStates.size(); ++i) {
      if (l.asteroidStates[i].asteroidId != s.asteroidStates[i].asteroidId ||
          !nearly(l.asteroidStates[i].remainingUnits, s.asteroidStates[i].remainingUnits)) {
        std::cerr << "[test_savegame] asteroidStates entry mismatch\n";
        ++fails;
        break;
      }
    }
  }

  // Logbook persistence.
  if (l.logbook.size() != s.logbook.size()) {
    std::cerr << "[test_savegame] logbook size mismatch\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < l.logbook.size(); ++i) {
      const auto& a = l.logbook[i];
      const auto& b = s.logbook[i];
      if (a.key != b.key || a.kind != b.kind || a.systemId != b.systemId || a.objectId != b.objectId) {
        std::cerr << "[test_savegame] logbook entry mismatch\n";
        ++fails;
        break;
      }
      if (a.commodity != b.commodity || !nearly(a.units, b.units) || !nearly(a.discoveredDay, b.discoveredDay) || a.sold != b.sold) {
        std::cerr << "[test_savegame] logbook entry fields mismatch\n";
        ++fails;
        break;
      }
    }
  }

  // Mission board persistence.
  if (l.missionOffersStationId != s.missionOffersStationId || l.missionOffersDayStamp != s.missionOffersDayStamp) {
    std::cerr << "[test_savegame] mission offers header mismatch\n";
    ++fails;
  }
  if (l.missionOffers.size() != s.missionOffers.size()) {
    std::cerr << "[test_savegame] mission offers size mismatch\n";
    ++fails;
  }

  // Traffic stamps.
  if (l.trafficStamps.size() != s.trafficStamps.size()) {
    std::cerr << "[test_savegame] trafficStamps size mismatch\n";
    ++fails;
  } else if (!l.trafficStamps.empty()) {
    if (l.trafficStamps.front().systemId != s.trafficStamps.front().systemId ||
        l.trafficStamps.front().dayStamp != s.trafficStamps.front().dayStamp) {
      std::cerr << "[test_savegame] trafficStamps entry mismatch\n";
      ++fails;
    }
  }

  // Traffic shipment replay.
  if (l.trafficShipments.size() != s.trafficShipments.size()) {
    std::cerr << "[test_savegame] trafficShipments size mismatch\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < l.trafficShipments.size(); ++i) {
      const auto& a = l.trafficShipments[i];
      const auto& b = s.trafficShipments[i];
      if (a.id != b.id || a.systemId != b.systemId || a.dayStamp != b.dayStamp || a.fromStation != b.fromStation || a.toStation != b.toStation) {
        std::cerr << "[test_savegame] trafficShipments entry mismatch\n";
        ++fails;
        break;
      }
      if (a.commodity != b.commodity || !nearly(a.units, b.units) || !nearly(a.departDay, b.departDay) || !nearly(a.arriveDay, b.arriveDay)) {
        std::cerr << "[test_savegame] trafficShipments fields mismatch\n";
        ++fails;
        break;
      }
    }
  }

  // Traffic convoy interdiction persistence.
  if (l.trafficInterdictions.size() != s.trafficInterdictions.size()) {
    std::cerr << "[test_savegame] trafficInterdictions size mismatch\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < l.trafficInterdictions.size(); ++i) {
      const auto& a = l.trafficInterdictions[i];
      const auto& b = s.trafficInterdictions[i];
      if (a.convoyId != b.convoyId || a.systemId != b.systemId || a.fromStation != b.fromStation || a.toStation != b.toStation) {
        std::cerr << "[test_savegame] trafficInterdictions entry mismatch\n";
        ++fails;
        break;
      }
      if (a.commodity != b.commodity || !nearly(a.units, b.units) || !nearly(a.expireDay, b.expireDay)) {
        std::cerr << "[test_savegame] trafficInterdictions fields mismatch\n";
        ++fails;
        break;
      }
    }
  }

  // Ambient traffic escort contract persistence.
  if (l.trafficEscort.active != s.trafficEscort.active ||
      l.trafficEscort.convoyId != s.trafficEscort.convoyId ||
      l.trafficEscort.payerFactionId != s.trafficEscort.payerFactionId ||
      l.trafficEscort.toStationId != s.trafficEscort.toStationId) {
    std::cerr << "[test_savegame] trafficEscort header mismatch\n";
    ++fails;
  } else {
    if (!nearly(l.trafficEscort.startDays, s.trafficEscort.startDays) ||
        !nearly(l.trafficEscort.untilDays, s.trafficEscort.untilDays) ||
        !nearly(l.trafficEscort.maxRangeKm, s.trafficEscort.maxRangeKm) ||
        !nearly(l.trafficEscort.tooFarSec, s.trafficEscort.tooFarSec) ||
        !nearly(l.trafficEscort.rewardCr, s.trafficEscort.rewardCr) ||
        !nearly(l.trafficEscort.bonusPerPirateCr, s.trafficEscort.bonusPerPirateCr) ||
        !nearly(l.trafficEscort.repReward, s.trafficEscort.repReward) ||
        l.trafficEscort.piratesKilled != s.trafficEscort.piratesKilled ||
        l.trafficEscort.piratesPresentAtStart != s.trafficEscort.piratesPresentAtStart) {
      std::cerr << "[test_savegame] trafficEscort fields mismatch\n";
      ++fails;
    }
  }

  if (l.trafficEscortSettlements.size() != s.trafficEscortSettlements.size()) {
    std::cerr << "[test_savegame] trafficEscortSettlements size mismatch\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < l.trafficEscortSettlements.size(); ++i) {
      const auto& a = l.trafficEscortSettlements[i];
      const auto& b = s.trafficEscortSettlements[i];
      if (a.convoyId != b.convoyId || !nearly(a.settledDay, b.settledDay)) {
        std::cerr << "[test_savegame] trafficEscortSettlements entry mismatch\n";
        ++fails;
        break;
      }
    }
  }

  // System security delta persistence.
  if (l.systemSecurityDeltas.size() != s.systemSecurityDeltas.size()) {
    std::cerr << "[test_savegame] systemSecurityDeltas size mismatch\n";
    ++fails;
  } else {
    for (std::size_t i = 0; i < l.systemSecurityDeltas.size(); ++i) {
      const auto& a = l.systemSecurityDeltas[i];
      const auto& b = s.systemSecurityDeltas[i];
      if (a.systemId != b.systemId) {
        std::cerr << "[test_savegame] systemSecurityDeltas systemId mismatch\n";
        ++fails;
        break;
      }
      if (!nearly(a.securityDelta, b.securityDelta) || !nearly(a.piracyDelta, b.piracyDelta) || !nearly(a.trafficDelta, b.trafficDelta)) {
        std::cerr << "[test_savegame] systemSecurityDeltas value mismatch\n";
        ++fails;
        break;
      }
      if (!nearly(a.lastUpdateDay, b.lastUpdateDay)) {
        std::cerr << "[test_savegame] systemSecurityDeltas lastUpdateDay mismatch\n";
        ++fails;
        break;
      }
    }
  }

  // Station overrides.
  if (l.stationOverrides.size() != s.stationOverrides.size()) {
    std::cerr << "[test_savegame] stationOverrides size mismatch\n";
    ++fails;
  } else {
    const auto& a = l.stationOverrides.front();
    const auto& b = s.stationOverrides.front();
    if (a.stationId != b.stationId) {
      std::cerr << "[test_savegame] stationId mismatch\n";
      ++fails;
    }
    if (!nearly(a.state.inventory[0], b.state.inventory[0])) {
      std::cerr << "[test_savegame] override inventory mismatch\n";
      ++fails;
    }
    if (a.state.history[0].empty() || b.state.history[0].empty()) {
      std::cerr << "[test_savegame] override history missing\n";
      ++fails;
    }
  }

  // Robustness: tolerate stale/corrupt section counts without breaking subsequent key parsing.
  {
    const std::string original = readAllText(path);
    if (original.empty()) {
      std::cerr << "[test_savegame] failed to read back saved file for robustness tests\n";
      ++fails;
    } else {
      // If the navRoute count is wrong, the loader should stop when it no longer sees "nav"
      // and continue parsing later keys (e.g. shipHull / trackedMissionId).
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "navRoute", s.navRoute.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_navroute.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt navRoute save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale navRoute count\n";
            ++fails;
          } else {
            if (x.navRoute.size() != s.navRoute.size()) {
              std::cerr << "[test_savegame] stale navRoute count changed parsed navRoute count\n";
              ++fails;
            }
            if (x.trackedMissionId != s.trackedMissionId) {
              std::cerr << "[test_savegame] stale navRoute count broke trackedMissionId parsing\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // If the missions count is wrong, the loader should stop when it no longer sees "mission"
      // and continue parsing top-level keys like trackedMissionId.
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "missions", s.missions.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_missions.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt missions save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale missions count\n";
            ++fails;
          } else {
            if (x.missions.size() != s.missions.size()) {
              std::cerr << "[test_savegame] stale missions count changed parsed mission count\n";
              ++fails;
            }
            if (x.trackedMissionId != s.trackedMissionId) {
              std::cerr << "[test_savegame] stale missions count broke trackedMissionId parsing\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // If the resolved_signals count is wrong, the loader should stop when it no longer sees
      // "signal_resolved" and continue parsing later sections (e.g. missions/trackedMissionId).
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "resolved_signals", s.resolvedSignalIds.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_resolved_signals.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt resolved_signals save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale resolved_signals count\n";
            ++fails;
          } else {
            if (x.resolvedSignalIds.size() != s.resolvedSignalIds.size()) {
              std::cerr << "[test_savegame] stale resolved_signals count changed parsed signal count\n";
              ++fails;
            }
            if (x.missions.size() != s.missions.size()) {
              std::cerr << "[test_savegame] stale resolved_signals count broke parsing of later keys (missions)\n";
              ++fails;
            }
            if (x.trackedMissionId != s.trackedMissionId) {
              std::cerr << "[test_savegame] stale resolved_signals count broke trackedMissionId parsing\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // If the asteroid_states count is wrong, the loader should stop when it no longer sees
      // "asteroid" and continue parsing keys that follow.
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "asteroid_states", s.asteroidStates.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_asteroid_states.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt asteroid_states save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale asteroid_states count\n";
            ++fails;
          } else {
            if (x.asteroidStates.size() != s.asteroidStates.size()) {
              std::cerr << "[test_savegame] stale asteroid_states count changed parsed asteroidStates count\n";
              ++fails;
            }
            if (x.missionOffers.size() != s.missionOffers.size()) {
              std::cerr << "[test_savegame] stale asteroid_states count broke parsing of later keys (mission_offers)\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // If the mission_offers count is wrong, the loader should stop when it no longer sees "offer"
      // and continue parsing keys that follow (e.g. reputation).
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "mission_offers", s.missionOffers.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_offers.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt mission_offers save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale mission_offers count\n";
            ++fails;
          } else {
            if (x.missionOffers.size() != s.missionOffers.size()) {
              std::cerr << "[test_savegame] stale mission_offers count changed parsed offers count\n";
              ++fails;
            }
            if (x.reputation.size() != s.reputation.size()) {
              std::cerr << "[test_savegame] stale mission_offers count broke parsing of later keys (reputation)\n";
              ++fails;
            } else if (!x.reputation.empty()) {
              if (!nearly(x.reputation.front().rep, s.reputation.front().rep)) {
                std::cerr << "[test_savegame] reputation value mismatch after stale mission_offers count\n";
                ++fails;
              }
            }
          }
          std::filesystem::remove(p);
        }
      }

      // If the trafficShipments count is wrong, the loader should stop when it no longer sees
      // "shipment" and continue parsing later keys (e.g. station_overrides).
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "trafficShipments", s.trafficShipments.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_traffic_shipments.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt trafficShipments save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale trafficShipments count\n";
            ++fails;
          } else {
            if (x.trafficShipments.size() != s.trafficShipments.size()) {
              std::cerr << "[test_savegame] stale trafficShipments count changed parsed shipment count\n";
              ++fails;
            }
            if (x.stationOverrides.size() != s.stationOverrides.size()) {
              std::cerr << "[test_savegame] stale trafficShipments count broke parsing of later keys (station_overrides)\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // Defensive: if trafficShipments lines are duplicated (same shipment id), the loader should
      // de-dup them instead of producing duplicated convoy signals after load.
      {
        bool duplicated = false;
        const std::string corrupt = duplicateSecondLineWithPrefix(original, "shipment ", duplicated);
        const std::string p = "savegame_test_corrupt_traffic_shipments_dup.sav";
        if (!duplicated || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write duplicate trafficShipments save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on duplicate trafficShipments entries\n";
            ++fails;
          } else {
            if (x.trafficShipments.size() != s.trafficShipments.size() - 1) {
              std::cerr << "[test_savegame] duplicate trafficShipments lines were not deduped\n";
              ++fails;
            }
            if (x.stationOverrides.size() != s.stationOverrides.size()) {
              std::cerr << "[test_savegame] duplicate trafficShipments broke parsing of later keys (station_overrides)\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // If the trafficInterdictions count is wrong, the loader should stop when it no longer sees
      // "interdict" and continue parsing later keys (e.g. station_overrides).
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "trafficInterdictions", s.trafficInterdictions.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_traffic_interdictions.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt trafficInterdictions save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale trafficInterdictions count\n";
            ++fails;
          } else {
            if (x.trafficInterdictions.size() != s.trafficInterdictions.size()) {
              std::cerr << "[test_savegame] stale trafficInterdictions count changed parsed entry count\n";
              ++fails;
            }
            if (x.stationOverrides.size() != s.stationOverrides.size()) {
              std::cerr << "[test_savegame] stale trafficInterdictions count broke parsing of later keys (station_overrides)\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // If the systemSecurityDeltas count is wrong, the loader should stop when it no longer sees
      // "secDelta" and continue parsing later keys (e.g. station_overrides).
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "systemSecurityDeltas", s.systemSecurityDeltas.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_system_security_deltas.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt systemSecurityDeltas save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale systemSecurityDeltas count\n";
            ++fails;
          } else {
            if (x.systemSecurityDeltas.size() != s.systemSecurityDeltas.size()) {
              std::cerr << "[test_savegame] stale systemSecurityDeltas count changed parsed entry count\n";
              ++fails;
            }
            if (x.stationOverrides.size() != s.stationOverrides.size()) {
              std::cerr << "[test_savegame] stale systemSecurityDeltas count broke parsing of later keys (station_overrides)\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // Defensive: if systemSecurityDeltas lines are duplicated (same system id), the loader should
      // de-dup them.
      {
        bool duplicated = false;
        const std::string corrupt = duplicateSecondLineWithPrefix(original, "secDelta ", duplicated);
        const std::string p = "savegame_test_corrupt_system_security_deltas_dup.sav";
        if (!duplicated || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write duplicate systemSecurityDeltas save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on duplicate systemSecurityDeltas entries\n";
            ++fails;
          } else {
            if (x.systemSecurityDeltas.size() != s.systemSecurityDeltas.size() - 1) {
              std::cerr << "[test_savegame] duplicate systemSecurityDeltas lines were not deduped\n";
              ++fails;
            }
            if (x.stationOverrides.size() != s.stationOverrides.size()) {
              std::cerr << "[test_savegame] duplicate systemSecurityDeltas broke parsing of later keys (station_overrides)\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // Defensive: if trafficInterdictions lines are duplicated (same convoy id), the loader should
      // de-dup them.
      {
        bool duplicated = false;
        const std::string corrupt = duplicateSecondLineWithPrefix(original, "interdict ", duplicated);
        const std::string p = "savegame_test_corrupt_traffic_interdictions_dup.sav";
        if (!duplicated || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write duplicate trafficInterdictions save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on duplicate trafficInterdictions entries\n";
            ++fails;
          } else {
            if (x.trafficInterdictions.size() != s.trafficInterdictions.size() - 1) {
              std::cerr << "[test_savegame] duplicate trafficInterdictions lines were not deduped\n";
              ++fails;
            }
            if (x.stationOverrides.size() != s.stationOverrides.size()) {
              std::cerr << "[test_savegame] duplicate trafficInterdictions broke parsing of later keys (station_overrides)\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // If station_overrides count is wrong (e.g., EOF before reaching the claimed count),
      // the loader should still succeed.
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "station_overrides", s.stationOverrides.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_overrides.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt station_overrides save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale station_overrides count\n";
            ++fails;
          } else if (x.stationOverrides.size() != s.stationOverrides.size()) {
            std::cerr << "[test_savegame] stale station_overrides count changed parsed overrides count\n";
            ++fails;
          }
          std::filesystem::remove(p);
        }
      }

      // If escort_convoys count is wrong, the loader should stop when it no longer sees
      // "convoy" and continue parsing later keys (e.g. trackedMissionId).
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "escort_convoys", s.escortConvoys.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_escort_convoys.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt escort_convoys save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale escort_convoys count\n";
            ++fails;
          } else {
            if (x.escortConvoys.size() != s.escortConvoys.size()) {
              std::cerr << "[test_savegame] stale escort_convoys count changed parsed escortConvoys count\n";
              ++fails;
            }
            if (x.trackedMissionId != s.trackedMissionId) {
              std::cerr << "[test_savegame] stale escort_convoys count broke trackedMissionId parsing\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

      // If bounty_targets count is wrong, the loader should stop when it no longer sees
      // "bounty_target" and continue parsing later keys (e.g. trackedMissionId).
      {
        bool replaced = false;
        const std::string corrupt = replaceCountLine(original, "bounty_targets", s.bountyTargets.size() + 5, replaced);
        const std::string p = "savegame_test_corrupt_bounty_targets.sav";
        if (!replaced || !writeAllText(p, corrupt)) {
          std::cerr << "[test_savegame] failed to write corrupt bounty_targets save\n";
          ++fails;
        } else {
          SaveGame x{};
          if (!loadFromFile(p, x)) {
            std::cerr << "[test_savegame] loadFromFile failed on stale bounty_targets count\n";
            ++fails;
          } else {
            if (x.bountyTargets.size() != s.bountyTargets.size()) {
              std::cerr << "[test_savegame] stale bounty_targets count changed parsed bountyTargets count\n";
              ++fails;
            }
            if (x.trackedMissionId != s.trackedMissionId) {
              std::cerr << "[test_savegame] stale bounty_targets count broke trackedMissionId parsing\n";
              ++fails;
            }
          }
          std::filesystem::remove(p);
        }
      }

    }
  }

  std::filesystem::remove(path);

  if (fails == 0) std::cout << "[test_savegame] pass\n";
  return fails;
}
