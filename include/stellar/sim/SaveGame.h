#pragma once

#include "stellar/econ/Economy.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h"
#include "stellar/sim/Industry.h"
#include "stellar/sim/Logbook.h"
#include "stellar/sim/TrafficEscort.h"
#include "stellar/sim/SystemSecurityDynamics.h"

#include <array>
#include <string>
#include <vector>

namespace stellar::sim {

struct StationEconomyOverride {
  StationId stationId{0};
  econ::StationEconomyState state{};
};

// Simple reputation record (per faction id).
// Stored in the save file so mission progression / market fees can persist.
struct FactionReputation {
  core::u32 factionId{0};
  double rep{0.0}; // roughly [-100,+100]
};

// Simple bounty record (per faction id).
// If non-zero, local police will treat the player as WANTED in that faction's space.
struct FactionBounty {
  core::u32 factionId{0};
  double bountyCr{0.0};
};

// Persistent fine record (per faction id).
// Fines are issued for minor offenses and can be paid later at stations.
// If a fine becomes overdue, it may convert into a bounty.
struct FactionFine {
  core::u32 factionId{0};
  double fineCr{0.0};
  double dueDay{0.0};
};

// Persistent day-stamp for deterministic "background" traffic simulation.
// Only systems the player has visited need to be tracked.
struct SystemTrafficStamp {
  SystemId systemId{0};
  int dayStamp{-1};
};

// Persisted record of a recent NPC-trade shipment.
//
// This is a small, save-friendly subset of sim::TrafficShipment used by
// stellar_game to restore the in-memory TrafficLedger after load.
//
// Design intent:
//  - Keep recent shipments only (TrafficLedger already prunes to a small window).
//  - Preserve stable IDs so TrafficConvoy signals remain consistent across save/load.
struct TrafficShipmentState {
  core::u64 id{0};

  SystemId systemId{0};
  int dayStamp{0};

  StationId fromStation{0};
  StationId toStation{0};
  core::u32 factionId{0};

  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  // Schedule metadata (used for convoy replay / visualization).
  double departDay{0.0};
  double arriveDay{0.0};
  double distKm{0.0};
  double speedKmS{0.0};
};

// Persisted record of a disrupted / interdicted traffic convoy.
//
// When a traffic convoy is destroyed in the simulation, stellar_game can record
// the convoy id here to prevent re-spawning it on system re-entry (anti-farming)
// and to allow markets to reflect lost cargo.
//
// This is intentionally small and only needs to cover a short window (similar to
// TrafficLedger keepDays).
struct TrafficInterdictionState {
  core::u64 convoyId{0};
  SystemId systemId{0};

  StationId fromStation{0};
  StationId toStation{0};

  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  // When this convoy would have naturally expired (usually arriveDay). Used for pruning.
  double expireDay{0.0};
};

// Player-owned cargo stored at a specific station.
//
// This is intentionally lightweight (no pointers/strings) and persists across runs.
// Storage fees are accrued lazily using lastFeeDay -> current timeDays.
struct StationStorage {
  StationId stationId{0};
  econ::StationType stationType{econ::StationType::Outpost};
  core::u32 factionId{0};

  // Snapshot of station tariff at the time the storage record was created.
  // Used by the fee model (see sim/Warehouse).
  double feeRate{0.0};

  // Storage fees are accrued from lastFeeDay up to the current timeDays when the
  // player interacts with the warehouse.
  double lastFeeDay{0.0};
  double feesDueCr{0.0};

  // Stored commodity units.
  std::array<double, econ::kCommodityCount> cargo{};
};

// Persisted asteroid depletion state (for deterministic in-system resource fields).
//
// The game can generate asteroid nodes deterministically per system and then
// override the remainingUnits from this table so mining progress persists across
// system transitions and save/load.
struct AsteroidState {
  core::u64 asteroidId{0};
  double remainingUnits{0.0};
};

// Persisted state for an escort-mission convoy NPC.
//
// Escort missions spawn a specific convoy ship (mission.targetNpcId). Without persisting
// minimal state, a save/load mid-escort would recreate the world without the convoy and
// the mission would immediately fail.
//
// This is intentionally small: just enough to respawn the convoy near its last known
// state and keep the escort mission's one-shot ambush scheduling stable.
struct EscortConvoyState {
  core::u64 convoyId{0};
  core::u64 missionId{0};

  SystemId systemId{0};
  StationId fromStation{0};
  StationId toStation{0};

  // Physics state (km / km/s).
  math::Vec3d posKm{0,0,0};
  math::Vec3d velKmS{0,0,0};
  math::Quatd orient{1,0,0,0};
  math::Vec3d angVelRadS{0,0,0};

  // Combat state (fractions of max; max values are derived from the NPC loadout).
  double hullFrac{1.0};
  double shieldFrac{1.0};

  // Best-effort haul value used for ambush scaling / UI.
  double cargoValueCr{0.0};

  // Escort mission runtime state (kept small so missions behave consistently across save/load).
  double tooFarSec{0.0};
  bool ambushSpawned{false};
  double nextAmbushDays{0.0};
};


// Mission-critical bounty targets (for bounty scan/kill missions) are persisted so they don't
// reset when the player saves/loads or leaves/returns to the system.
struct BountyTargetState {
  core::u64 targetId{0};
  core::u64 missionId{0};

  SystemId systemId{0};
  StationId hideoutStation{0}; // station the target tends to linger near (best-effort)

  // Physics state (km / km/s).
  math::Vec3d posKm{0,0,0};
  math::Vec3d velKmS{0,0,0};
  math::Quatd orient{1,0,0,0};
  math::Vec3d angVelRadS{0,0,0};

  // Combat state (fractions of max; max values are derived from the NPC loadout).
  double hullFrac{1.0};
  double shieldFrac{1.0};
};


// Lightweight "gameplay" mission representation.
// Stored in the save file so early progression loops (cargo delivery/courier/bounties)
// persist across runs.
// IMPORTANT: MissionType is persisted in save files as an integer.
//
// Always use explicit, stable numeric values here.
// Do not reorder existing members, or older saves can silently remap mission types.
enum class MissionType : core::u8 {
  Courier       = 0,
  Delivery      = 1,
  BountyScan    = 2,
  BountyKill    = 3,
  MultiDelivery = 4,
  Passenger     = 5,
  Smuggle       = 6,

  // Recover goods from a mission derelict signal and return to a station.
  Salvage       = 7,

  // Escort a convoy from one station/system to another.
  Escort        = 8,
};

struct Mission {
  core::u64 id{0};
  MissionType type{MissionType::Courier};

  // Issuing faction (used for reputation effects).
  core::u32 factionId{0};

  SystemId fromSystem{0};
  StationId fromStation{0};

  SystemId toSystem{0};
  StationId toStation{0};

  // Delivery / smuggling missions use a commodity + units.
  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  // Multi-hop deliveries: optional intermediate stop.
  // If viaSystem/viaStation are non-zero, the player must first deliver to the via stop,
  // then proceed to the final destination (toSystem/toStation).
  SystemId viaSystem{0};
  StationId viaStation{0};
  core::u8 leg{0}; // 0 = heading to via (if any) else final; 1 = heading to final

  // Bounty scan missions can attach to a specific NPC id (spawned in the target system).
  // 0 means "no specific target".
  core::u64 targetNpcId{0};

  // Bounty scan progress flag.
  bool scanned{false};

  double reward{0.0};
  double deadlineDay{0.0};

  bool completed{false};
  bool failed{false};

  // If true, the station provided the cargo at acceptance time (taking from station inventory).
  // If false, the player must source the cargo themselves.
  bool cargoProvided{true};
};

struct SaveGame {
  int version{29};

  core::u64 seed{0};
  double timeDays{0.0};

  SystemId currentSystem{0};
  StationId dockedStation{0};

  // Player ship
  math::Vec3d shipPosKm{0,0,0};
  math::Vec3d shipVelKmS{0,0,0};
  math::Quatd shipOrient{1,0,0,0};
  math::Vec3d shipAngVelRadS{0,0,0};

  // Economy
  double credits{1000.0};
  // Outstanding insurance loan debt (paid down manually at stations).
  double insuranceDebtCr{0.0};
  std::array<double, econ::kCommodityCount> cargo{}; // units

  // Exploration
  double explorationDataCr{0.0};
  std::vector<core::u64> scannedKeys{};
  std::vector<LogbookEntry> logbook{};

  // World state (in-system signals / resource depletion)
  //
  // These are intentionally lightweight: only IDs (signals) and per-asteroid
  // remaining units (resource fields). The game uses deterministic IDs for
  // procedurally generated objects so it can safely persist partial depletion.
  std::vector<core::u64> resolvedSignalIds{};
  std::vector<AsteroidState> asteroidStates{};

  // Ship meta/progression
  double fuel{45.0};
  double fuelMax{45.0};
  double fsdRangeLy{18.0};
  double hull{1.0}; // 0..1
  double shield{1.0}; // 0..1
  double heat{0.0}; // gameplay heat (0..~120)

  // Countermeasures (ammo / consumables)
  //
  // These are intentionally simple counts of "bursts" available.
  // Gameplay code clamps these to the per-hull maximum on load.
  core::u8 cmFlares{6};
  core::u8 cmChaff{4};
  core::u8 cmHeatSinks{1};
  // Power distributor state
  int pipsEng{2};
  int pipsWep{2};
  int pipsSys{2};
  double capEngFrac{1.0};
  double capWepFrac{1.0};
  double capSysFrac{1.0};

  double cargoCapacityKg{420.0};

  // Passenger capacity (seats) for cabin-style missions.
  // Kept intentionally simple for the prototype loop.
  int passengerSeats{2};
  double fsdReadyDay{0.0}; // timeDays when the next hyperspace jump is allowed

  // Navigation UI state (quality-of-life).
  //
  // This is persisted so plotted routes + auto-run preferences survive save/load.
  std::vector<SystemId> navRoute{};
  core::u32 navRouteHop{0};
  bool navAutoRun{false};
  core::u8 navRouteMode{0};
  bool navConstrainToCurrentFuelRange{true};
  StationId pendingArrivalStation{0};


  // Loadout / progression (kept simple for now: small ints, interpreted by gameplay code).
  // These are *not* physics-critical; they tune HUD/combat feel and basic progression loops.
  core::u8 shipHull{0};        // 0 = Scout (starter), 1 = Hauler, 2 = Fighter
  core::u8 thrusterMk{1};      // 1..3
  core::u8 shieldMk{1};        // 1..3
  core::u8 distributorMk{1};   // 1..3
  core::u8 weaponPrimary{0};   // WeaponType enum index (see ShipLoadout.h)
  core::u8 weaponSecondary{2}; // default cannon

  // Guided weapons consume ammo; 255 is treated as "full" (backwards compat for older saves).
  core::u8 weaponAmmoPrimary{255};
  core::u8 weaponAmmoSecondary{255};

  // Smuggling / stealth: reduces chance of cargo scans when carrying contraband.
  core::u8 smuggleHoldMk{0}; // 0..3

  // Fuel scoop module (0 = none, 1..3 = Mk1..Mk3)
  core::u8 fuelScoopMk{0};

  // Missions
  core::u64 nextMissionId{1};
  std::vector<Mission> missions{};

  // Escort mission persistence (mission-critical NPCs).
  std::vector<EscortConvoyState> escortConvoys{};

  // Mission-critical bounty targets (bounty scan/kill missions).
  std::vector<BountyTargetState> bountyTargets{};

  // Industry / fabrication orders (station processing queues).
  core::u64 nextIndustryOrderId{1};
  std::vector<IndustryOrder> industryOrders{};

  // Mission tracker UI (quality-of-life): remember which active mission is "tracked".
  core::u64 trackedMissionId{0};

  // Mission board (cached offers) - persisted so boards don't reroll on UI refresh / reload.
  StationId missionOffersStationId{0};
  int missionOffersDayStamp{-1};
  std::vector<Mission> missionOffers{};

  // Reputation
  std::vector<FactionReputation> reputation{};

  // Law / bounties
  std::vector<FactionBounty> bounties{};

  // Outstanding fines (minor offenses), payable later.
  std::vector<FactionFine> fines{};
  // Bounty vouchers earned for destroying criminals (redeem at stations).
  std::vector<FactionBounty> bountyVouchers{};

  // Background NPC traffic simulation (for markets). Stores last simulated day per system.
  std::vector<SystemTrafficStamp> trafficStamps{};

  // Recent NPC trade shipments captured from sim::simulateNpcTradeTraffic(...).
  //
  // This is used to restore the in-memory TrafficLedger on load so that
  // TrafficConvoy signals remain consistent across save/load.
  std::vector<TrafficShipmentState> trafficShipments{};

  // Recently disrupted traffic convoys (anti-farm + economy impact persistence).
  std::vector<TrafficInterdictionState> trafficInterdictions{};

  // Ambient traffic escort contracts (accepted on-grid via the System Scanner).
  //
  // These were previously purely runtime state. Persisting them makes the
  // traffic layer feel more like a real, consistent game system:
  //  - Mid-contract quicksave/quickload no longer insta-fails the escort.
  //  - Completed convoys can't be repeatedly re-escorted after reload.
  TrafficEscortContractState trafficEscort{};
  std::vector<TrafficEscortSettlementState> trafficEscortSettlements{};

  // Per-system dynamic deltas applied to sim::SystemSecurityProfile.
  // These decay over time back toward the deterministic baseline.
  std::vector<SystemSecurityDeltaState> systemSecurityDeltas{};

  // Player station storage / warehouse entries.
  std::vector<StationStorage> stationStorage{};

  std::vector<StationEconomyOverride> stationOverrides{};
};

bool saveToFile(const SaveGame& s, const std::string& path);
bool loadFromFile(const std::string& path, SaveGame& out);

} // namespace stellar::sim
