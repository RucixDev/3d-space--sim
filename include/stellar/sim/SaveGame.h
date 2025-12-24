#pragma once

#include "stellar/econ/Economy.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Celestial.h"

#include <array>
#include <string>
#include <vector>

namespace stellar::sim {

struct StationEconomyOverride {
  StationId stationId{0};
  econ::StationEconomyState state{};
};

// Lightweight "gameplay" mission representation.
// Stored in the save file so early progression loops (cargo delivery/courier/bounties)
// persist across runs.
enum class MissionType : core::u8 {
  Courier = 0,
  Delivery,
  BountyScan,
};

struct Mission {
  core::u64 id{0};
  MissionType type{MissionType::Courier};

  SystemId fromSystem{0};
  StationId fromStation{0};

  SystemId toSystem{0};
  StationId toStation{0};

  // Delivery missions use a commodity + units.
  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  // Bounty scan missions can attach to a specific NPC id (spawned in the target system).
  // 0 means "no specific target".
  core::u64 targetNpcId{0};

  double reward{0.0};
  double deadlineDay{0.0};

  bool completed{false};
  bool failed{false};

  // If true, the station provided the cargo at acceptance time (taking from station inventory).
  // If false, the player must source the cargo themselves.
  bool cargoProvided{true};
};

struct SaveGame {
  int version{2};

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
  std::array<double, econ::kCommodityCount> cargo{}; // units

  // Ship meta/progression
  double fuel{45.0};
  double fuelMax{45.0};
  double hull{1.0}; // 0..1
  double cargoCapacityKg{420.0};
  double fsdReadyDay{0.0}; // timeDays when the next hyperspace jump is allowed

  // Missions
  core::u64 nextMissionId{1};
  std::vector<Mission> missions{};

  std::vector<StationEconomyOverride> stationOverrides{};
};

bool saveToFile(const SaveGame& s, const std::string& path);
bool loadFromFile(const std::string& path, SaveGame& out);

} // namespace stellar::sim
