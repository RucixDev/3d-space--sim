#pragma once

#include "stellar/sim/SaveGame.h"

namespace stellar::sim {

// Small helper utilities for UI/logging.
// (Keep this header light and inline-only.)

inline const char* missionTypeName(MissionType t) {
  switch (t) {
    case MissionType::Courier: return "Courier";
    case MissionType::Delivery: return "Delivery";
    case MissionType::MultiDelivery: return "MultiDelivery";
    case MissionType::Escort: return "Escort";
    case MissionType::Salvage: return "Salvage";
    case MissionType::Passenger: return "Passenger";
    case MissionType::Smuggle: return "Smuggle";
    case MissionType::BountyScan: return "BountyScan";
    case MissionType::BountyKill: return "BountyKill";
    default: return "Unknown";
  }
}

} // namespace stellar::sim
