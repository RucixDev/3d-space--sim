#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/sim/Celestial.h"

#include <vector>

namespace stellar::sim {

enum class MissionType : core::u8 {
  Courier = 0,
  Delivery,
  BountyScan,
  Count
};

// Minimal mission record suitable for savegames.
// UI-friendly strings should be generated at runtime.
struct Mission {
  core::u64 id{0};
  MissionType type{MissionType::Courier};

  SystemId originSystem{0};
  StationId originStation{0};

  SystemId destSystem{0};
  StationId destStation{0};

  // Delivery payload (used only for Delivery missions).
  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  // Bounty scan flag (used only for BountyScan missions).
  bool scanned{false};

  // Rewards / time
  double rewardCredits{0.0};
  double expiryDay{0.0}; // 0 = no expiry
};

} // namespace stellar::sim
