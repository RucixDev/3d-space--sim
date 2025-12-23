#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Market.h"

#include <filesystem>
#include <optional>
#include <string>

namespace stellar::sim {

struct ShipState {
  stellar::math::Vec3d positionAU{0.0, 0.0, 5.0};
  stellar::math::Vec3d velocityAUPerDay{0.0, 0.0, 0.0};
  double yawDeg = 0.0;
  double pitchDeg = 0.0;
};

struct SaveGame {
  int version = 1;

  stellar::core::u64 universeSeed = 1;
  std::size_t currentSystemIndex = 0;
  double simTimeDays = 0.0;

  ShipState ship;

  double credits = 1000.0;
  CargoHold cargo;

  std::size_t selectedCommodity = 0;
};

bool writeSave(const SaveGame& save, const std::filesystem::path& path, std::string* outError = nullptr);
std::optional<SaveGame> readSave(const std::filesystem::path& path, std::string* outError = nullptr);

} // namespace stellar::sim