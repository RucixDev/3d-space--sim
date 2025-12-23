#include "stellar/sim/SaveGame.h"

#include "stellar/sim/Market.h"

#include <charconv>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string_view>
#include <system_error>
#include <type_traits>

namespace stellar::sim {
namespace {

std::string trim(std::string s) {
  auto isSpace = [](unsigned char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; };
  while (!s.empty() && isSpace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
  while (!s.empty() && isSpace(static_cast<unsigned char>(s.back()))) s.pop_back();
  return s;
}

template <typename T>
bool parseInt(std::string_view s, T& out) {
  static_assert(std::is_integral_v<T>);
  const char* first = s.data();
  const char* last = s.data() + s.size();
  auto [ptr, ec] = std::from_chars(first, last, out);
  return ec == std::errc{} && ptr == last;
}

bool parseDouble(std::string_view s, double& out) {
  std::string tmp(s);
  char* end = nullptr;
  out = std::strtod(tmp.c_str(), &end);
  return end && end == tmp.c_str() + tmp.size();
}

bool parseVec3(std::string_view s, stellar::math::Vec3d& out) {
  // Expect x,y,z
  const std::string tmp(s);
  std::istringstream iss(tmp);
  char comma1 = 0;
  char comma2 = 0;
  double x = 0.0, y = 0.0, z = 0.0;
  if (!(iss >> x)) return false;
  if (!(iss >> comma1) || comma1 != ',') return false;
  if (!(iss >> y)) return false;
  if (!(iss >> comma2) || comma2 != ',') return false;
  if (!(iss >> z)) return false;
  out = {x, y, z};
  return true;
}

std::string vec3ToString(const stellar::math::Vec3d& v) {
  std::ostringstream oss;
  oss << v.x << "," << v.y << "," << v.z;
  return oss.str();
}

}

bool writeSave(const SaveGame& save, const std::filesystem::path& path, std::string* outError) {
  std::ofstream out(path, std::ios::binary);
  if (!out) {
    if (outError) *outError = "Failed to open save file for writing: " + path.string();
    return false;
  }

  out << "version=" << save.version << "\n";
  out << "universe_seed=" << save.universeSeed << "\n";
  out << "system_index=" << save.currentSystemIndex << "\n";
  out << "sim_time_days=" << save.simTimeDays << "\n";
  out << "ship.pos_au=" << vec3ToString(save.ship.positionAU) << "\n";
  out << "ship.vel_au_per_day=" << vec3ToString(save.ship.velocityAUPerDay) << "\n";
  out << "ship.yaw_deg=" << save.ship.yawDeg << "\n";
  out << "ship.pitch_deg=" << save.ship.pitchDeg << "\n";
  out << "credits=" << save.credits << "\n";
  out << "selected_commodity=" << save.selectedCommodity << "\n";

  for (std::size_t i = 0; i < kCommodityCount; ++i) {
    const auto c = static_cast<Commodity>(i);
    out << "cargo." << commodityKey(c) << "=" << save.cargo.units[i] << "\n";
  }

  return true;
}

std::optional<SaveGame> readSave(const std::filesystem::path& path, std::string* outError) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    if (outError) *outError = "Failed to open save file for reading: " + path.string();
    return std::nullopt;
  }

  SaveGame save;

  std::string line;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;

    const auto eq = line.find('=');
    if (eq == std::string::npos) continue;

    const std::string key = trim(line.substr(0, eq));
    const std::string val = trim(line.substr(eq + 1));

    if (key == "version") {
      parseInt(val, save.version);
    } else if (key == "universe_seed") {
      parseInt(val, save.universeSeed);
    } else if (key == "system_index") {
      parseInt(val, save.currentSystemIndex);
    } else if (key == "sim_time_days") {
      parseDouble(val, save.simTimeDays);
    } else if (key == "ship.pos_au") {
      parseVec3(val, save.ship.positionAU);
    } else if (key == "ship.vel_au_per_day") {
      parseVec3(val, save.ship.velocityAUPerDay);
    } else if (key == "ship.yaw_deg") {
      parseDouble(val, save.ship.yawDeg);
    } else if (key == "ship.pitch_deg") {
      parseDouble(val, save.ship.pitchDeg);
    } else if (key == "credits") {
      parseDouble(val, save.credits);
    } else if (key == "selected_commodity") {
      parseInt(val, save.selectedCommodity);
    } else if (key.rfind("cargo.", 0) == 0) {
      const auto commodityStr = key.substr(6);
      for (std::size_t i = 0; i < kCommodityCount; ++i) {
        const auto c = static_cast<Commodity>(i);
        if (commodityStr == commodityKey(c)) {
          parseInt(val, save.cargo.units[i]);
          break;
        }
      }
    }
  }

  return save;
}

} // namespace stellar::sim
