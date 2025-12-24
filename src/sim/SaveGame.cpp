#include "stellar/sim/SaveGame.h"

#include "stellar/core/Log.h"

#include <fstream>
#include <sstream>

namespace stellar::sim {

bool saveToFile(const SaveGame& s, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Error, "SaveGame: failed to open file for writing: " + path);
    return false;
  }

  f.setf(std::ios::fixed);
  f.precision(8);

  f << "StellarForgeSave " << s.version << "\n";
  f << "seed " << s.seed << "\n";
  f << "timeDays " << s.timeDays << "\n";
  f << "currentSystem " << s.currentSystem << "\n";
  f << "dockedStation " << s.dockedStation << "\n";

  f << "shipPosKm " << s.shipPosKm.x << " " << s.shipPosKm.y << " " << s.shipPosKm.z << "\n";
  f << "shipVelKmS " << s.shipVelKmS.x << " " << s.shipVelKmS.y << " " << s.shipVelKmS.z << "\n";
  f << "shipOrient " << s.shipOrient.w << " " << s.shipOrient.x << " " << s.shipOrient.y << " " << s.shipOrient.z << "\n";
  f << "shipAngVel " << s.shipAngVelRadS.x << " " << s.shipAngVelRadS.y << " " << s.shipAngVelRadS.z << "\n";

  f << "credits " << s.credits << "\n";

  f << "cargo";
  for (double u : s.cargo) f << " " << u;
  f << "\n";

  // Player state / progression
  f << "hull " << s.hull << " " << s.hullMax << "\n";
  f << "fuel " << s.fuel << " " << s.fuelMax << "\n";
  f << "cargoCapacity " << s.cargoCapacity << "\n";
  f << "fsdCooldownSec " << s.fsdCooldownSec << "\n";

  // Missions
  f << "missions " << s.missions.size() << "\n";
  for (const auto& m : s.missions) {
    f << "mission "
      << m.id << " "
      << static_cast<int>(m.type) << " "
      << m.originSystem << " "
      << m.originStation << " "
      << m.destSystem << " "
      << m.destStation << " "
      << static_cast<int>(m.commodity) << " "
      << m.units << " "
      << (m.scanned ? 1 : 0) << " "
      << m.rewardCredits << " "
      << m.expiryDay
      << "\n";
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
    } else if (key == "cargo") {
      for (std::size_t i = 0; i < econ::kCommodityCount; ++i) f >> out.cargo[i];
    } else if (key == "hull") {
      f >> out.hull >> out.hullMax;
    } else if (key == "fuel") {
      f >> out.fuel >> out.fuelMax;
    } else if (key == "cargoCapacity") {
      f >> out.cargoCapacity;
    } else if (key == "fsdCooldownSec") {
      f >> out.fsdCooldownSec;
    } else if (key == "missions") {
      std::size_t n = 0;
      f >> n;
      out.missions.clear();
      out.missions.reserve(n);

      for (std::size_t i = 0; i < n; ++i) {
        if (!expectToken(f, "mission")) return false;
        Mission m{};
        int type = 0;
        int comm = 0;
        int scanned = 0;
        f >> m.id >> type >> m.originSystem >> m.originStation >> m.destSystem >> m.destStation;
        f >> comm >> m.units >> scanned >> m.rewardCredits >> m.expiryDay;
        m.type = static_cast<MissionType>(type);
        m.commodity = static_cast<econ::CommodityId>(comm);
        m.scanned = (scanned != 0);
        out.missions.push_back(std::move(m));
      }
    } else if (key == "station_overrides") {
      std::size_t n = 0;
      f >> n;
      out.stationOverrides.clear();
      out.stationOverrides.reserve(n);

      for (std::size_t si = 0; si < n; ++si) {
        StationEconomyOverride ov{};
        if (!expectToken(f, "station")) return false;
        f >> ov.stationId;

        // Parse station block until endstation
        std::string sk;
        while (f >> sk) {
          if (sk == "endstation") break;
          if (sk == "lastUpdateDay") {
            f >> ov.state.lastUpdateDay;
          } else if (sk == "lastSampleDay") {
            f >> ov.state.lastSampleDay;
          } else if (sk == "inventory") {
            for (std::size_t i = 0; i < econ::kCommodityCount; ++i) f >> ov.state.inventory[i];
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
