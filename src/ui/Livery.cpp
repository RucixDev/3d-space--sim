#include "stellar/ui/Livery.h"

#include "stellar/core/Log.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

namespace stellar::ui {

static std::string lowerAscii(std::string s) {
  for (char& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

std::string defaultLiveryPath() { return "livery.txt"; }

LiveryConfig makeDefaultLivery() { return LiveryConfig{}; }

bool saveLiveryToFile(const LiveryConfig& cfg, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Warn, "Livery: failed to open for writing: " + path);
    return false;
  }

  f.setf(std::ios::fixed);
  f.precision(6);

  f << "StellarForgeLivery " << cfg.version << "\n";
  f << "pattern " << stellar::render::toString(cfg.pattern) << "\n";
  f << "seed " << (unsigned long long)cfg.seed << "\n";
  f << "base " << cfg.base[0] << " " << cfg.base[1] << " " << cfg.base[2] << "\n";
  f << "accent1 " << cfg.accent1[0] << " " << cfg.accent1[1] << " " << cfg.accent1[2] << "\n";
  f << "accent2 " << cfg.accent2[0] << " " << cfg.accent2[1] << " " << cfg.accent2[2] << "\n";
  f << "scale " << cfg.scale << "\n";
  f << "angleDeg " << cfg.angleDeg << "\n";
  f << "detail " << cfg.detail << "\n";
  f << "wear " << cfg.wear << "\n";
  f << "contrast " << cfg.contrast << "\n";
  f << "decal " << (cfg.decal ? 1 : 0) << "\n";
  f << "textureSizePx " << cfg.textureSizePx << "\n";
  f << "applyInWorld " << (cfg.applyInWorld ? 1 : 0) << "\n";
  return true;
}

bool loadLiveryFromFile(const std::string& path, LiveryConfig& out) {
  std::ifstream f(path);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Debug, "Livery: file not found: " + path);
    return false;
  }

  std::string header;
  int version = 0;
  if (!(f >> header >> version)) {
    stellar::core::log(stellar::core::LogLevel::Warn, "Livery: failed to read header");
    return false;
  }
  if (header != "StellarForgeLivery") {
    stellar::core::log(stellar::core::LogLevel::Warn, "Livery: bad header");
    return false;
  }

  LiveryConfig cfg = makeDefaultLivery();
  cfg.version = version;

  std::string line;
  std::getline(f, line); // consume remainder of header line
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream ss(line);
    std::string key;
    if (!(ss >> key)) continue;
    key = lowerAscii(key);

    if (key == "pattern") {
      std::string v;
      ss >> v;
      if (auto p = stellar::render::patternFromString(v)) cfg.pattern = *p;
    } else if (key == "seed") {
      unsigned long long s = 0;
      ss >> s;
      cfg.seed = (core::u64)s;
    } else if (key == "base") {
      ss >> cfg.base[0] >> cfg.base[1] >> cfg.base[2];
    } else if (key == "accent1") {
      ss >> cfg.accent1[0] >> cfg.accent1[1] >> cfg.accent1[2];
    } else if (key == "accent2") {
      ss >> cfg.accent2[0] >> cfg.accent2[1] >> cfg.accent2[2];
    } else if (key == "scale") {
      ss >> cfg.scale;
    } else if (key == "angledeg") {
      ss >> cfg.angleDeg;
    } else if (key == "detail") {
      ss >> cfg.detail;
    } else if (key == "wear") {
      ss >> cfg.wear;
    } else if (key == "contrast") {
      ss >> cfg.contrast;
    } else if (key == "decal") {
      int v = 1;
      ss >> v;
      cfg.decal = (v != 0);
    } else if (key == "texturesizepx") {
      ss >> cfg.textureSizePx;
    } else if (key == "applyinworld") {
      int v = 1;
      ss >> v;
      cfg.applyInWorld = (v != 0);
    }
  }

  // Clamp to sane values.
  cfg.scale = std::clamp(cfg.scale, 0.25f, 4.0f);
  cfg.angleDeg = std::clamp(cfg.angleDeg, -180.0f, 180.0f);
  cfg.detail = std::clamp(cfg.detail, 0.0f, 1.0f);
  cfg.wear = std::clamp(cfg.wear, 0.0f, 1.0f);
  cfg.contrast = std::clamp(cfg.contrast, 0.0f, 1.0f);
  cfg.textureSizePx = std::clamp(cfg.textureSizePx, 64, 2048);

  out = cfg;
  return true;
}

} // namespace stellar::ui
