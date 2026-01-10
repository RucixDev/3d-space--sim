#include "stellar/ui/AudioSettings.h"

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

std::string defaultAudioSettingsPath() {
  return "audio_settings.txt";
}

AudioSettings makeDefaultAudioSettings() {
  return AudioSettings{};
}

static float clamp01(float v) {
  return std::clamp(v, 0.0f, 1.0f);
}

bool saveAudioSettingsToFile(const AudioSettings& s, const std::string& path) {
  std::ofstream f(path);
  if (!f) {
    core::log(core::LogLevel::Error, "AudioSettings: unable to write " + path);
    return false;
  }

  f << "enabled " << (s.enabled ? 1 : 0) << "\n";
  f << "master " << clamp01(s.master) << "\n";
  f << "engine " << clamp01(s.engine) << "\n";
  f << "weapons " << clamp01(s.weapons) << "\n";
  f << "ui " << clamp01(s.ui) << "\n";
  f << "world " << clamp01(s.world) << "\n";
  f << "spatialize " << (s.spatialize ? 1 : 0) << "\n";

  return true;
}

bool loadAudioSettingsFromFile(const std::string& path, AudioSettings& out) {
  std::ifstream f(path);
  if (!f) return false;

  AudioSettings s = out;

  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    std::string key;
    if (!(iss >> key)) continue;
    key = lowerAscii(key);

    if (key == "enabled") {
      int v = 1;
      iss >> v;
      s.enabled = (v != 0);
    } else if (key == "master") {
      float v = s.master;
      iss >> v;
      s.master = clamp01(v);
    } else if (key == "engine") {
      float v = s.engine;
      iss >> v;
      s.engine = clamp01(v);
    } else if (key == "weapons") {
      float v = s.weapons;
      iss >> v;
      s.weapons = clamp01(v);
    } else if (key == "ui") {
      float v = s.ui;
      iss >> v;
      s.ui = clamp01(v);
    } else if (key == "world") {
      float v = s.world;
      iss >> v;
      s.world = clamp01(v);
    } else if (key == "spatialize") {
      int v = 1;
      iss >> v;
      s.spatialize = (v != 0);
    } else {
      // Unknown key: ignore for forward compatibility.
      continue;
    }
  }

  out = s;
  return true;
}

} // namespace stellar::ui
