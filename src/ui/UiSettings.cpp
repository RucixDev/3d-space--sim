#include "stellar/ui/UiSettings.h"

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

const char* toString(UiTheme t) {
  switch (t) {
    case UiTheme::Dark: return "Dark";
    case UiTheme::Light: return "Light";
    case UiTheme::Classic: return "Classic";
    case UiTheme::HighContrast: return "HighContrast";
    default: return "Dark";
  }
}

std::optional<UiTheme> themeFromString(const std::string& s) {
  const std::string k = lowerAscii(s);
  if (k == "dark") return UiTheme::Dark;
  if (k == "light") return UiTheme::Light;
  if (k == "classic") return UiTheme::Classic;
  if (k == "highcontrast" || k == "high_contrast" || k == "high-contrast") return UiTheme::HighContrast;
  return std::nullopt;
}

std::string defaultUiSettingsPath() {
  return "ui_settings.txt";
}

UiSettings makeDefaultUiSettings() {
  return UiSettings{};
}

bool saveUiSettingsToFile(const UiSettings& s, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Warn, "UiSettings: failed to open for writing: " + path);
    return false;
  }

  f.setf(std::ios::fixed);
  f.precision(6);

  f << "StellarForgeUiSettings " << s.version << "\n";
  f << "imguiIniFile " << (s.imguiIniFile.empty() ? std::string("imgui.ini") : s.imguiIniFile) << "\n";
  f << "autoScaleFromDpi " << (s.autoScaleFromDpi ? 1 : 0) << "\n";
  f << "scaleUser " << s.scaleUser << "\n";
  f << "theme " << toString(s.theme) << "\n";
  f << "autoSaveOnExit " << (s.autoSaveOnExit ? 1 : 0) << "\n";

  f << "dockEnabled " << (s.dock.dockingEnabled ? 1 : 0) << "\n";
  f << "dockPassthruCentral " << (s.dock.passthruCentral ? 1 : 0) << "\n";
  f << "dockLockCentralView " << (s.dock.lockCentralView ? 1 : 0) << "\n";
  f << "dockLeftRatio " << s.dock.leftRatio << "\n";
  f << "dockRightRatio " << s.dock.rightRatio << "\n";
  f << "dockBottomRatio " << s.dock.bottomRatio << "\n";

  return true;
}

bool loadUiSettingsFromFile(const std::string& path, UiSettings& out) {
  std::ifstream f(path);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Debug, "UiSettings: file not found: " + path);
    return false;
  }

  std::string header;
  int version = 0;
  if (!(f >> header >> version)) {
    stellar::core::log(stellar::core::LogLevel::Warn, "UiSettings: failed to read header");
    return false;
  }
  if (header != "StellarForgeUiSettings") {
    stellar::core::log(stellar::core::LogLevel::Warn, "UiSettings: bad header");
    return false;
  }

  UiSettings s = makeDefaultUiSettings();
  s.version = version;

  std::string line;
  std::getline(f, line); // consume remainder of header line
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream ss(line);
    std::string key;
    if (!(ss >> key)) continue;
    key = lowerAscii(key);

    if (key == "autoscalefromdpi") {
      int v = 1;
      ss >> v;
      s.autoScaleFromDpi = (v != 0);
    } else if (key == "imguiinifile" || key == "inifile") {
      ss >> s.imguiIniFile;
    } else if (key == "scaleuser" || key == "scale") {
      ss >> s.scaleUser;
    } else if (key == "theme") {
      std::string v;
      ss >> v;
      if (auto t = themeFromString(v)) s.theme = *t;
    } else if (key == "autosaveonexit") {
      int v = 1;
      ss >> v;
      s.autoSaveOnExit = (v != 0);
    } else if (key == "dockenabled" || key == "dockingenabled") {
      int v = 1;
      ss >> v;
      s.dock.dockingEnabled = (v != 0);
    } else if (key == "dockpassthrucentral" || key == "dockpassthroughcentral") {
      int v = 1;
      ss >> v;
      s.dock.passthruCentral = (v != 0);
    } else if (key == "docklockcentralview") {
      int v = 1;
      ss >> v;
      s.dock.lockCentralView = (v != 0);
    } else if (key == "dockleftratio") {
      ss >> s.dock.leftRatio;
    } else if (key == "dockrightratio") {
      ss >> s.dock.rightRatio;
    } else if (key == "dockbottomratio") {
      ss >> s.dock.bottomRatio;
    }
  }

  // Clamp to sane values.
  if (s.imguiIniFile.empty()) s.imguiIniFile = "imgui.ini";
  s.scaleUser = std::clamp(s.scaleUser, 0.50f, 3.00f);
  s.dock.leftRatio = std::clamp(s.dock.leftRatio, 0.10f, 0.45f);
  s.dock.rightRatio = std::clamp(s.dock.rightRatio, 0.10f, 0.45f);
  s.dock.bottomRatio = std::clamp(s.dock.bottomRatio, 0.10f, 0.45f);

  out = s;
  return true;
}

} // namespace stellar::ui
