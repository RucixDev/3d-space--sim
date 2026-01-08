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

  // Fonts
  // Empty fontFile => use Dear ImGui embedded default.
  f << "fontFile " << (s.font.file.empty() ? std::string("(default)") : s.font.file) << "\n";
  f << "fontSizePx " << s.font.sizePx << "\n";
  f << "fontCrispScaling " << (s.font.crispScaling ? 1 : 0) << "\n";

  // Multi-viewport (platform windows)
  f << "viewportsEnabled " << (s.viewports.enabled ? 1 : 0) << "\n";
  f << "viewportsNoTaskBarIcon " << (s.viewports.noTaskBarIcon ? 1 : 0) << "\n";
  f << "viewportsNoAutoMerge " << (s.viewports.noAutoMerge ? 1 : 0) << "\n";
  f << "viewportsNoDecoration " << (s.viewports.noDecoration ? 1 : 0) << "\n";

  f << "dockEnabled " << (s.dock.dockingEnabled ? 1 : 0) << "\n";
  f << "dockPassthruCentral " << (s.dock.passthruCentral ? 1 : 0) << "\n";
  f << "dockLockCentralView " << (s.dock.lockCentralView ? 1 : 0) << "\n";
  f << "dockLeftRatio " << s.dock.leftRatio << "\n";
  f << "dockRightRatio " << s.dock.rightRatio << "\n";
  f << "dockBottomRatio " << s.dock.bottomRatio << "\n";

  // Style overrides (optional fine-tuning)
  f << "styleOverridesEnabled " << (s.style.enabled ? 1 : 0) << "\n";
  f << "styleAlpha " << s.style.globalAlpha << "\n";
  f << "styleDensity " << s.style.density << "\n";
  f << "styleRounding " << s.style.rounding << "\n";
  f << "styleBorderScale " << s.style.borderScale << "\n";
  f << "styleAccentEnabled " << (s.style.accentEnabled ? 1 : 0) << "\n";
  f << "styleAccentColor " << s.style.accentColor[0] << " " << s.style.accentColor[1] << " " << s.style.accentColor[2] << "\n";
  f << "styleAccentStrength " << s.style.accentStrength << "\n";

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
  // Auto-upgrade version so saving writes the latest schema.
  s.version = std::max(version, s.version);

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
    } else if (key == "fontfile" || key == "font" || key == "fontpath") {
      std::string v;
      ss >> v;
      if (v == "(default)" || v == "default") v.clear();
      s.font.file = v;
    } else if (key == "fontsizepx" || key == "fontsize" || key == "font_size_px") {
      ss >> s.font.sizePx;
    } else if (key == "fontcrispscaling" || key == "fontcrispscale" || key == "fontbuildatscale") {
      int v = 1;
      ss >> v;
      s.font.crispScaling = (v != 0);
    } else if (key == "viewportsenabled" || key == "multiviewport" || key == "multiviewports" || key == "viewports") {
      int v = 0;
      ss >> v;
      s.viewports.enabled = (v != 0);
    } else if (key == "viewportsnotaskbaricon" || key == "viewportnotaskbaricon") {
      int v = 1;
      ss >> v;
      s.viewports.noTaskBarIcon = (v != 0);
    } else if (key == "viewportsnoautomerge") {
      int v = 0;
      ss >> v;
      s.viewports.noAutoMerge = (v != 0);
    } else if (key == "viewportsnodecoration") {
      int v = 0;
      ss >> v;
      s.viewports.noDecoration = (v != 0);
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
    } else if (key == "styleoverridesenabled" || key == "styleenabled" || key == "style_overrides_enabled") {
      int v = 0;
      ss >> v;
      s.style.enabled = (v != 0);
    } else if (key == "stylealpha" || key == "style_global_alpha" || key == "styleopacity") {
      ss >> s.style.globalAlpha;
    } else if (key == "styledensity" || key == "style_density") {
      ss >> s.style.density;
    } else if (key == "stylerounding" || key == "style_rounding") {
      ss >> s.style.rounding;
    } else if (key == "styleborderscale" || key == "style_border_scale" || key == "styleborder") {
      ss >> s.style.borderScale;
    } else if (key == "styleaccentenabled" || key == "style_accent_enabled") {
      int v = 0;
      ss >> v;
      s.style.accentEnabled = (v != 0);
    } else if (key == "styleaccentcolor" || key == "style_accent_color") {
      float r = s.style.accentColor[0];
      float g = s.style.accentColor[1];
      float b = s.style.accentColor[2];
      ss >> r >> g >> b;
      s.style.accentColor[0] = r;
      s.style.accentColor[1] = g;
      s.style.accentColor[2] = b;
    } else if (key == "styleaccentstrength" || key == "style_accent_strength") {
      ss >> s.style.accentStrength;
    }
  }

  // Clamp to sane values.
  if (s.imguiIniFile.empty()) s.imguiIniFile = "imgui.ini";
  s.scaleUser = std::clamp(s.scaleUser, 0.50f, 3.00f);
  s.font.sizePx = std::clamp(s.font.sizePx, 10.0f, 32.0f);
  s.dock.leftRatio = std::clamp(s.dock.leftRatio, 0.10f, 0.45f);
  s.dock.rightRatio = std::clamp(s.dock.rightRatio, 0.10f, 0.45f);
  s.dock.bottomRatio = std::clamp(s.dock.bottomRatio, 0.10f, 0.45f);

  // Style overrides
  s.style.globalAlpha = std::clamp(s.style.globalAlpha, 0.20f, 1.00f);
  s.style.density = std::clamp(s.style.density, 0.70f, 1.40f);
  s.style.rounding = std::clamp(s.style.rounding, 0.00f, 3.00f);
  s.style.borderScale = std::clamp(s.style.borderScale, 0.00f, 3.00f);
  s.style.accentStrength = std::clamp(s.style.accentStrength, 0.00f, 1.00f);
  for (float& c : s.style.accentColor) c = std::clamp(c, 0.0f, 1.0f);

  out = s;
  return true;
}

} // namespace stellar::ui
