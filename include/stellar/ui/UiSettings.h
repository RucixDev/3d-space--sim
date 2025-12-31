#pragma once

#include <optional>
#include <string>

namespace stellar::ui {

// User-facing UI settings that persist across runs.
//
// This keeps high-level preferences (theme, scaling, docking knobs) separate from
// Dear ImGui's own imgui.ini state (window positions, sizes, docking layout).
//
// The goal is to allow players to:
//  - tune UI readability (scale)
//  - pick an accessible theme (e.g. high contrast)
//  - tweak the default dockspace ratios
// without having to delete their entire imgui.ini.

enum class UiTheme : int {
  Dark = 0,
  Light,
  Classic,
  HighContrast,
};

const char* toString(UiTheme t);
std::optional<UiTheme> themeFromString(const std::string& s);

struct DockSettings {
  bool dockingEnabled{true};
  bool passthruCentral{true};
  bool lockCentralView{true};
  float leftRatio{0.26f};
  float rightRatio{0.30f};
  float bottomRatio{0.25f};
};

struct UiSettings {
  int version{1};

  // Dear ImGui .ini file used for window positions/docking layout.
  // By default this is "imgui.ini" in the working directory.
  std::string imguiIniFile{"imgui.ini"};

  // If true, the app derives a DPI scale factor from SDL and multiplies it with
  // `scaleUser`. This is useful on HiDPI displays.
  bool autoScaleFromDpi{true};

  // User multiplier applied on top of the optional DPI factor.
  float scaleUser{1.0f};

  UiTheme theme{UiTheme::Dark};

  // Save this file on exit even if the user didn't click "Save".
  bool autoSaveOnExit{true};

  DockSettings dock{};
};

// Default config file path (relative to the working directory).
std::string defaultUiSettingsPath();

UiSettings makeDefaultUiSettings();

bool saveUiSettingsToFile(const UiSettings& s, const std::string& path);
bool loadUiSettingsFromFile(const std::string& path, UiSettings& out);

} // namespace stellar::ui
