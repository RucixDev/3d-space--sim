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
//  - optionally fine-tune padding/rounding/accent colors (style overrides)
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

// Dear ImGui "multi-viewport" support (aka platform windows).
//
// When enabled, Dear ImGui can create additional OS windows for detaching
// UI panels outside the main app window (useful for multi-monitor setups).
//
// NOTE: This is separate from docking: docking is about arranging windows
// *inside* a main viewport. Viewports allow windows to leave the main OS
// window entirely.
struct ViewportSettings {
  bool enabled{false};

  // If true, platform windows won't create taskbar icons (recommended to avoid
  // clutter when you detach many windows).
  bool noTaskBarIcon{true};

  // If true, viewports won't automatically merge into the main viewport when
  // dragged over it. This can feel more "window-manager-like".
  bool noAutoMerge{false};

  // If true, platform windows are created without OS decorations.
  bool noDecoration{false};
};

// Font configuration.
//
// Dear ImGui's default embedded font is tuned for 13px with no scaling.
// When using UI scaling (HiDPI or user-controlled scaling), rebuilding the
// font atlas at the target pixel size yields sharper text than relying on
// FontGlobalScale.
struct FontSettings {
  // Optional .ttf/.otf file path. Empty => use Dear ImGui default font.
  std::string file;

  // Base font size at 1.0x UI scale.
  float sizePx{13.0f};

  // If true, build the font atlas at (sizePx * uiScale) and keep
  // io.FontGlobalScale at 1.0 for sharper text.
  // If false, build at sizePx and rely on io.FontGlobalScale.
  bool crispScaling{true};
};

// Optional fine-tuning knobs for Dear ImGui style.
//
// These apply on top of the base theme and the effective UI scale.
// They are intentionally "simple" so they can be persisted without having to
// serialize all of ImGuiStyle.
struct StyleOverrides {
  bool enabled{false};

  // Global alpha multiplier applied to the whole UI.
  // NOTE: When multi-viewport is enabled, the app may ignore this to avoid
  // semi-transparent OS windows.
  float globalAlpha{1.0f};

  // Scales padding/spacing/rounding/borders relative to the base theme.
  // 1.0 = default, <1 = compact, >1 = roomy.
  float density{1.0f};

  // Multiplier applied to rounding values (WindowRounding, FrameRounding, ...)
  // after density/scale.
  float rounding{1.0f};

  // Multiplier applied to border thickness (WindowBorderSize, FrameBorderSize, ...).
  float borderScale{1.0f};

  // Accent tint (blends selected interactive colors toward accentColor).
  bool accentEnabled{false};
  float accentColor[3]{0.26f, 0.59f, 0.98f};
  float accentStrength{0.35f};
};

struct UiSettings {
  int version{3};

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

  // Multi-monitor / detachable windows.
  ViewportSettings viewports{};

  // Font / text settings.
  FontSettings font{};

  DockSettings dock{};

  // Optional fine tuning for style.
  StyleOverrides style{};
};

// Default config file path (relative to the working directory).
std::string defaultUiSettingsPath();

UiSettings makeDefaultUiSettings();

bool saveUiSettingsToFile(const UiSettings& s, const std::string& path);
bool loadUiSettingsFromFile(const std::string& path, UiSettings& out);

} // namespace stellar::ui
