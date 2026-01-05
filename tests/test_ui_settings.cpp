#include "test_harness.h"

#include "stellar/ui/UiSettings.h"

#include <cmath>
#include <fstream>

using stellar::ui::UiSettings;

static bool feq(float a, float b, float eps = 1e-5f) {
  return std::fabs(a - b) <= eps;
}

int test_ui_settings() {
  int failures = 0;

  const std::string path = "test_ui_settings_tmp.txt";

  // Roundtrip save/load.
  {
    UiSettings s = stellar::ui::makeDefaultUiSettings();
    s.imguiIniFile = "ui_layouts/combat.ini";
    s.autoScaleFromDpi = false;
    s.scaleUser = 1.234f;
    s.theme = stellar::ui::UiTheme::HighContrast;
    s.autoSaveOnExit = false;
    s.viewports.enabled = true;
    s.viewports.noTaskBarIcon = false;
    s.viewports.noAutoMerge = true;
    s.viewports.noDecoration = true;
    s.dock.dockingEnabled = false;
    s.dock.passthruCentral = false;
    s.dock.lockCentralView = false;
    s.dock.leftRatio = 0.33f;
    s.dock.rightRatio = 0.22f;
    s.dock.bottomRatio = 0.18f;

    CHECK(stellar::ui::saveUiSettingsToFile(s, path));

    UiSettings out;
    CHECK(stellar::ui::loadUiSettingsFromFile(path, out));
    CHECK(out.imguiIniFile == s.imguiIniFile);
    CHECK(out.autoScaleFromDpi == s.autoScaleFromDpi);
    CHECK(feq(out.scaleUser, s.scaleUser));
    CHECK(out.theme == s.theme);
    CHECK(out.autoSaveOnExit == s.autoSaveOnExit);
    CHECK(out.viewports.enabled == s.viewports.enabled);
    CHECK(out.viewports.noTaskBarIcon == s.viewports.noTaskBarIcon);
    CHECK(out.viewports.noAutoMerge == s.viewports.noAutoMerge);
    CHECK(out.viewports.noDecoration == s.viewports.noDecoration);
    CHECK(out.dock.dockingEnabled == s.dock.dockingEnabled);
    CHECK(out.dock.passthruCentral == s.dock.passthruCentral);
    CHECK(out.dock.lockCentralView == s.dock.lockCentralView);
    CHECK(feq(out.dock.leftRatio, s.dock.leftRatio));
    CHECK(feq(out.dock.rightRatio, s.dock.rightRatio));
    CHECK(feq(out.dock.bottomRatio, s.dock.bottomRatio));
  }

  // Clamp & resilience: values outside ranges should be clamped to sane limits.
  {
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    f << "StellarForgeUiSettings 1\n";
    f << "imguiIniFile \n";        // empty -> default
    f << "scaleUser 99\n";         // clamp
    f << "dockLeftRatio 0.01\n";   // clamp
    f << "dockRightRatio 0.99\n";  // clamp
    f << "dockBottomRatio -5\n";   // clamp
    f.close();

    UiSettings out;
    CHECK(stellar::ui::loadUiSettingsFromFile(path, out));
    CHECK(out.imguiIniFile == "imgui.ini");
    CHECK(out.scaleUser <= 3.0f + 1e-4f);
    CHECK(out.dock.leftRatio >= 0.10f - 1e-4f);
    CHECK(out.dock.rightRatio <= 0.45f + 1e-4f);
    CHECK(out.dock.bottomRatio >= 0.10f - 1e-4f);
  }

  return failures;
}
