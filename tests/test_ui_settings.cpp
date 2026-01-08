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

    s.style.enabled = true;
    s.style.globalAlpha = 0.86f;
    s.style.density = 0.92f;
    s.style.rounding = 1.25f;
    s.style.borderScale = 1.10f;
    s.style.accentEnabled = true;
    s.style.accentColor[0] = 0.90f;
    s.style.accentColor[1] = 0.40f;
    s.style.accentColor[2] = 0.20f;
    s.style.accentStrength = 0.42f;

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

    CHECK(out.style.enabled == s.style.enabled);
    CHECK(feq(out.style.globalAlpha, s.style.globalAlpha));
    CHECK(feq(out.style.density, s.style.density));
    CHECK(feq(out.style.rounding, s.style.rounding));
    CHECK(feq(out.style.borderScale, s.style.borderScale));
    CHECK(out.style.accentEnabled == s.style.accentEnabled);
    CHECK(feq(out.style.accentColor[0], s.style.accentColor[0]));
    CHECK(feq(out.style.accentColor[1], s.style.accentColor[1]));
    CHECK(feq(out.style.accentColor[2], s.style.accentColor[2]));
    CHECK(feq(out.style.accentStrength, s.style.accentStrength));
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
    f << "styleAlpha -1\n";        // clamp
    f << "styleDensity 99\n";      // clamp
    f << "styleRounding 99\n";     // clamp
    f << "styleBorderScale -5\n"; // clamp
    f << "styleAccentColor 2 -1 0.5\n"; // clamp
    f << "styleAccentStrength 9\n"; // clamp
    f.close();

    UiSettings out;
    CHECK(stellar::ui::loadUiSettingsFromFile(path, out));
    CHECK(out.imguiIniFile == "imgui.ini");
    CHECK(out.scaleUser <= 3.0f + 1e-4f);
    CHECK(out.dock.leftRatio >= 0.10f - 1e-4f);
    CHECK(out.dock.rightRatio <= 0.45f + 1e-4f);
    CHECK(out.dock.bottomRatio >= 0.10f - 1e-4f);

    CHECK(out.style.globalAlpha >= 0.20f - 1e-4f);
    CHECK(out.style.density <= 1.40f + 1e-4f);
    CHECK(out.style.rounding <= 3.00f + 1e-4f);
    CHECK(out.style.borderScale >= 0.00f - 1e-4f);
    CHECK(out.style.accentColor[0] <= 1.00f + 1e-4f);
    CHECK(out.style.accentColor[1] >= 0.00f - 1e-4f);
    CHECK(feq(out.style.accentColor[2], 0.50f));
    CHECK(out.style.accentStrength <= 1.00f + 1e-4f);
  }

  return failures;
}
