#include "test_harness.h"

#include "stellar/ui/HudSettings.h"

#include <cmath>
#include <fstream>

using stellar::ui::HudSettings;
using stellar::ui::Color4f;

static bool feq(double a, double b, double eps = 1e-6) {
  return std::fabs(a - b) <= eps;
}

static bool feqf(float a, float b, float eps = 1e-4f) {
  return std::fabs(a - b) <= eps;
}

static bool ceq(const Color4f& a, const Color4f& b, float eps = 1e-4f) {
  return feqf(a.r, b.r, eps) && feqf(a.g, b.g, eps) && feqf(a.b, b.b, eps) && feqf(a.a, b.a, eps);
}

int test_hud_settings() {
  int failures = 0;

  const std::string path = "test_hud_settings_tmp.txt";

  // Roundtrip save/load.
  {
    HudSettings s = stellar::ui::makeDefaultHudSettings();
    s.autoSaveOnExit = false;
    s.showRadarHud = false;
    s.objectiveHudEnabled = false;
    s.threatHudEnabled = true;
    s.jumpHudEnabled = true;
    s.radarRangeKm = 555555.0;
    s.radarMaxBlips = 123;
    s.offscreenTargetIndicator = false;

    s.combatHudEnabled = true;
    s.useProceduralReticle = false;
    s.showWeaponRings = false;
    s.showHeatRing = true;
    s.showDistributorRings = false;
    s.reticleSizePx = 77.0f;
    s.reticleAlpha = 0.25f;

    s.showLeadIndicator = false;
    s.leadUseLastFiredWeapon = false;
    s.leadSizePx = 44.0f;
    s.leadMaxTimeSec = 9.0;

    s.showFlightPathMarker = true;
    s.flightMarkerUseLocalFrame = false;
    s.flightMarkerClampToEdge = false;
    s.flightMarkerSizePx = 33.0f;

    s.tacticalOverlayEnabled = false;
    s.tacticalShowLabels = false;
    s.tacticalRangeKm = 999999.0;
    s.tacticalMaxMarkers = 222;
    s.tacticalShowStations = false;
    s.tacticalShowPlanets = false;
    s.tacticalShowContacts = false;
    s.tacticalShowCargo = false;
    s.tacticalShowAsteroids = false;
    s.tacticalShowSignals = false;

    // Style / colors (v2)
    s.overlayBgAlpha = 0.12f;
    s.overlayBgAlphaEdit = 0.97f;
    s.tintRadarIcons = true;
    s.tintTacticalIcons = true;
    s.colorPrimary = Color4f{0.10f, 0.20f, 0.30f, 0.40f};
    s.colorAccent = Color4f{0.20f, 0.30f, 0.40f, 0.50f};
    s.colorDanger = Color4f{0.90f, 0.10f, 0.10f, 0.80f};
    s.colorGrid = Color4f{0.01f, 0.02f, 0.03f, 1.00f};
    s.colorText = Color4f{0.80f, 0.81f, 0.82f, 0.83f};
    s.colorBackground = Color4f{0.00f, 0.00f, 0.00f, 0.25f};

    CHECK(stellar::ui::saveHudSettingsToFile(s, path));

    HudSettings out;
    CHECK(stellar::ui::loadHudSettingsFromFile(path, out));

    CHECK(out.autoSaveOnExit == s.autoSaveOnExit);
    CHECK(out.showRadarHud == s.showRadarHud);
    CHECK(out.objectiveHudEnabled == s.objectiveHudEnabled);
    CHECK(out.threatHudEnabled == s.threatHudEnabled);
    CHECK(out.jumpHudEnabled == s.jumpHudEnabled);
    CHECK(feq(out.radarRangeKm, s.radarRangeKm));
    CHECK(out.radarMaxBlips == s.radarMaxBlips);
    CHECK(out.offscreenTargetIndicator == s.offscreenTargetIndicator);

    CHECK(out.combatHudEnabled == s.combatHudEnabled);
    CHECK(out.useProceduralReticle == s.useProceduralReticle);
    CHECK(out.showWeaponRings == s.showWeaponRings);
    CHECK(out.showHeatRing == s.showHeatRing);
    CHECK(out.showDistributorRings == s.showDistributorRings);
    CHECK(std::fabs(out.reticleSizePx - s.reticleSizePx) <= 1e-4f);
    CHECK(std::fabs(out.reticleAlpha - s.reticleAlpha) <= 1e-4f);

    CHECK(out.showLeadIndicator == s.showLeadIndicator);
    CHECK(out.leadUseLastFiredWeapon == s.leadUseLastFiredWeapon);
    CHECK(std::fabs(out.leadSizePx - s.leadSizePx) <= 1e-4f);
    CHECK(feq(out.leadMaxTimeSec, s.leadMaxTimeSec));

    CHECK(out.showFlightPathMarker == s.showFlightPathMarker);
    CHECK(out.flightMarkerUseLocalFrame == s.flightMarkerUseLocalFrame);
    CHECK(out.flightMarkerClampToEdge == s.flightMarkerClampToEdge);
    CHECK(std::fabs(out.flightMarkerSizePx - s.flightMarkerSizePx) <= 1e-4f);

    CHECK(out.tacticalOverlayEnabled == s.tacticalOverlayEnabled);
    CHECK(out.tacticalShowLabels == s.tacticalShowLabels);
    CHECK(feq(out.tacticalRangeKm, s.tacticalRangeKm));
    CHECK(out.tacticalMaxMarkers == s.tacticalMaxMarkers);
    CHECK(out.tacticalShowStations == s.tacticalShowStations);
    CHECK(out.tacticalShowPlanets == s.tacticalShowPlanets);
    CHECK(out.tacticalShowContacts == s.tacticalShowContacts);
    CHECK(out.tacticalShowCargo == s.tacticalShowCargo);
    CHECK(out.tacticalShowAsteroids == s.tacticalShowAsteroids);
    CHECK(out.tacticalShowSignals == s.tacticalShowSignals);

    CHECK(feqf(out.overlayBgAlpha, s.overlayBgAlpha));
    CHECK(feqf(out.overlayBgAlphaEdit, s.overlayBgAlphaEdit));
    CHECK(out.tintRadarIcons == s.tintRadarIcons);
    CHECK(out.tintTacticalIcons == s.tintTacticalIcons);
    CHECK(ceq(out.colorPrimary, s.colorPrimary));
    CHECK(ceq(out.colorAccent, s.colorAccent));
    CHECK(ceq(out.colorDanger, s.colorDanger));
    CHECK(ceq(out.colorGrid, s.colorGrid));
    CHECK(ceq(out.colorText, s.colorText));
    CHECK(ceq(out.colorBackground, s.colorBackground));
  }

  // Clamp & resilience: values outside ranges should be clamped.
  {
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    f << "StellarForgeHudSettings 1\n";
    f << "radarRangeKm 1\n";           // clamp up
    f << "radarMaxBlips 9999\n";        // clamp down
    f << "reticleSizePx -5\n";          // clamp
    f << "reticleAlpha 123\n";          // clamp
    f << "leadMaxTimeSec -50\n";        // clamp
    f << "tacticalRangeKm 99999999\n";  // clamp
    f << "tacticalMaxMarkers -2\n";     // clamp
    f << "overlayBgAlpha -10\n";         // clamp
    f << "overlayBgAlphaEdit 10\n";      // clamp
    f << "colorPrimary 2 -3 0.5 9\n";    // clamp
    f << "unknownKey 42\n";             // ignored
    f.close();

    HudSettings out;
    CHECK(stellar::ui::loadHudSettingsFromFile(path, out));

    CHECK(out.radarRangeKm >= 25000.0 - 1e-4);
    CHECK(out.radarMaxBlips <= 256);
    CHECK(out.reticleSizePx >= 8.0f - 1e-4f);
    CHECK(out.reticleAlpha <= 1.0f + 1e-4f);
    CHECK(out.leadMaxTimeSec >= 1.0 - 1e-4);
    CHECK(out.tacticalRangeKm <= 2000000.0 + 1e-4);
    CHECK(out.tacticalMaxMarkers >= 8);

    CHECK(out.overlayBgAlpha >= 0.0f - 1e-4f);
    CHECK(out.overlayBgAlphaEdit <= 1.0f + 1e-4f);
    CHECK(out.colorPrimary.r <= 1.0f + 1e-4f);
    CHECK(out.colorPrimary.g >= 0.0f - 1e-4f);
    CHECK(out.colorPrimary.a <= 1.0f + 1e-4f);
  }

  return failures;
}
