#include "stellar/ui/HudSettings.h"

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

static bool parseBool(const std::string& s, bool def) {
  const std::string k = lowerAscii(s);
  if (k == "1" || k == "true" || k == "yes" || k == "on" || k == "enabled") return true;
  if (k == "0" || k == "false" || k == "no" || k == "off" || k == "disabled") return false;
  return def;
}

static void writeColor(std::ostream& f, const char* key, const Color4f& c) {
  f << key << " " << c.r << " " << c.g << " " << c.b << " " << c.a << "\n";
}

static bool readColor(std::istringstream& ss, Color4f& out) {
  // Accept either RGB or RGBA.
  float r = out.r, g = out.g, b = out.b, a = out.a;
  if (!(ss >> r >> g >> b)) return false;
  if (!(ss >> a)) a = out.a;

  // Heuristic: if values exceed 1.0, assume author used 0..255 components.
  // (Supports hand-editing + copying colors from external tools.)
  const float maxv = std::max(std::max(r, g), std::max(b, a));
  if (maxv > 1.0001f) {
    r /= 255.0f;
    g /= 255.0f;
    b /= 255.0f;
    a /= 255.0f;
  }
  out = Color4f{r, g, b, a};
  return true;
}

std::string defaultHudSettingsPath() {
  return "hud_settings.txt";
}

HudSettings makeDefaultHudSettings() {
  return HudSettings{};
}

bool saveHudSettingsToFile(const HudSettings& s, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Warn, "HudSettings: failed to open for writing: " + path);
    return false;
  }

  f.setf(std::ios::fixed);
  f.precision(6);

  f << "StellarForgeHudSettings " << s.version << "\n";
  f << "autoSaveOnExit " << (s.autoSaveOnExit ? 1 : 0) << "\n";

  // Primary overlays
  f << "showRadarHud " << (s.showRadarHud ? 1 : 0) << "\n";
  f << "objectiveHudEnabled " << (s.objectiveHudEnabled ? 1 : 0) << "\n";
  f << "threatHudEnabled " << (s.threatHudEnabled ? 1 : 0) << "\n";
  f << "jumpHudEnabled " << (s.jumpHudEnabled ? 1 : 0) << "\n";

  // Radar
  f << "radarRangeKm " << s.radarRangeKm << "\n";
  f << "radarMaxBlips " << s.radarMaxBlips << "\n";

  // Misc
  f << "offscreenTargetIndicator " << (s.offscreenTargetIndicator ? 1 : 0) << "\n";

  // Combat symbology
  f << "combatHudEnabled " << (s.combatHudEnabled ? 1 : 0) << "\n";
  f << "useProceduralReticle " << (s.useProceduralReticle ? 1 : 0) << "\n";
  f << "showWeaponRings " << (s.showWeaponRings ? 1 : 0) << "\n";
  f << "showHeatRing " << (s.showHeatRing ? 1 : 0) << "\n";
  f << "showDistributorRings " << (s.showDistributorRings ? 1 : 0) << "\n";
  f << "reticleSizePx " << s.reticleSizePx << "\n";
  f << "reticleAlpha " << s.reticleAlpha << "\n";

  f << "showLeadIndicator " << (s.showLeadIndicator ? 1 : 0) << "\n";
  f << "leadUseLastFiredWeapon " << (s.leadUseLastFiredWeapon ? 1 : 0) << "\n";
  f << "leadSizePx " << s.leadSizePx << "\n";
  f << "leadMaxTimeSec " << s.leadMaxTimeSec << "\n";

  f << "showFlightPathMarker " << (s.showFlightPathMarker ? 1 : 0) << "\n";
  f << "flightMarkerUseLocalFrame " << (s.flightMarkerUseLocalFrame ? 1 : 0) << "\n";
  f << "flightMarkerClampToEdge " << (s.flightMarkerClampToEdge ? 1 : 0) << "\n";
  f << "flightMarkerSizePx " << s.flightMarkerSizePx << "\n";

  // Tactical overlay
  f << "tacticalOverlayEnabled " << (s.tacticalOverlayEnabled ? 1 : 0) << "\n";
  f << "tacticalShowLabels " << (s.tacticalShowLabels ? 1 : 0) << "\n";
  f << "tacticalRangeKm " << s.tacticalRangeKm << "\n";
  f << "tacticalMaxMarkers " << s.tacticalMaxMarkers << "\n";
  f << "tacticalShowStations " << (s.tacticalShowStations ? 1 : 0) << "\n";
  f << "tacticalShowPlanets " << (s.tacticalShowPlanets ? 1 : 0) << "\n";
  f << "tacticalShowContacts " << (s.tacticalShowContacts ? 1 : 0) << "\n";
  f << "tacticalShowCargo " << (s.tacticalShowCargo ? 1 : 0) << "\n";
  f << "tacticalShowAsteroids " << (s.tacticalShowAsteroids ? 1 : 0) << "\n";
  f << "tacticalShowSignals " << (s.tacticalShowSignals ? 1 : 0) << "\n";

  // Style / colors (v2+)
  f << "overlayBgAlpha " << s.overlayBgAlpha << "\n";
  f << "overlayBgAlphaEdit " << s.overlayBgAlphaEdit << "\n";
  f << "tintRadarIcons " << (s.tintRadarIcons ? 1 : 0) << "\n";
  f << "tintTacticalIcons " << (s.tintTacticalIcons ? 1 : 0) << "\n";

  writeColor(f, "colorPrimary", s.colorPrimary);
  writeColor(f, "colorAccent", s.colorAccent);
  writeColor(f, "colorDanger", s.colorDanger);
  writeColor(f, "colorGrid", s.colorGrid);
  writeColor(f, "colorText", s.colorText);
  writeColor(f, "colorBackground", s.colorBackground);

  return true;
}

bool loadHudSettingsFromFile(const std::string& path, HudSettings& out) {
  std::ifstream f(path);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Debug, "HudSettings: file not found: " + path);
    return false;
  }

  std::string header;
  int version = 0;
  if (!(f >> header >> version)) {
    stellar::core::log(stellar::core::LogLevel::Warn, "HudSettings: failed to read header");
    return false;
  }
  if (header != "StellarForgeHudSettings") {
    stellar::core::log(stellar::core::LogLevel::Warn, "HudSettings: bad header");
    return false;
  }

  HudSettings s = makeDefaultHudSettings();
  // Preserve the file's version but clamp it upward so saving writes a modern schema.
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

    // For boolean keys accept either integer or string tokens.
    auto readBool = [&](bool current) {
      std::string tok;
      if (!(ss >> tok)) return current;
      return parseBool(tok, current);
    };

    if (key == "autosaveonexit") {
      s.autoSaveOnExit = readBool(s.autoSaveOnExit);

    } else if (key == "showradarhud" || key == "radar" || key == "radarenabled") {
      s.showRadarHud = readBool(s.showRadarHud);
    } else if (key == "objectivehudenabled" || key == "objectivehud" || key == "objective") {
      s.objectiveHudEnabled = readBool(s.objectiveHudEnabled);
    } else if (key == "threathudenabled" || key == "threathud" || key == "threat") {
      s.threatHudEnabled = readBool(s.threatHudEnabled);
    } else if (key == "jumphudenabled" || key == "jumphud" || key == "jump") {
      s.jumpHudEnabled = readBool(s.jumpHudEnabled);

    } else if (key == "radarrangekm" || key == "radarrange") {
      ss >> s.radarRangeKm;
    } else if (key == "radarmaxblips" || key == "radarblips") {
      ss >> s.radarMaxBlips;

    } else if (key == "offscreentargetindicator" || key == "offscreentarget" || key == "targetindicator") {
      s.offscreenTargetIndicator = readBool(s.offscreenTargetIndicator);

    } else if (key == "combathudenabled" || key == "combathud") {
      s.combatHudEnabled = readBool(s.combatHudEnabled);
    } else if (key == "useproceduralreticle" || key == "proceduralreticle") {
      s.useProceduralReticle = readBool(s.useProceduralReticle);
    } else if (key == "showweaponrings" || key == "weaponrings") {
      s.showWeaponRings = readBool(s.showWeaponRings);
    } else if (key == "showheatring" || key == "heatring") {
      s.showHeatRing = readBool(s.showHeatRing);
    } else if (key == "showdistributorrings" || key == "distributorrings") {
      s.showDistributorRings = readBool(s.showDistributorRings);
    } else if (key == "reticlesizepx" || key == "reticlesize") {
      ss >> s.reticleSizePx;
    } else if (key == "reticlealpha") {
      ss >> s.reticleAlpha;

    } else if (key == "showleadindicator" || key == "leadindicator") {
      s.showLeadIndicator = readBool(s.showLeadIndicator);
    } else if (key == "leaduselastfiredweapon" || key == "leaduselastweapon") {
      s.leadUseLastFiredWeapon = readBool(s.leadUseLastFiredWeapon);
    } else if (key == "leadsizepx" || key == "leadsize") {
      ss >> s.leadSizePx;
    } else if (key == "leadmaxtimesec" || key == "leadmaxtime") {
      ss >> s.leadMaxTimeSec;

    } else if (key == "showflightpathmarker" || key == "flightpathmarker" || key == "velocitymarker") {
      s.showFlightPathMarker = readBool(s.showFlightPathMarker);
    } else if (key == "flightmarkeruselocalframe" || key == "flightmarkerlocal") {
      s.flightMarkerUseLocalFrame = readBool(s.flightMarkerUseLocalFrame);
    } else if (key == "flightmarkerclamptoedge" || key == "flightmarkerclampetoedge" || key == "flightmarkerclamp") {
      s.flightMarkerClampToEdge = readBool(s.flightMarkerClampToEdge);
    } else if (key == "flightmarkersizepx" || key == "flightmarkersize") {
      ss >> s.flightMarkerSizePx;

    } else if (key == "tacticaloverlayenabled" || key == "tacticaloverlay" || key == "tactical") {
      s.tacticalOverlayEnabled = readBool(s.tacticalOverlayEnabled);
    } else if (key == "tacticalshowlabels" || key == "tacticallabels") {
      s.tacticalShowLabels = readBool(s.tacticalShowLabels);
    } else if (key == "tacticalrangekm" || key == "tacticalrange") {
      ss >> s.tacticalRangeKm;
    } else if (key == "tacticalmaxmarkers" || key == "tacticalmarkers") {
      ss >> s.tacticalMaxMarkers;
    } else if (key == "tacticalshowstations") {
      s.tacticalShowStations = readBool(s.tacticalShowStations);
    } else if (key == "tacticalshowplanets") {
      s.tacticalShowPlanets = readBool(s.tacticalShowPlanets);
    } else if (key == "tacticalshowcontacts") {
      s.tacticalShowContacts = readBool(s.tacticalShowContacts);
    } else if (key == "tacticalshowcargo") {
      s.tacticalShowCargo = readBool(s.tacticalShowCargo);
    } else if (key == "tacticalshowasteroids") {
      s.tacticalShowAsteroids = readBool(s.tacticalShowAsteroids);
    } else if (key == "tacticalshowsignals") {
      s.tacticalShowSignals = readBool(s.tacticalShowSignals);

    } else if (key == "overlaybgalpha" || key == "overlayalpha" || key == "overlaybg") {
      ss >> s.overlayBgAlpha;
    } else if (key == "overlaybgalphaedit" || key == "overlayalphaedit" || key == "overlaybgedit") {
      ss >> s.overlayBgAlphaEdit;
    } else if (key == "tintradaricons" || key == "tintradar") {
      s.tintRadarIcons = readBool(s.tintRadarIcons);
    } else if (key == "tinttacticalicons" || key == "tinttactical") {
      s.tintTacticalIcons = readBool(s.tintTacticalIcons);

    } else if (key == "colorprimary" || key == "primary") {
      readColor(ss, s.colorPrimary);
    } else if (key == "coloraccent" || key == "accent") {
      readColor(ss, s.colorAccent);
    } else if (key == "colordanger" || key == "danger") {
      readColor(ss, s.colorDanger);
    } else if (key == "colorgrid" || key == "grid") {
      readColor(ss, s.colorGrid);
    } else if (key == "colortext" || key == "text") {
      readColor(ss, s.colorText);
    } else if (key == "colorbackground" || key == "background") {
      readColor(ss, s.colorBackground);
    }
  }

  // Clamp to sane ranges.
  s.radarRangeKm = std::clamp(s.radarRangeKm, 25000.0, 1200000.0);
  s.radarMaxBlips = std::clamp(s.radarMaxBlips, 8, 256);

  s.reticleSizePx = std::clamp(s.reticleSizePx, 8.0f, 160.0f);
  s.reticleAlpha = std::clamp(s.reticleAlpha, 0.0f, 1.0f);

  s.leadSizePx = std::clamp(s.leadSizePx, 6.0f, 120.0f);
  s.leadMaxTimeSec = std::clamp(s.leadMaxTimeSec, 1.0, 120.0);

  s.flightMarkerSizePx = std::clamp(s.flightMarkerSizePx, 6.0f, 120.0f);

  s.tacticalRangeKm = std::clamp(s.tacticalRangeKm, 20000.0, 2000000.0);
  s.tacticalMaxMarkers = std::clamp(s.tacticalMaxMarkers, 8, 512);

  s.overlayBgAlpha = std::clamp(s.overlayBgAlpha, 0.0f, 1.0f);
  s.overlayBgAlphaEdit = std::clamp(s.overlayBgAlphaEdit, 0.0f, 1.0f);

  auto clampC = [](Color4f& c) {
    c.r = std::clamp(c.r, 0.0f, 1.0f);
    c.g = std::clamp(c.g, 0.0f, 1.0f);
    c.b = std::clamp(c.b, 0.0f, 1.0f);
    c.a = std::clamp(c.a, 0.0f, 1.0f);
  };
  clampC(s.colorPrimary);
  clampC(s.colorAccent);
  clampC(s.colorDanger);
  clampC(s.colorGrid);
  clampC(s.colorText);
  clampC(s.colorBackground);

  out = s;
  return true;
}

} // namespace stellar::ui
