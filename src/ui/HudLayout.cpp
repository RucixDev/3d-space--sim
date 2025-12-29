#include "stellar/ui/HudLayout.h"

#include "stellar/core/Log.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

namespace stellar::ui {

static std::string lowerAscii(std::string s) {
  for (char& c : s) {
    c = (char)std::tolower((unsigned char)c);
  }
  return s;
}

std::string defaultHudLayoutPath() {
  return "hud_layout.txt";
}

HudLayout makeDefaultHudLayout() {
  HudLayout l{};
  l.version = 1;
  l.safeMarginPx = 14.0f;

  // Defaults are expressed in normalized coordinates (0..1) to stay stable across
  // resolutions.
  //
  // Use a ~2% inset from the edges. Pivot points match the window anchor.
  {
    auto& w = l.widget(HudWidgetId::Radar);
    w.posNormX = 0.985f;
    w.posNormY = 0.985f;
    w.pivotX = 1.0f;
    w.pivotY = 1.0f;
    w.scale = 1.0f;
    w.enabled = true;
  }
  {
    auto& w = l.widget(HudWidgetId::Objective);
    w.posNormX = 0.985f;
    w.posNormY = 0.02f;
    w.pivotX = 1.0f;
    w.pivotY = 0.0f;
    w.scale = 1.0f;
    w.enabled = true;
  }
  {
    auto& w = l.widget(HudWidgetId::Threat);
    w.posNormX = 0.02f;
    w.posNormY = 0.02f;
    w.pivotX = 0.0f;
    w.pivotY = 0.0f;
    w.scale = 1.0f;
    w.enabled = true;
  }
  {
    auto& w = l.widget(HudWidgetId::Jump);
    w.posNormX = 0.5f;
    w.posNormY = 0.02f;
    w.pivotX = 0.5f;
    w.pivotY = 0.0f;
    w.scale = 1.0f;
    w.enabled = true;
  }

  return l;
}

const char* toString(HudWidgetId id) {
  switch (id) {
    case HudWidgetId::Radar: return "Radar";
    case HudWidgetId::Objective: return "Objective";
    case HudWidgetId::Threat: return "Threat";
    case HudWidgetId::Jump: return "Jump";
    default: return "Unknown";
  }
}

std::optional<HudWidgetId> widgetFromString(const std::string& s) {
  const std::string k = lowerAscii(s);
  if (k == "radar") return HudWidgetId::Radar;
  if (k == "objective" || k == "objectives") return HudWidgetId::Objective;
  if (k == "threat" || k == "pirate" || k == "pirates") return HudWidgetId::Threat;
  if (k == "jump" || k == "fsd" || k == "hyperspace") return HudWidgetId::Jump;
  return std::nullopt;
}

bool saveToFile(const HudLayout& layout, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Warn, "HUDLayout: failed to open for writing: " + path);
    return false;
  }

  f.setf(std::ios::fixed);
  f.precision(6);

  f << "StellarForgeHudLayout " << layout.version << "\n";
  f << "global safeMarginPx " << layout.safeMarginPx << "\n";

  for (std::size_t i = 0; i < layout.widgets.size(); ++i) {
    const auto id = (HudWidgetId)i;
    const auto& w = layout.widgets[i];
    f << "widget " << toString(id)
      << " pos " << w.posNormX << " " << w.posNormY
      << " pivot " << w.pivotX << " " << w.pivotY
      << " scale " << w.scale
      << " enabled " << (w.enabled ? 1 : 0)
      << "\n";
  }
  return true;
}

bool loadFromFile(const std::string& path, HudLayout& out) {
  std::ifstream f(path);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Debug, "HUDLayout: file not found: " + path);
    return false;
  }

  std::string header;
  int version = 0;
  if (!(f >> header >> version)) {
    stellar::core::log(stellar::core::LogLevel::Warn, "HUDLayout: failed to read header");
    return false;
  }
  if (header != "StellarForgeHudLayout") {
    stellar::core::log(stellar::core::LogLevel::Warn, "HUDLayout: bad header");
    return false;
  }

  HudLayout loaded = makeDefaultHudLayout();
  loaded.version = version;

  // Parse line-by-line so future versions can add keys without breaking older builds.
  std::string line;
  std::getline(f, line); // consume remainder of header line
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream ss(line);
    std::string tok;
    if (!(ss >> tok)) continue;

    tok = lowerAscii(tok);
    if (tok == "global") {
      std::string key;
      while (ss >> key) {
        key = lowerAscii(key);
        if (key == "safemarginpx") {
          ss >> loaded.safeMarginPx;
        } else {
          // Unknown global key: best-effort skip one value.
          std::string skip;
          ss >> skip;
        }
      }
      continue;
    }

    if (tok == "widget") {
      std::string name;
      if (!(ss >> name)) continue;
      const auto idOpt = widgetFromString(name);
      if (!idOpt) continue;
      auto& w = loaded.widget(*idOpt);

      std::string key;
      while (ss >> key) {
        key = lowerAscii(key);
        if (key == "pos") {
          ss >> w.posNormX >> w.posNormY;
        } else if (key == "pivot") {
          ss >> w.pivotX >> w.pivotY;
        } else if (key == "scale") {
          ss >> w.scale;
        } else if (key == "enabled") {
          int v = 1;
          ss >> v;
          w.enabled = (v != 0);
        } else {
          // Unknown key: best-effort skip one value.
          std::string skip;
          ss >> skip;
        }
      }

      // Clamp to sane values.
      w.posNormX = std::clamp(w.posNormX, 0.0f, 1.0f);
      w.posNormY = std::clamp(w.posNormY, 0.0f, 1.0f);
      w.pivotX = std::clamp(w.pivotX, 0.0f, 1.0f);
      w.pivotY = std::clamp(w.pivotY, 0.0f, 1.0f);
      w.scale = std::clamp(w.scale, 0.50f, 2.50f);
    }
  }

  out = loaded;
  return true;
}

} // namespace stellar::ui
