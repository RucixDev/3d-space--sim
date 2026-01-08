#pragma once

#include <SDL.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace stellar::game {

// Simple, app-local keybinding + control configuration.
//
// Design goals:
//  - Human-editable text file (controls.txt by default).
//  - Uses SDL scancodes so bindings are layout/locale stable.
//  - Supports Ctrl/Alt/Shift modifiers for discrete actions.
//  - Continuous axes (6DOF) use + / - scancode pairs.
//
// This intentionally lives in the game app (not the core library) so headless builds
// don't gain an SDL dependency.

struct KeyChord {
  SDL_Scancode scancode{SDL_SCANCODE_UNKNOWN};
  SDL_Keymod mods{KMOD_NONE};

  bool bound() const { return scancode != SDL_SCANCODE_UNKNOWN; }
};

struct AxisPair {
  SDL_Scancode positive{SDL_SCANCODE_UNKNOWN};
  SDL_Scancode negative{SDL_SCANCODE_UNKNOWN};
};

// Continuous (held) controls that use SDL_GetKeyboardState.
struct FlightAxes {
  AxisPair thrustForward; // + forward, - backward
  AxisPair thrustRight;   // + right, - left
  AxisPair thrustUp;      // + up, - down

  AxisPair pitch; // + pitch down, - pitch up (matches existing arrow mapping)
  AxisPair yaw;   // + yaw right, - yaw left
  AxisPair roll;  // + roll right, - roll left
};

struct FlightHolds {
  SDL_Scancode boost{SDL_SCANCODE_UNKNOWN};
  SDL_Scancode brake{SDL_SCANCODE_UNKNOWN};
  SDL_Scancode dampersEnable{SDL_SCANCODE_UNKNOWN};
  SDL_Scancode dampersDisable{SDL_SCANCODE_UNKNOWN};
};

// Discrete actions (triggered on KEYDOWN).
struct ActionBinds {
  KeyChord quit;
  KeyChord quicksave;
  KeyChord quickload;

  // Primary UI windows
  KeyChord commandPalette;
  KeyChord actionWheel;
  KeyChord toggleGalaxy;
  KeyChord toggleShip;
  KeyChord toggleMarket;
  KeyChord toggleContacts;
  KeyChord toggleMissions;
  KeyChord toggleScanner;
  KeyChord toggleTrade;
  KeyChord toggleGuide;
  KeyChord toggleHangar;
  KeyChord toggleWorldVisuals;
  KeyChord toggleSpriteLab;
  KeyChord toggleVfxLab;
  KeyChord togglePostFx;
  KeyChord toggleControlsWindow;

  // HUD/layout
  KeyChord hudLayoutToggleEdit; // Ctrl+H by default
  KeyChord hudLayoutSave;       // Ctrl+S
  KeyChord hudLayoutLoad;       // Ctrl+L
  KeyChord hudLayoutReset;      // Ctrl+R
  KeyChord toggleRadarHud;
  KeyChord toggleTacticalOverlay;

  // Flight/gameplay
  KeyChord pause;
  KeyChord toggleAutopilot;
  KeyChord navAssistApproach;
  KeyChord navAssistMatchVelocity;
  KeyChord toggleMouseSteer;
  KeyChord supercruise;
  KeyChord fsdJump;
  KeyChord scannerAction;

  KeyChord requestDockingClearance;
  KeyChord dockOrUndock;

  KeyChord cycleTargets;
  KeyChord targetStationCycle;
  KeyChord targetPlanetCycle;
  KeyChord targetContactCycle;
  KeyChord targetStar;
  KeyChord clearTarget;

  KeyChord complyOrSubmit;
  KeyChord bribe;
  KeyChord toggleCargoScoop;
  KeyChord deployCountermeasure;
};

struct ControlsConfig {
  int version{1};
  FlightAxes axes{};
  FlightHolds holds{};
  ActionBinds actions{};
};

inline std::string defaultControlsPath() { return "controls.txt"; }

inline SDL_Keymod normalizeMods(SDL_Keymod m) {
  SDL_Keymod out = KMOD_NONE;
  if (m & KMOD_CTRL) out = (SDL_Keymod)(out | KMOD_CTRL);
  if (m & KMOD_SHIFT) out = (SDL_Keymod)(out | KMOD_SHIFT);
  if (m & KMOD_ALT) out = (SDL_Keymod)(out | KMOD_ALT);
  if (m & KMOD_GUI) out = (SDL_Keymod)(out | KMOD_GUI);
  return out;
}

inline bool chordMatchesEvent(const KeyChord& chord, const SDL_KeyboardEvent& ev) {
  if (!chord.bound()) return false;
  if (ev.keysym.scancode != chord.scancode) return false;
  return normalizeMods((SDL_Keymod)ev.keysym.mod) == normalizeMods(chord.mods);
}

inline std::string lowerAscii(std::string s) {
  for (char& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

inline std::string scancodeToken(SDL_Scancode sc) {
  const char* name = SDL_GetScancodeName(sc);
  std::string s = (name && name[0]) ? std::string(name) : std::string("Unknown");
  for (char& c : s) {
    if (c == ' ') c = '_';
  }
  return s;
}

inline std::string scancodeLabel(SDL_Scancode sc) {
  const char* name = SDL_GetScancodeName(sc);
  return (name && name[0]) ? std::string(name) : std::string("(unbound)");
}

inline SDL_Scancode tokenToScancode(std::string token) {
  // Stored tokens use underscores instead of spaces.
  for (char& c : token) {
    if (c == '_') c = ' ';
  }

  // SDL's name lookup is reasonably flexible, but add a few aliases.
  const std::string k = lowerAscii(token);
  if (k == "backquote" || k == "grave" || k == "`") return SDL_SCANCODE_GRAVE;
  if (k == "lshift" || k == "leftshift") return SDL_SCANCODE_LSHIFT;
  if (k == "rshift" || k == "rightshift") return SDL_SCANCODE_RSHIFT;
  if (k == "lctrl" || k == "leftctrl" || k == "lcontrol") return SDL_SCANCODE_LCTRL;
  if (k == "rctrl" || k == "rightctrl" || k == "rcontrol") return SDL_SCANCODE_RCTRL;
  if (k == "lalt" || k == "leftalt") return SDL_SCANCODE_LALT;
  if (k == "ralt" || k == "rightalt") return SDL_SCANCODE_RALT;
  if (k == "escape" || k == "esc") return SDL_SCANCODE_ESCAPE;

  SDL_Scancode sc = SDL_GetScancodeFromName(token.c_str());
  return sc;
}

inline std::string modsToString(SDL_Keymod mods) {
  std::string s;
  mods = normalizeMods(mods);
  if (mods & KMOD_CTRL)  s += "Ctrl+";
  if (mods & KMOD_SHIFT) s += "Shift+";
  if (mods & KMOD_ALT)   s += "Alt+";
  if (mods & KMOD_GUI)   s += "Gui+";
  return s;
}

inline std::string chordToString(const KeyChord& chord) {
  if (!chord.bound()) return "(unbound)";
  return modsToString(chord.mods) + scancodeToken(chord.scancode);
}

inline std::string chordLabel(const KeyChord& chord) {
  if (!chord.bound()) return "(unbound)";
  return modsToString(chord.mods) + scancodeLabel(chord.scancode);
}

inline std::optional<KeyChord> parseChord(std::string_view chordStr) {
  std::string s(chordStr);
  if (s.empty()) return std::nullopt;
  if (lowerAscii(s) == "none" || lowerAscii(s) == "unbound") {
    return KeyChord{};
  }

  KeyChord out{};
  out.mods = KMOD_NONE;

  // Split on '+' into parts.
  std::vector<std::string> parts;
  {
    std::string cur;
    for (char c : s) {
      if (c == '+') {
        if (!cur.empty()) parts.push_back(cur);
        cur.clear();
      } else {
        cur.push_back(c);
      }
    }
    if (!cur.empty()) parts.push_back(cur);
  }
  if (parts.empty()) return std::nullopt;

  // Mods are everything except the last part.
  for (std::size_t i = 0; i + 1 < parts.size(); ++i) {
    const std::string p = lowerAscii(parts[i]);
    if (p == "ctrl" || p == "control") out.mods = (SDL_Keymod)(out.mods | KMOD_CTRL);
    else if (p == "shift") out.mods = (SDL_Keymod)(out.mods | KMOD_SHIFT);
    else if (p == "alt") out.mods = (SDL_Keymod)(out.mods | KMOD_ALT);
    else if (p == "gui" || p == "meta" || p == "cmd" || p == "win") out.mods = (SDL_Keymod)(out.mods | KMOD_GUI);
    // Unknown modifier tokens are ignored for forward compatibility.
  }

  out.scancode = tokenToScancode(parts.back());
  if (out.scancode == SDL_SCANCODE_UNKNOWN) return std::nullopt;
  return out;
}

inline double axisValue(const AxisPair& a, const Uint8* keys) {
  const double pos = (a.positive != SDL_SCANCODE_UNKNOWN && keys[a.positive]) ? 1.0 : 0.0;
  const double neg = (a.negative != SDL_SCANCODE_UNKNOWN && keys[a.negative]) ? 1.0 : 0.0;
  return pos - neg;
}

inline bool holdDown(SDL_Scancode sc, const Uint8* keys) {
  return sc != SDL_SCANCODE_UNKNOWN && keys[sc] != 0;
}

inline ControlsConfig makeDefaultControls() {
  ControlsConfig c{};
  c.version = 1;

  // 6DOF flight
  c.axes.thrustForward = {SDL_SCANCODE_W, SDL_SCANCODE_S};
  c.axes.thrustRight   = {SDL_SCANCODE_D, SDL_SCANCODE_A};
  c.axes.thrustUp      = {SDL_SCANCODE_R, SDL_SCANCODE_F};
  c.axes.pitch         = {SDL_SCANCODE_UP, SDL_SCANCODE_DOWN};
  c.axes.yaw           = {SDL_SCANCODE_RIGHT, SDL_SCANCODE_LEFT};
  c.axes.roll          = {SDL_SCANCODE_E, SDL_SCANCODE_Q};

  c.holds.boost = SDL_SCANCODE_LSHIFT;
  c.holds.brake = SDL_SCANCODE_X;
  c.holds.dampersEnable  = SDL_SCANCODE_Z;
  c.holds.dampersDisable = SDL_SCANCODE_C;

  // Discrete actions
  c.actions.quit = {SDL_SCANCODE_ESCAPE, KMOD_NONE};
  c.actions.quicksave = {SDL_SCANCODE_F5, KMOD_NONE};
  c.actions.quickload = {SDL_SCANCODE_F9, KMOD_NONE};

  // UI helpers
  c.actions.commandPalette = {SDL_SCANCODE_P, KMOD_CTRL};
  c.actions.actionWheel    = {SDL_SCANCODE_V, KMOD_ALT};

  c.actions.toggleGalaxy   = {SDL_SCANCODE_TAB, KMOD_NONE};
  c.actions.toggleShip     = {SDL_SCANCODE_F1, KMOD_NONE};
  c.actions.toggleMarket   = {SDL_SCANCODE_F2, KMOD_NONE};
  c.actions.toggleContacts = {SDL_SCANCODE_F3, KMOD_NONE};
  c.actions.toggleMissions = {SDL_SCANCODE_F4, KMOD_NONE};
  c.actions.toggleScanner  = {SDL_SCANCODE_F6, KMOD_NONE};
  c.actions.toggleTrade    = {SDL_SCANCODE_F7, KMOD_NONE};
  c.actions.toggleGuide    = {SDL_SCANCODE_F8, KMOD_NONE};
  c.actions.toggleHangar   = {SDL_SCANCODE_F9, KMOD_SHIFT};
  c.actions.toggleWorldVisuals = {SDL_SCANCODE_F10, KMOD_SHIFT};
  c.actions.toggleSpriteLab     = {SDL_SCANCODE_F10, KMOD_NONE};
  c.actions.toggleVfxLab        = {SDL_SCANCODE_F11, KMOD_NONE};
  c.actions.togglePostFx        = {SDL_SCANCODE_F12, KMOD_NONE};
  c.actions.toggleControlsWindow = {SDL_SCANCODE_K, KMOD_CTRL};

  c.actions.hudLayoutToggleEdit = {SDL_SCANCODE_H, KMOD_CTRL};
  c.actions.hudLayoutSave       = {SDL_SCANCODE_S, KMOD_CTRL};
  c.actions.hudLayoutLoad       = {SDL_SCANCODE_L, KMOD_CTRL};
  c.actions.hudLayoutReset      = {SDL_SCANCODE_R, (SDL_Keymod)(KMOD_CTRL | KMOD_SHIFT)};
  c.actions.toggleRadarHud      = {SDL_SCANCODE_R, KMOD_CTRL};
  c.actions.toggleTacticalOverlay = {SDL_SCANCODE_GRAVE, KMOD_NONE};

  c.actions.pause = {SDL_SCANCODE_SPACE, KMOD_NONE};
  c.actions.toggleAutopilot = {SDL_SCANCODE_P, KMOD_NONE};
  c.actions.navAssistApproach = {SDL_SCANCODE_PAGEUP, KMOD_NONE};
  c.actions.navAssistMatchVelocity = {SDL_SCANCODE_PAGEDOWN, KMOD_NONE};
  c.actions.toggleMouseSteer = {SDL_SCANCODE_M, KMOD_NONE};
  c.actions.supercruise = {SDL_SCANCODE_H, KMOD_NONE};
  c.actions.fsdJump = {SDL_SCANCODE_J, KMOD_NONE};
  c.actions.scannerAction = {SDL_SCANCODE_K, KMOD_NONE};

  c.actions.requestDockingClearance = {SDL_SCANCODE_L, KMOD_NONE};
  c.actions.dockOrUndock = {SDL_SCANCODE_G, KMOD_NONE};

  c.actions.cycleTargets = {SDL_SCANCODE_V, KMOD_NONE};
  c.actions.targetStationCycle = {SDL_SCANCODE_T, KMOD_NONE};
  c.actions.targetPlanetCycle  = {SDL_SCANCODE_B, KMOD_NONE};
  c.actions.targetContactCycle = {SDL_SCANCODE_N, KMOD_NONE};
  c.actions.targetStar         = {SDL_SCANCODE_U, KMOD_NONE};
  c.actions.clearTarget  = {SDL_SCANCODE_Y, KMOD_NONE};

  c.actions.complyOrSubmit = {SDL_SCANCODE_I, KMOD_NONE};
  c.actions.bribe = {SDL_SCANCODE_C, KMOD_NONE};
  c.actions.toggleCargoScoop = {SDL_SCANCODE_O, KMOD_NONE};
  c.actions.deployCountermeasure = {SDL_SCANCODE_BACKSPACE, KMOD_NONE};

  return c;
}

// ---- Serialization ----

inline void saveKeyChord(std::ostream& o, const char* name, const KeyChord& c) {
  o << "bind " << name << " " << chordToString(c) << "\n";
}

inline void saveAxis(std::ostream& o, const char* name, const AxisPair& a) {
  o << "axis " << name << " " << scancodeToken(a.positive) << " " << scancodeToken(a.negative) << "\n";
}

inline void saveHold(std::ostream& o, const char* name, SDL_Scancode sc) {
  o << "hold " << name << " " << scancodeToken(sc) << "\n";
}

inline bool saveToFile(const ControlsConfig& cfg, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) return false;
  f << "StellarForgeControls " << cfg.version << "\n";
  f << "# Chords are stored as Ctrl+Shift+Key. Key names use '_' instead of spaces.\n";
  f << "\n";

  // Axes
  f << "# Flight axes (+ -)\n";
  saveAxis(f, "ThrustForward", cfg.axes.thrustForward);
  saveAxis(f, "ThrustRight", cfg.axes.thrustRight);
  saveAxis(f, "ThrustUp", cfg.axes.thrustUp);
  saveAxis(f, "Pitch", cfg.axes.pitch);
  saveAxis(f, "Yaw", cfg.axes.yaw);
  saveAxis(f, "Roll", cfg.axes.roll);
  f << "\n";

  // Holds
  f << "# Flight holds (held keys)\n";
  saveHold(f, "Boost", cfg.holds.boost);
  saveHold(f, "Brake", cfg.holds.brake);
  saveHold(f, "DampersEnable", cfg.holds.dampersEnable);
  saveHold(f, "DampersDisable", cfg.holds.dampersDisable);
  f << "\n";

  // Actions
  f << "# Actions (trigger on key press)\n";
  const ActionBinds& a = cfg.actions;
  saveKeyChord(f, "Quit", a.quit);
  saveKeyChord(f, "QuickSave", a.quicksave);
  saveKeyChord(f, "QuickLoad", a.quickload);
  saveKeyChord(f, "CommandPalette", a.commandPalette);
  saveKeyChord(f, "ActionWheel", a.actionWheel);
  saveKeyChord(f, "ToggleGalaxy", a.toggleGalaxy);
  saveKeyChord(f, "ToggleShip", a.toggleShip);
  saveKeyChord(f, "ToggleMarket", a.toggleMarket);
  saveKeyChord(f, "ToggleContacts", a.toggleContacts);
  saveKeyChord(f, "ToggleMissions", a.toggleMissions);
  saveKeyChord(f, "ToggleScanner", a.toggleScanner);
  saveKeyChord(f, "ToggleTrade", a.toggleTrade);
  saveKeyChord(f, "ToggleGuide", a.toggleGuide);
  saveKeyChord(f, "ToggleHangar", a.toggleHangar);
  saveKeyChord(f, "ToggleWorldVisuals", a.toggleWorldVisuals);
  saveKeyChord(f, "ToggleSpriteLab", a.toggleSpriteLab);
  saveKeyChord(f, "ToggleVfxLab", a.toggleVfxLab);
  saveKeyChord(f, "TogglePostFx", a.togglePostFx);
  saveKeyChord(f, "ToggleControlsWindow", a.toggleControlsWindow);
  saveKeyChord(f, "HudEditToggle", a.hudLayoutToggleEdit);
  saveKeyChord(f, "HudSave", a.hudLayoutSave);
  saveKeyChord(f, "HudLoad", a.hudLayoutLoad);
  saveKeyChord(f, "HudReset", a.hudLayoutReset);
  saveKeyChord(f, "ToggleRadarHud", a.toggleRadarHud);
  saveKeyChord(f, "ToggleTacticalOverlay", a.toggleTacticalOverlay);
  saveKeyChord(f, "Pause", a.pause);
  saveKeyChord(f, "ToggleAutopilot", a.toggleAutopilot);
  saveKeyChord(f, "NavAssistApproach", a.navAssistApproach);
  saveKeyChord(f, "NavAssistMatchVelocity", a.navAssistMatchVelocity);
  saveKeyChord(f, "ToggleMouseSteer", a.toggleMouseSteer);
  saveKeyChord(f, "Supercruise", a.supercruise);
  saveKeyChord(f, "FsdJump", a.fsdJump);
  saveKeyChord(f, "ScannerAction", a.scannerAction);
  saveKeyChord(f, "RequestClearance", a.requestDockingClearance);
  saveKeyChord(f, "DockOrUndock", a.dockOrUndock);
  saveKeyChord(f, "CycleTargets", a.cycleTargets);
  saveKeyChord(f, "TargetStationCycle", a.targetStationCycle);
  saveKeyChord(f, "TargetPlanetCycle", a.targetPlanetCycle);
  saveKeyChord(f, "TargetContactCycle", a.targetContactCycle);
  saveKeyChord(f, "TargetStar", a.targetStar);
  saveKeyChord(f, "ClearTarget", a.clearTarget);
  saveKeyChord(f, "ComplyOrSubmit", a.complyOrSubmit);
  saveKeyChord(f, "Bribe", a.bribe);
  saveKeyChord(f, "ToggleCargoScoop", a.toggleCargoScoop);
  saveKeyChord(f, "DeployCountermeasure", a.deployCountermeasure);

  return true;
}

inline bool loadFromFile(const std::string& path, ControlsConfig& out) {
  std::ifstream f(path);
  if (!f) return false;

  std::string header;
  int version = 0;
  if (!(f >> header >> version)) return false;
  if (header != "StellarForgeControls") return false;

  ControlsConfig cfg = makeDefaultControls();
  cfg.version = version;

  // Build lookup maps for name->pointer.
  std::unordered_map<std::string, KeyChord*> bindMap;
  std::unordered_map<std::string, AxisPair*> axisMap;
  std::unordered_map<std::string, SDL_Scancode*> holdMap;
  {
    ActionBinds& a = cfg.actions;
    auto addB = [&](const char* n, KeyChord* p) { bindMap[lowerAscii(n)] = p; };
    addB("Quit", &a.quit);
    addB("QuickSave", &a.quicksave);
    addB("QuickLoad", &a.quickload);
    addB("CommandPalette", &a.commandPalette);
    addB("ActionWheel", &a.actionWheel);
    addB("ToggleGalaxy", &a.toggleGalaxy);
    addB("ToggleShip", &a.toggleShip);
    addB("ToggleMarket", &a.toggleMarket);
    addB("ToggleContacts", &a.toggleContacts);
    addB("ToggleMissions", &a.toggleMissions);
    addB("ToggleScanner", &a.toggleScanner);
    addB("ToggleTrade", &a.toggleTrade);
    addB("ToggleGuide", &a.toggleGuide);
    addB("ToggleHangar", &a.toggleHangar);
    addB("ToggleWorldVisuals", &a.toggleWorldVisuals);
    addB("ToggleSpriteLab", &a.toggleSpriteLab);
    addB("ToggleVfxLab", &a.toggleVfxLab);
    addB("TogglePostFx", &a.togglePostFx);
    addB("ToggleControlsWindow", &a.toggleControlsWindow);
    addB("HudEditToggle", &a.hudLayoutToggleEdit);
    addB("HudSave", &a.hudLayoutSave);
    addB("HudLoad", &a.hudLayoutLoad);
    addB("HudReset", &a.hudLayoutReset);
    addB("ToggleRadarHud", &a.toggleRadarHud);
    addB("ToggleTacticalOverlay", &a.toggleTacticalOverlay);
    addB("Pause", &a.pause);
    addB("ToggleAutopilot", &a.toggleAutopilot);
    addB("NavAssistApproach", &a.navAssistApproach);
    addB("NavAssistMatchVelocity", &a.navAssistMatchVelocity);
    addB("ToggleMouseSteer", &a.toggleMouseSteer);
    addB("Supercruise", &a.supercruise);
    addB("FsdJump", &a.fsdJump);
    addB("ScannerAction", &a.scannerAction);
    addB("RequestClearance", &a.requestDockingClearance);
    addB("DockOrUndock", &a.dockOrUndock);
    addB("CycleTargets", &a.cycleTargets);
    addB("TargetStationCycle", &a.targetStationCycle);
    addB("TargetPlanetCycle", &a.targetPlanetCycle);
    addB("TargetContactCycle", &a.targetContactCycle);
    addB("TargetStar", &a.targetStar);
    addB("ClearTarget", &a.clearTarget);
    addB("ComplyOrSubmit", &a.complyOrSubmit);
    addB("Bribe", &a.bribe);
    addB("ToggleCargoScoop", &a.toggleCargoScoop);
    addB("DeployCountermeasure", &a.deployCountermeasure);

    auto addA = [&](const char* n, AxisPair* p) { axisMap[lowerAscii(n)] = p; };
    addA("ThrustForward", &cfg.axes.thrustForward);
    addA("ThrustRight", &cfg.axes.thrustRight);
    addA("ThrustUp", &cfg.axes.thrustUp);
    addA("Pitch", &cfg.axes.pitch);
    addA("Yaw", &cfg.axes.yaw);
    addA("Roll", &cfg.axes.roll);

    auto addH = [&](const char* n, SDL_Scancode* p) { holdMap[lowerAscii(n)] = p; };
    addH("Boost", &cfg.holds.boost);
    addH("Brake", &cfg.holds.brake);
    addH("DampersEnable", &cfg.holds.dampersEnable);
    addH("DampersDisable", &cfg.holds.dampersDisable);
  }

  std::string line;
  std::getline(f, line); // consume rest of header line
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream iss(line);
    std::string kind;
    if (!(iss >> kind)) continue;
    kind = lowerAscii(kind);

    if (kind == "bind") {
      std::string name;
      std::string chordStr;
      if (!(iss >> name >> chordStr)) continue;
      auto it = bindMap.find(lowerAscii(name));
      if (it == bindMap.end()) continue;
      const auto parsed = parseChord(chordStr);
      if (!parsed) continue;
      *it->second = *parsed;
    } else if (kind == "axis") {
      std::string name, posTok, negTok;
      if (!(iss >> name >> posTok >> negTok)) continue;
      auto it = axisMap.find(lowerAscii(name));
      if (it == axisMap.end()) continue;
      it->second->positive = tokenToScancode(posTok);
      it->second->negative = tokenToScancode(negTok);
    } else if (kind == "hold") {
      std::string name, tok;
      if (!(iss >> name >> tok)) continue;
      auto it = holdMap.find(lowerAscii(name));
      if (it == holdMap.end()) continue;
      *it->second = tokenToScancode(tok);
    }
  }

  out = cfg;
  return true;
}

} // namespace stellar::game
