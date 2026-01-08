#include "ControlsWindow.h"

#include "stellar/ui/FuzzySearch.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace stellar::game {

static bool anyNonSpace(const char* s) {
  if (!s) return false;
  for (const unsigned char* p = (const unsigned char*)s; *p; ++p) {
    if (std::isspace(*p) == 0) return true;
  }
  return false;
}

static void drawHighlightedText(std::string_view text, const std::vector<int>& positions, ImU32 colText, ImU32 colHi) {
  if (positions.empty()) {
    ImGui::TextUnformatted(text.data(), text.data() + text.size());
    return;
  }

  // Build a mask for contiguous range coalescing.
  std::vector<unsigned char> mark(text.size(), 0);
  for (int p : positions) {
    if (p >= 0 && (std::size_t)p < text.size()) mark[(std::size_t)p] = 1;
  }

  ImDrawList* dl = ImGui::GetWindowDrawList();
  ImVec2 cur = ImGui::GetCursorScreenPos();
  const float lineH = ImGui::GetTextLineHeight();

  const char* base = text.data();

  std::size_t i = 0;
  while (i < text.size()) {
    const bool hi = mark[i] != 0;
    std::size_t j = i + 1;
    while (j < text.size() && (mark[j] != 0) == hi) ++j;

    const char* segB = base + i;
    const char* segE = base + j;
    dl->AddText(cur, hi ? colHi : colText, segB, segE);

    const ImVec2 sz = ImGui::CalcTextSize(segB, segE);
    cur.x += sz.x;
    i = j;
  }

  // Advance cursor.
  ImGui::Dummy(ImVec2(cur.x - ImGui::GetCursorScreenPos().x, lineH));
}

static void drawHighlightedText(std::string_view text, const std::vector<int>& positions) {
  drawHighlightedText(text, positions,
                      ImGui::GetColorU32(ImGuiCol_Text),
                      ImGui::GetColorU32(ImGuiCol_PlotHistogram));
}

// A stable key for action chords.
struct ChordKey {
  SDL_Scancode sc{SDL_SCANCODE_UNKNOWN};
  SDL_Keymod mods{KMOD_NONE};
  bool operator==(const ChordKey& o) const { return sc == o.sc && mods == o.mods; }
};

struct ChordKeyHash {
  std::size_t operator()(const ChordKey& k) const noexcept {
    // SDL_Scancode is an enum; SDL_Keymod is bitmask.
    return (std::size_t)k.sc * 1315423911u ^ (std::size_t)k.mods;
  }
};

void openControlsWindow(ControlsWindowState& st, bool focusFilter) {
  st.open = true;
  st.focusFilter = focusFilter;
}

bool handleControlsRebindKeydown(ControlsWindowState& st,
                                const SDL_KeyboardEvent& ev,
                                ControlsConfig& controls,
                                bool& controlsDirty,
                                const ToastFn& toastFn) {
  (void)controls;

  if (st.rebind.kind == ControlsWindowState::RebindKind::None) return false;

  const SDL_Scancode sc = ev.keysym.scancode;
  const SDL_Keymod mods = normalizeMods((SDL_Keymod)ev.keysym.mod);

  // ESC cancels the capture (bind Escape via file-edit if you really want it).
  if (sc == SDL_SCANCODE_ESCAPE && mods == KMOD_NONE) {
    if (toastFn) toastFn("Keybind capture canceled.", 1.2);
    st.rebind = ControlsWindowState::RebindRequest{};
    return true;
  }

  // Backspace/Delete clears the binding.
  if (sc == SDL_SCANCODE_BACKSPACE || sc == SDL_SCANCODE_DELETE) {
    if (st.rebind.kind == ControlsWindowState::RebindKind::Action && st.rebind.chord) {
      *st.rebind.chord = KeyChord{};
    } else if (st.rebind.kind == ControlsWindowState::RebindKind::Hold && st.rebind.hold) {
      *st.rebind.hold = SDL_SCANCODE_UNKNOWN;
    } else if (st.rebind.kind == ControlsWindowState::RebindKind::AxisPositive && st.rebind.axis) {
      st.rebind.axis->positive = SDL_SCANCODE_UNKNOWN;
    } else if (st.rebind.kind == ControlsWindowState::RebindKind::AxisNegative && st.rebind.axis) {
      st.rebind.axis->negative = SDL_SCANCODE_UNKNOWN;
    }

    controlsDirty = true;
    if (toastFn) toastFn(std::string("Unbound: ") + st.rebind.label, 1.6);
    st.rebind = ControlsWindowState::RebindRequest{};
    return true;
  }

  // Apply new binding.
  if (st.rebind.kind == ControlsWindowState::RebindKind::Action && st.rebind.chord) {
    *st.rebind.chord = KeyChord{sc, mods};
    if (toastFn) toastFn("Bound " + st.rebind.label + " -> " + chordLabel(*st.rebind.chord), 1.8);
  } else if (st.rebind.kind == ControlsWindowState::RebindKind::Hold && st.rebind.hold) {
    *st.rebind.hold = sc;
    if (toastFn) toastFn("Bound " + st.rebind.label + " -> " + scancodeLabel(sc), 1.8);
  } else if (st.rebind.kind == ControlsWindowState::RebindKind::AxisPositive && st.rebind.axis) {
    st.rebind.axis->positive = sc;
    if (toastFn) toastFn("Bound " + st.rebind.label + " (+) -> " + scancodeLabel(sc), 1.8);
  } else if (st.rebind.kind == ControlsWindowState::RebindKind::AxisNegative && st.rebind.axis) {
    st.rebind.axis->negative = sc;
    if (toastFn) toastFn("Bound " + st.rebind.label + " (-) -> " + scancodeLabel(sc), 1.8);
  }

  controlsDirty = true;
  st.rebind = ControlsWindowState::RebindRequest{};
  return true;
}

static std::string sanitizeProfileName(const char* in) {
  std::string out;
  for (const unsigned char* p = (const unsigned char*)in; *p; ++p) {
    const char c = (char)*p;
    if (std::isalnum(*p) || c == '_' || c == '-') out.push_back(c);
  }
  if (out.empty()) out = "profile";
  return out;
}

// Returns (score, labelPositions). If query is empty, returns score 0.
static std::pair<int, std::vector<int>> matchLabel(std::string_view query,
                                                   std::string_view label,
                                                   std::string_view detail) {
  if (!anyNonSpace(query.data())) return {0, {}};

  std::string hay;
  hay.reserve(label.size() + detail.size() + 2);
  hay += label;
  if (!detail.empty()) {
    hay.push_back(' ');
    hay += detail;
  }

  const auto r = stellar::ui::fuzzyMatch(query, hay);
  if (r.score < 0) return {-1, {}};

  std::vector<int> labelPos;
  const int split = (int)label.size();
  for (int p : r.positions) {
    if (p >= 0 && p < split) labelPos.push_back(p);
  }
  return {r.score, std::move(labelPos)};
}

void drawControlsWindow(ControlsWindowState& st,
                        ControlsConfig& controls,
                        const ControlsConfig& defaults,
                        bool& controlsDirty,
                        bool& autoSaveOnExit,
                        std::string& controlsPath,
                        const ToastFn& toastFn) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(760.0f, 610.0f), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Controls", &st.open)) {
    ImGui::End();
    return;
  }

  // Ctrl+F focuses the filter when the window is focused.
  if (st.rebind.kind == ControlsWindowState::RebindKind::None) {
    const bool focused = ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows);
    if (focused && ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_F)) {
      st.focusFilter = true;
    }
  }

  ImGui::TextDisabled(
      "Rebind controls & keybinds. Click Rebind, then press a key. Esc cancels capture; Backspace/Delete clears.");
  ImGui::TextDisabled("Tip: Ctrl+F focuses the filter. Use profiles to keep multiple key layouts.");

  if (st.rebind.kind != ControlsWindowState::RebindKind::None) {
    ImGui::Separator();
    ImGui::Text("Capturing: %s", st.rebind.label.c_str());
    ImGui::SameLine();
    if (ImGui::Button("Cancel capture")) {
      st.rebind = ControlsWindowState::RebindRequest{};
    }
  }

  ImGui::Separator();

  // ---- Profiles ----
  if (ImGui::CollapsingHeader("Profiles", ImGuiTreeNodeFlags_DefaultOpen)) {
    namespace fs = std::filesystem;
    const fs::path dir("controls_profiles");
    struct Profile { std::string name; std::string path; bool isDefault; };

    std::vector<Profile> profiles;
    profiles.push_back(Profile{"Default", defaultControlsPath(), true});
    if (fs::exists(dir) && fs::is_directory(dir)) {
      for (const auto& ent : fs::directory_iterator(dir)) {
        if (!ent.is_regular_file()) continue;
        const fs::path p = ent.path();
        if (p.extension() != ".txt") continue;
        profiles.push_back(Profile{p.stem().string(), p.string(), false});
      }
    }

    std::sort(profiles.begin() + 1, profiles.end(), [](const Profile& a, const Profile& b) {
      return a.name < b.name;
    });

    // Keep selection stable.
    int curIdx = 0;
    for (int i = 0; i < (int)profiles.size(); ++i) {
      if (profiles[(std::size_t)i].path == controlsPath) { curIdx = i; break; }
    }
    if (st.selectedProfile < 0 || st.selectedProfile >= (int)profiles.size()) st.selectedProfile = curIdx;

    ImGui::TextDisabled("Active file: %s", controlsPath.c_str());

    if (ImGui::BeginListBox("##controls_profiles", ImVec2(-FLT_MIN, 120.0f))) {
      for (int i = 0; i < (int)profiles.size(); ++i) {
        const bool sel = (st.selectedProfile == i);
        std::string label = profiles[(std::size_t)i].name;
        if (profiles[(std::size_t)i].path == controlsPath) label += "  (active)";
        if (ImGui::Selectable(label.c_str(), sel)) st.selectedProfile = i;
      }
      ImGui::EndListBox();
    }

    const bool selOk = (st.selectedProfile >= 0 && st.selectedProfile < (int)profiles.size());
    const bool selIsDefault = selOk ? profiles[(std::size_t)st.selectedProfile].isDefault : true;
    if (!selOk) ImGui::BeginDisabled();
    if (ImGui::Button("Activate selected") && selOk) {
      ControlsConfig loaded = makeDefaultControls();
      const std::string path = profiles[(std::size_t)st.selectedProfile].path;
      if (loadFromFile(path, loaded)) {
        controls = loaded;
        controlsDirty = false;
        controlsPath = path;
        if (toastFn) toastFn("Activated controls profile: " + profiles[(std::size_t)st.selectedProfile].name, 1.8);
      } else {
        if (toastFn) toastFn("Failed to load controls profile.", 2.0);
      }
    }
    ImGui::SameLine();
    if (ImGui::Button("Save current")) {
      const bool ok = saveToFile(controls, controlsPath);
      if (ok) controlsDirty = false;
      if (toastFn) toastFn(ok ? "Saved controls." : "Failed to save controls.", 1.8);
    }
    if (!selOk) ImGui::EndDisabled();

    ImGui::SameLine();
    ImGui::Checkbox("Auto-save on exit", &autoSaveOnExit);
    if (controlsDirty) {
      ImGui::SameLine();
      ImGui::TextColored(ImVec4(1.0f, 0.85f, 0.25f, 1.0f), "*unsaved*");
    }

    // Save-as workflow.
    ImGui::Separator();
    ImGui::InputTextWithHint("##new_controls_profile", "New profile name (letters/numbers/-/_)", st.newProfileName, sizeof(st.newProfileName));
    ImGui::SameLine();
    if (ImGui::Button("Save As")) {
      const std::string safe = sanitizeProfileName(st.newProfileName);
      fs::create_directories(dir);
      fs::path p = dir / (safe + ".txt");
      const bool ok = saveToFile(controls, p.string());
      if (ok) {
        controlsPath = p.string();
        controlsDirty = false;
        if (toastFn) toastFn("Saved new controls profile.", 1.8);
      } else {
        if (toastFn) toastFn("Failed to save new profile.", 1.8);
      }
    }
    ImGui::SameLine();
    if (selIsDefault) ImGui::BeginDisabled();
    if (ImGui::Button("Delete selected") && selOk && !selIsDefault) {
      const std::string path = profiles[(std::size_t)st.selectedProfile].path;
      // Safety: only delete files inside controls_profiles.
      if (path.rfind("controls_profiles", 0) == 0) {
        std::error_code ec;
        fs::remove(fs::path(path), ec);
        if (!ec) {
          if (toastFn) toastFn("Deleted controls profile.", 1.6);
          if (controlsPath == path) {
            controlsPath = defaultControlsPath();
          }
        } else {
          if (toastFn) toastFn("Failed to delete profile.", 1.6);
        }
      } else {
        if (toastFn) toastFn("Refusing to delete outside controls_profiles/", 2.0);
      }
    }
    if (selIsDefault) ImGui::EndDisabled();
  }

  ImGui::Separator();

  // ---- Filter / diagnostics ----
  if (st.focusFilter) {
    ImGui::SetKeyboardFocusHere();
    st.focusFilter = false;
  }
  ImGui::PushItemWidth(-1);
  if (ImGui::InputTextWithHint("##controls_filter", "Filter controls... (fuzzy)", st.filter, sizeof(st.filter), ImGuiInputTextFlags_AutoSelectAll)) {
    // keep
  }
  ImGui::PopItemWidth();

  ImGui::Checkbox("Changed only", &st.showChangedOnly);
  ImGui::SameLine();
  ImGui::Checkbox("Conflicts only", &st.showConflictsOnly);
  ImGui::SameLine();
  ImGui::Checkbox("Sort by relevance", &st.sortByRelevance);

  // Build conflict maps.
  std::unordered_map<ChordKey, std::vector<std::string>, ChordKeyHash> actionConf;
  std::unordered_map<SDL_Scancode, std::vector<std::string>> holdConf;
  std::unordered_map<SDL_Scancode, std::vector<std::string>> axisConf;

  auto addAction = [&](const char* name, const KeyChord& chord) {
    if (!chord.bound()) return;
    ChordKey k{chord.scancode, normalizeMods(chord.mods)};
    actionConf[k].push_back(name);
  };
  const auto& a = controls.actions;
  addAction("Quit", a.quit);
  addAction("Quick save", a.quicksave);
  addAction("Quick load", a.quickload);
  addAction("Command Palette", a.commandPalette);
  addAction("Action Wheel", a.actionWheel);
  addAction("Toggle Galaxy window", a.toggleGalaxy);
  addAction("Toggle Ship window", a.toggleShip);
  addAction("Toggle Market window", a.toggleMarket);
  addAction("Toggle Contacts window", a.toggleContacts);
  addAction("Toggle Missions window", a.toggleMissions);
  addAction("Toggle Scanner window", a.toggleScanner);
  addAction("Toggle Trade window", a.toggleTrade);
  addAction("Toggle Guide window", a.toggleGuide);
  addAction("Toggle Hangar window", a.toggleHangar);
  addAction("Toggle World Visuals window", a.toggleWorldVisuals);
  addAction("Toggle Sprite Lab", a.toggleSpriteLab);
  addAction("Toggle VFX Lab", a.toggleVfxLab);
  addAction("Toggle Post FX", a.togglePostFx);
  addAction("Toggle Controls window", a.toggleControlsWindow);

  addAction("HUD Layout: edit mode", a.hudLayoutToggleEdit);
  addAction("HUD Layout: save", a.hudLayoutSave);
  addAction("HUD Layout: load", a.hudLayoutLoad);
  addAction("HUD Layout: reset", a.hudLayoutReset);
  addAction("Toggle Radar HUD", a.toggleRadarHud);
  addAction("Toggle Tactical overlay", a.toggleTacticalOverlay);

  addAction("Pause", a.pause);
  addAction("Toggle Autopilot", a.toggleAutopilot);
  addAction("Nav Assist: Approach", a.navAssistApproach);
  addAction("Nav Assist: Match Velocity", a.navAssistMatchVelocity);
  addAction("Toggle Mouse steer", a.toggleMouseSteer);
  addAction("Supercruise", a.supercruise);
  addAction("FSD Jump", a.fsdJump);
  addAction("Scanner action", a.scannerAction);
  addAction("Docking clearance", a.requestDockingClearance);
  addAction("Dock/Undock", a.dockOrUndock);
  addAction("Cycle targets", a.cycleTargets);
  addAction("Target: station", a.targetStationCycle);
  addAction("Target: planet", a.targetPlanetCycle);
  addAction("Target: contact", a.targetContactCycle);
  addAction("Target: star", a.targetStar);
  addAction("Clear target", a.clearTarget);
  addAction("Comply / submit", a.complyOrSubmit);
  addAction("Bribe", a.bribe);
  addAction("Toggle Cargo Scoop", a.toggleCargoScoop);
  addAction("Deploy Countermeasure", a.deployCountermeasure);

  auto addHold = [&](const char* name, SDL_Scancode sc) {
    if (sc == SDL_SCANCODE_UNKNOWN) return;
    holdConf[sc].push_back(name);
  };
  addHold("Boost", controls.holds.boost);
  addHold("Brake", controls.holds.brake);
  addHold("Dampers enable", controls.holds.dampersEnable);
  addHold("Dampers disable", controls.holds.dampersDisable);

  auto addAxisUse = [&](const char* name, SDL_Scancode sc) {
    if (sc == SDL_SCANCODE_UNKNOWN) return;
    axisConf[sc].push_back(name);
  };
  addAxisUse("Thrust forward +", controls.axes.thrustForward.positive);
  addAxisUse("Thrust forward -", controls.axes.thrustForward.negative);
  addAxisUse("Thrust right +", controls.axes.thrustRight.positive);
  addAxisUse("Thrust right -", controls.axes.thrustRight.negative);
  addAxisUse("Thrust up +", controls.axes.thrustUp.positive);
  addAxisUse("Thrust up -", controls.axes.thrustUp.negative);
  addAxisUse("Pitch +", controls.axes.pitch.positive);
  addAxisUse("Pitch -", controls.axes.pitch.negative);
  addAxisUse("Yaw +", controls.axes.yaw.positive);
  addAxisUse("Yaw -", controls.axes.yaw.negative);
  addAxisUse("Roll +", controls.axes.roll.positive);
  addAxisUse("Roll -", controls.axes.roll.negative);

  auto countConflicts = [&]() {
    int c = 0;
    for (const auto& kv : actionConf) if (kv.second.size() > 1) ++c;
    for (const auto& kv : holdConf) if (kv.second.size() > 1) ++c;
    for (const auto& kv : axisConf) if (kv.second.size() > 1) ++c;
    return c;
  };

  const int conflictGroups = countConflicts();
  if (conflictGroups > 0) {
    ImGui::SameLine();
    ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.45f, 1.0f), "Conflicts: %d", conflictGroups);
  } else {
    ImGui::SameLine();
    ImGui::TextDisabled("Conflicts: 0");
  }

  ImGui::Separator();

  const std::string query = st.filter;
  const bool hasQuery = anyNonSpace(query.c_str());

  // ---- Draw tables ----
  auto requestActionRebind = [&](const char* label, KeyChord& chord) {
    st.rebind.kind = ControlsWindowState::RebindKind::Action;
    st.rebind.chord = &chord;
    st.rebind.hold = nullptr;
    st.rebind.axis = nullptr;
    st.rebind.label = label;
  };

  auto requestHoldRebind = [&](const char* label, SDL_Scancode& sc) {
    st.rebind.kind = ControlsWindowState::RebindKind::Hold;
    st.rebind.chord = nullptr;
    st.rebind.hold = &sc;
    st.rebind.axis = nullptr;
    st.rebind.label = label;
  };

  auto requestAxisRebind = [&](const char* label, AxisPair& axis, bool positive) {
    st.rebind.kind = positive ? ControlsWindowState::RebindKind::AxisPositive : ControlsWindowState::RebindKind::AxisNegative;
    st.rebind.chord = nullptr;
    st.rebind.hold = nullptr;
    st.rebind.axis = &axis;
    st.rebind.label = label;
  };

  auto rowShouldShow = [&](bool changed, bool conflict, int score) {
    if (st.showChangedOnly && !changed) return false;
    if (st.showConflictsOnly && !conflict) return false;
    if (hasQuery && score < 0) return false;
    return true;
  };

  // ------- Axes -------
  if (ImGui::CollapsingHeader("Flight Axes", ImGuiTreeNodeFlags_DefaultOpen)) {
    struct AxisRow {
      const char* label;
      AxisPair* axis;
      const AxisPair* def;
      int score;
      std::vector<int> pos;
      bool changed;
      bool conflict;
    };
    std::vector<AxisRow> rows;
    rows.reserve(8);

    auto add = [&](const char* label, AxisPair& axis, const AxisPair& def) {
      const std::string curKeys = std::string(scancodeLabel(axis.positive)) + " / " + scancodeLabel(axis.negative);
      const auto [score, pos] = matchLabel(query, label, curKeys);
      const bool changed = axis.positive != def.positive || axis.negative != def.negative;
      bool conflict = false;
      if (axis.positive != SDL_SCANCODE_UNKNOWN) {
        const auto it = axisConf.find(axis.positive);
        if (it != axisConf.end() && it->second.size() > 1) conflict = true;
      }
      if (axis.negative != SDL_SCANCODE_UNKNOWN) {
        const auto it = axisConf.find(axis.negative);
        if (it != axisConf.end() && it->second.size() > 1) conflict = true;
      }
      if (axis.positive != SDL_SCANCODE_UNKNOWN && axis.positive == axis.negative) conflict = true;
      rows.push_back({label, &axis, &def, score, pos, changed, conflict});
    };

    add("Thrust forward", controls.axes.thrustForward, defaults.axes.thrustForward);
    add("Thrust right", controls.axes.thrustRight, defaults.axes.thrustRight);
    add("Thrust up", controls.axes.thrustUp, defaults.axes.thrustUp);
    add("Pitch", controls.axes.pitch, defaults.axes.pitch);
    add("Yaw", controls.axes.yaw, defaults.axes.yaw);
    add("Roll", controls.axes.roll, defaults.axes.roll);

    // Filter
    std::vector<AxisRow> filtered;
    filtered.reserve(rows.size());
    for (auto& r : rows) {
      if (!rowShouldShow(r.changed, r.conflict, r.score)) continue;
      filtered.push_back(std::move(r));
    }
    if (st.sortByRelevance && hasQuery) {
      std::stable_sort(filtered.begin(), filtered.end(), [](const AxisRow& a, const AxisRow& b) {
        return a.score > b.score;
      });
    }

    if (ImGui::BeginTable("##controls_axes", 6,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_ScrollY,
                          ImVec2(0.0f, 180.0f))) {
      ImGui::TableSetupColumn("Axis");
      ImGui::TableSetupColumn("Keys");
      ImGui::TableSetupColumn("Default");
      ImGui::TableSetupColumn("Rebind");
      ImGui::TableSetupColumn("Reset");
      ImGui::TableSetupColumn("Clear");
      ImGui::TableHeadersRow();

      for (auto& r : filtered) {
        AxisPair& axis = *r.axis;
        const AxisPair& def = *r.def;

        ImGui::TableNextRow();

        ImGui::TableSetColumnIndex(0);
        if (!st.scrollToLabel.empty() && st.scrollToLabel == r.label) {
          ImGui::SetScrollHereY(0.35f);
          st.scrollToLabel.clear();
        }
        if (r.conflict) {
          ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.45f, 1.0f), "!");
          ImGui::SameLine();
        }
        if (hasQuery) {
          drawHighlightedText(r.label, r.pos);
        } else {
          ImGui::TextUnformatted(r.label);
        }

        ImGui::TableSetColumnIndex(1);
        ImGui::Text("%s / %s", scancodeLabel(axis.positive).c_str(), scancodeLabel(axis.negative).c_str());

        ImGui::TableSetColumnIndex(2);
        ImGui::TextDisabled("%s / %s", scancodeLabel(def.positive).c_str(), scancodeLabel(def.negative).c_str());

        ImGui::TableSetColumnIndex(3);
        if (ImGui::SmallButton((std::string("Rebind +##") + r.label).c_str())) {
          requestAxisRebind(r.label, axis, true);
        }
        ImGui::SameLine();
        if (ImGui::SmallButton((std::string("Rebind -##") + r.label).c_str())) {
          requestAxisRebind(r.label, axis, false);
        }

        ImGui::TableSetColumnIndex(4);
        if (!r.changed) ImGui::BeginDisabled();
        if (ImGui::SmallButton((std::string("Reset##") + r.label).c_str())) {
          axis = def;
          controlsDirty = true;
          if (toastFn) toastFn(std::string("Reset: ") + r.label, 1.2);
        }
        if (!r.changed) ImGui::EndDisabled();

        ImGui::TableSetColumnIndex(5);
        if (ImGui::SmallButton((std::string("Clear##") + r.label).c_str())) {
          axis.positive = SDL_SCANCODE_UNKNOWN;
          axis.negative = SDL_SCANCODE_UNKNOWN;
          controlsDirty = true;
          if (toastFn) toastFn(std::string("Unbound: ") + r.label, 1.6);
        }
      }

      ImGui::EndTable();
    }
  }

  // ------- Holds -------
  if (ImGui::CollapsingHeader("Flight Holds", ImGuiTreeNodeFlags_DefaultOpen)) {
    struct HoldRow {
      const char* label;
      SDL_Scancode* sc;
      SDL_Scancode def;
      int score;
      std::vector<int> pos;
      bool changed;
      bool conflict;
    };
    std::vector<HoldRow> rows;
    rows.reserve(8);

    auto add = [&](const char* label, SDL_Scancode& sc, SDL_Scancode def) {
      const std::string cur = scancodeLabel(sc);
      const auto [score, pos] = matchLabel(query, label, cur);
      const bool changed = sc != def;
      bool conflict = false;
      if (sc != SDL_SCANCODE_UNKNOWN) {
        const auto it = holdConf.find(sc);
        if (it != holdConf.end() && it->second.size() > 1) conflict = true;
      }
      rows.push_back({label, &sc, def, score, pos, changed, conflict});
    };

    add("Boost", controls.holds.boost, defaults.holds.boost);
    add("Brake", controls.holds.brake, defaults.holds.brake);
    add("Dampers enable", controls.holds.dampersEnable, defaults.holds.dampersEnable);
    add("Dampers disable", controls.holds.dampersDisable, defaults.holds.dampersDisable);

    std::vector<HoldRow> filtered;
    filtered.reserve(rows.size());
    for (auto& r : rows) {
      if (!rowShouldShow(r.changed, r.conflict, r.score)) continue;
      filtered.push_back(std::move(r));
    }
    if (st.sortByRelevance && hasQuery) {
      std::stable_sort(filtered.begin(), filtered.end(), [](const HoldRow& a, const HoldRow& b) {
        return a.score > b.score;
      });
    }

    if (ImGui::BeginTable("##controls_holds", 6,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_ScrollY,
                          ImVec2(0.0f, 140.0f))) {
      ImGui::TableSetupColumn("Hold");
      ImGui::TableSetupColumn("Key");
      ImGui::TableSetupColumn("Default");
      ImGui::TableSetupColumn("Rebind");
      ImGui::TableSetupColumn("Reset");
      ImGui::TableSetupColumn("Clear");
      ImGui::TableHeadersRow();

      for (auto& r : filtered) {
        ImGui::TableNextRow();

        ImGui::TableSetColumnIndex(0);
        if (!st.scrollToLabel.empty() && st.scrollToLabel == r.label) {
          ImGui::SetScrollHereY(0.35f);
          st.scrollToLabel.clear();
        }
        if (r.conflict) {
          ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.45f, 1.0f), "!");
          ImGui::SameLine();
        }
        if (hasQuery) {
          drawHighlightedText(r.label, r.pos);
        } else {
          ImGui::TextUnformatted(r.label);
        }

        ImGui::TableSetColumnIndex(1);
        ImGui::TextUnformatted(scancodeLabel(*r.sc).c_str());

        ImGui::TableSetColumnIndex(2);
        ImGui::TextDisabled("%s", scancodeLabel(r.def).c_str());

        ImGui::TableSetColumnIndex(3);
        if (ImGui::SmallButton((std::string("Rebind##") + r.label).c_str())) {
          requestHoldRebind(r.label, *r.sc);
        }

        ImGui::TableSetColumnIndex(4);
        if (!r.changed) ImGui::BeginDisabled();
        if (ImGui::SmallButton((std::string("Reset##") + r.label).c_str())) {
          *r.sc = r.def;
          controlsDirty = true;
          if (toastFn) toastFn(std::string("Reset: ") + r.label, 1.2);
        }
        if (!r.changed) ImGui::EndDisabled();

        ImGui::TableSetColumnIndex(5);
        if (ImGui::SmallButton((std::string("Clear##") + r.label).c_str())) {
          *r.sc = SDL_SCANCODE_UNKNOWN;
          controlsDirty = true;
          if (toastFn) toastFn(std::string("Unbound: ") + r.label, 1.6);
        }
      }

      ImGui::EndTable();
    }
  }

  // ------- Actions -------
  if (ImGui::CollapsingHeader("Actions", ImGuiTreeNodeFlags_DefaultOpen)) {
    struct ActionRow {
      const char* label;
      KeyChord* chord;
      const KeyChord* def;
      int score;
      std::vector<int> pos;
      bool changed;
      bool conflict;
    };

    std::vector<ActionRow> rows;
    rows.reserve(128);

    auto add = [&](const char* label, KeyChord& chord, const KeyChord& def) {
      const std::string cur = chordLabel(chord);
      const auto [score, pos] = matchLabel(query, label, cur);
      const bool changed = chord.scancode != def.scancode || normalizeMods(chord.mods) != normalizeMods(def.mods);
      bool conflict = false;
      if (chord.bound()) {
        ChordKey k{chord.scancode, normalizeMods(chord.mods)};
        const auto it = actionConf.find(k);
        if (it != actionConf.end() && it->second.size() > 1) conflict = true;
      }
      rows.push_back({label, &chord, &def, score, pos, changed, conflict});
    };

    const auto& ad = defaults.actions;
    add("Quit", controls.actions.quit, ad.quit);
    add("Quick save", controls.actions.quicksave, ad.quicksave);
    add("Quick load", controls.actions.quickload, ad.quickload);
    add("Command Palette", controls.actions.commandPalette, ad.commandPalette);

    add("Action Wheel", controls.actions.actionWheel, ad.actionWheel);

    add("Toggle Galaxy window", controls.actions.toggleGalaxy, ad.toggleGalaxy);
    add("Toggle Ship window", controls.actions.toggleShip, ad.toggleShip);
    add("Toggle Market window", controls.actions.toggleMarket, ad.toggleMarket);
    add("Toggle Contacts window", controls.actions.toggleContacts, ad.toggleContacts);
    add("Toggle Missions window", controls.actions.toggleMissions, ad.toggleMissions);
    add("Toggle Scanner window", controls.actions.toggleScanner, ad.toggleScanner);
    add("Toggle Trade window", controls.actions.toggleTrade, ad.toggleTrade);
    add("Toggle Guide window", controls.actions.toggleGuide, ad.toggleGuide);
    add("Toggle Hangar window", controls.actions.toggleHangar, ad.toggleHangar);
    add("Toggle World Visuals window", controls.actions.toggleWorldVisuals, ad.toggleWorldVisuals);
    add("Toggle Sprite Lab", controls.actions.toggleSpriteLab, ad.toggleSpriteLab);
    add("Toggle VFX Lab", controls.actions.toggleVfxLab, ad.toggleVfxLab);
    add("Toggle Post FX", controls.actions.togglePostFx, ad.togglePostFx);
    add("Toggle Controls window", controls.actions.toggleControlsWindow, ad.toggleControlsWindow);

    add("Toggle Radar HUD", controls.actions.toggleRadarHud, ad.toggleRadarHud);
    add("Toggle Tactical overlay", controls.actions.toggleTacticalOverlay, ad.toggleTacticalOverlay);
    add("HUD Layout: edit mode", controls.actions.hudLayoutToggleEdit, ad.hudLayoutToggleEdit);
    add("HUD Layout: save", controls.actions.hudLayoutSave, ad.hudLayoutSave);
    add("HUD Layout: load", controls.actions.hudLayoutLoad, ad.hudLayoutLoad);
    add("HUD Layout: reset", controls.actions.hudLayoutReset, ad.hudLayoutReset);

    add("Pause", controls.actions.pause, ad.pause);
    add("Toggle Autopilot", controls.actions.toggleAutopilot, ad.toggleAutopilot);
    add("Nav Assist: Approach", controls.actions.navAssistApproach, ad.navAssistApproach);
    add("Nav Assist: Match Velocity", controls.actions.navAssistMatchVelocity, ad.navAssistMatchVelocity);
    add("Mouse steer", controls.actions.toggleMouseSteer, ad.toggleMouseSteer);
    add("Supercruise", controls.actions.supercruise, ad.supercruise);
    add("FSD Jump", controls.actions.fsdJump, ad.fsdJump);
    add("Scanner action", controls.actions.scannerAction, ad.scannerAction);

    add("Docking clearance", controls.actions.requestDockingClearance, ad.requestDockingClearance);
    add("Dock/Undock", controls.actions.dockOrUndock, ad.dockOrUndock);

    add("Cycle targets", controls.actions.cycleTargets, ad.cycleTargets);
    add("Target: station", controls.actions.targetStationCycle, ad.targetStationCycle);
    add("Target: planet", controls.actions.targetPlanetCycle, ad.targetPlanetCycle);
    add("Target: contact", controls.actions.targetContactCycle, ad.targetContactCycle);
    add("Target: star", controls.actions.targetStar, ad.targetStar);
    add("Clear target", controls.actions.clearTarget, ad.clearTarget);

    add("Comply / submit", controls.actions.complyOrSubmit, ad.complyOrSubmit);
    add("Bribe", controls.actions.bribe, ad.bribe);
    add("Toggle Cargo Scoop", controls.actions.toggleCargoScoop, ad.toggleCargoScoop);
    add("Deploy Countermeasure", controls.actions.deployCountermeasure, ad.deployCountermeasure);

    std::vector<ActionRow> filtered;
    filtered.reserve(rows.size());
    for (auto& r : rows) {
      if (!rowShouldShow(r.changed, r.conflict, r.score)) continue;
      filtered.push_back(std::move(r));
    }
    if (st.sortByRelevance && hasQuery) {
      std::stable_sort(filtered.begin(), filtered.end(), [](const ActionRow& a, const ActionRow& b) {
        return a.score > b.score;
      });
    }

    if (ImGui::BeginTable("##controls_actions", 7,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_ScrollY,
                          ImVec2(0.0f, 280.0f))) {
      ImGui::TableSetupColumn("Action");
      ImGui::TableSetupColumn("Key");
      ImGui::TableSetupColumn("Default");
      ImGui::TableSetupColumn("Rebind");
      ImGui::TableSetupColumn("Reset");
      ImGui::TableSetupColumn("Clear");
      ImGui::TableSetupColumn("Copy");
      ImGui::TableHeadersRow();

      for (auto& r : filtered) {
        KeyChord& chord = *r.chord;
        const KeyChord& def = *r.def;

        ImGui::TableNextRow();

        ImGui::TableSetColumnIndex(0);
        if (!st.scrollToLabel.empty() && st.scrollToLabel == r.label) {
          ImGui::SetScrollHereY(0.35f);
          st.scrollToLabel.clear();
        }
        if (r.conflict) {
          ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.45f, 1.0f), "!");
          ImGui::SameLine();
        }
        if (hasQuery) {
          drawHighlightedText(r.label, r.pos);
        } else {
          ImGui::TextUnformatted(r.label);
        }

        ImGui::TableSetColumnIndex(1);
        if (r.conflict) ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.55f, 0.45f, 1.0f));
        ImGui::TextUnformatted(chordLabel(chord).c_str());
        if (r.conflict) {
          ImGui::PopStyleColor();
          if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
            ImGui::BeginTooltip();
            ImGui::TextUnformatted("Conflict: multiple actions share this binding:");
            if (chord.bound()) {
              ChordKey k{chord.scancode, normalizeMods(chord.mods)};
              const auto it = actionConf.find(k);
              if (it != actionConf.end()) {
                for (const auto& name : it->second) ImGui::BulletText("%s", name.c_str());
              }
            }
            ImGui::EndTooltip();
          }
        }

        ImGui::TableSetColumnIndex(2);
        ImGui::TextDisabled("%s", chordLabel(def).c_str());

        ImGui::TableSetColumnIndex(3);
        if (ImGui::SmallButton((std::string("Rebind##") + r.label).c_str())) {
          requestActionRebind(r.label, chord);
        }

        ImGui::TableSetColumnIndex(4);
        if (!r.changed) ImGui::BeginDisabled();
        if (ImGui::SmallButton((std::string("Reset##") + r.label).c_str())) {
          chord = def;
          controlsDirty = true;
          if (toastFn) toastFn(std::string("Reset: ") + r.label, 1.2);
        }
        if (!r.changed) ImGui::EndDisabled();

        ImGui::TableSetColumnIndex(5);
        if (ImGui::SmallButton((std::string("Clear##") + r.label).c_str())) {
          chord = KeyChord{};
          controlsDirty = true;
          if (toastFn) toastFn(std::string("Unbound: ") + r.label, 1.6);
        }

        ImGui::TableSetColumnIndex(6);
        if (ImGui::SmallButton((std::string("Copy##") + r.label).c_str())) {
          const std::string s = chordToString(chord);
          ImGui::SetClipboardText(s.c_str());
          if (toastFn) toastFn("Copied binding to clipboard.", 1.0);
        }
      }

      ImGui::EndTable();
    }
  }

  ImGui::End();
}

} // namespace stellar::game
