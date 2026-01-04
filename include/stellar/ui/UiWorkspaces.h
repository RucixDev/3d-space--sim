#pragma once

#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::ui {

// -----------------------------------------------------------------------------
// UI Workspaces
// -----------------------------------------------------------------------------
//
// A "workspace" is a named snapshot of:
//   - which UI windows are open/closed (by string key)
//   - which Dear ImGui ini file (layout/docking) is active
//
// The game prototype stores a LOT of UI state as loose booleans in main.cpp.
// Workspaces provide a way to quickly switch between different UI arrangements
// (e.g. Trading vs Combat vs Cinematic) without editing config files by hand.
//
// This module is intentionally headless: it has no ImGui dependency and simply
// persists a human-editable text file.

struct UiWorkspace {
  std::string name;
  std::string imguiIniFile; // e.g. "imgui.ini" or "ui_layouts/combat.ini"

  // Key -> open state.
  // Keys are application-defined but should remain stable across versions.
  // std::less<> enables heterogeneous lookup with std::string_view.
  std::map<std::string, bool, std::less<>> windows;
};

struct UiWorkspaces {
  int version{1};
  bool autoSaveOnExit{true};

  // Name of the active workspace.
  std::string active;

  std::vector<UiWorkspace> items;
};

// Default file location (relative to the working directory).
std::string defaultUiWorkspacesPath();

// Create an empty workspace set with sane defaults.
UiWorkspaces makeDefaultUiWorkspaces();

bool saveUiWorkspacesToFile(const UiWorkspaces& s, const std::string& path);
bool loadUiWorkspacesFromFile(const std::string& path, UiWorkspaces& out);

UiWorkspace* findWorkspace(UiWorkspaces& s, std::string_view name);
const UiWorkspace* findWorkspace(const UiWorkspaces& s, std::string_view name);

// Returns the active workspace, or nullptr if active is unset / not found.
const UiWorkspace* activeWorkspace(const UiWorkspaces& s);

// Insert or replace an entry by name. Returns a reference to the stored entry.
UiWorkspace& addOrUpdateWorkspace(UiWorkspaces& s, UiWorkspace ws);

// Remove a workspace by name. Returns true if removed.
bool removeWorkspace(UiWorkspaces& s, std::string_view name);

// Convenience helpers for window flags.
bool getWindowOpen(const UiWorkspace& ws, std::string_view key, bool fallback);
void setWindowOpen(UiWorkspace& ws, std::string key, bool open);

// Sanitize user-provided names (trim whitespace / remove control chars).
std::string sanitizeWorkspaceName(std::string_view name);

} // namespace stellar::ui
