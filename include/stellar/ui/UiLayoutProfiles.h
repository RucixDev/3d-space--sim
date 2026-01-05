#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace stellar::ui {

// -----------------------------------------------------------------------------
// UI Layout Profiles (Dear ImGui ini files)
// -----------------------------------------------------------------------------
//
// Dear ImGui stores window positions, sizes and docking layout in an ini file.
// This helper provides a small, headless API to manage multiple ini profiles
// stored under a dedicated folder (default: "ui_layouts/").
//
// IMPORTANT: This module does NOT depend on Dear ImGui. It only deals with file
// names/paths and discovery. The game app is responsible for calling
// ImGui::LoadIniSettingsFromDisk / SaveIniSettingsToDisk.

struct UiLayoutProfile {
  std::string name;     // e.g. "combat"
  std::string path;     // e.g. "ui_layouts/combat.ini" or "imgui.ini"
  bool isDefault{false};
};

// Default directory (relative to working directory).
std::string defaultUiLayoutsDir();

// Sanitize an arbitrary string into a filename stem suitable for "<stem>.ini".
//
// Rules:
//  - keeps ASCII [a-z0-9], '_' and '-'
//  - converts whitespace to '_'
//  - lowercases A-Z
//  - collapses repeated '_' and trims leading/trailing '_'
//  - returns "profile" if the result is empty
std::string sanitizeLayoutStem(std::string_view name);

// Returns "<dir>/<sanitize(name)>.ini". Does not create directories.
std::string layoutIniPath(std::string_view name, std::string_view dir = {});

// Returns the available profiles.
//
// The first entry is always the default (name "Default", path defaultIniPath).
// The remaining entries are .ini files directly under `dir`, sorted by name.
std::vector<UiLayoutProfile> listLayoutProfiles(const std::string& dir = defaultUiLayoutsDir(),
                                                const std::string& defaultIniPath = "imgui.ini");

// Best-effort: create the layouts dir if missing.
bool ensureLayoutDir(const std::string& dir = defaultUiLayoutsDir());

// Returns true iff `path` refers to a direct child ".ini" under `dir`.
// This is used as a safety check before delete/rename operations.
bool isSafeLayoutIniPath(std::string_view dir, std::string_view path);

// Returns a unique path "<dir>/<stem>.ini", adding suffixes "_2", "_3", ...
// if needed.
std::string uniqueLayoutIniPath(std::string_view desiredStem,
                                const std::string& dir = defaultUiLayoutsDir());

// Basic file ops. These do not enforce safety checks; call isSafeLayoutIniPath
// before operating on user-provided paths.
bool copyLayoutIniFile(const std::string& srcPath, const std::string& dstPath, bool overwrite);
bool renameLayoutIniFile(const std::string& srcPath, const std::string& dstPath, bool overwrite);
bool deleteLayoutIniFile(const std::string& path);

} // namespace stellar::ui
