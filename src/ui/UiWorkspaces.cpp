#include "stellar/ui/UiWorkspaces.h"

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

static std::string trimAscii(std::string s) {
  std::size_t a = 0;
  while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
  std::size_t b = s.size();
  while (b > a && std::isspace((unsigned char)s[b - 1])) --b;
  return s.substr(a, b - a);
}

static std::string sanitizeKey(std::string s) {
  // Window keys are expected to be simple tokens, but we still clean up.
  for (char& c : s) {
    if (c == '\n' || c == '\r' || c == '\t') c = ' ';
  }
  s = trimAscii(std::move(s));
  // Keys are used for lookups, so keep them compact.
  // Replace spaces with underscores to avoid ambiguous parsing.
  for (char& c : s) {
    if (c == ' ') c = '_';
  }
  return s;
}

std::string sanitizeWorkspaceName(std::string_view name) {
  std::string s(name);
  for (char& c : s) {
    if (c == '\n' || c == '\r' || c == '\t') c = ' ';
  }
  s = trimAscii(std::move(s));
  if (s.empty()) s = "Workspace";
  return s;
}

std::string defaultUiWorkspacesPath() {
  return "ui_workspaces.txt";
}

UiWorkspaces makeDefaultUiWorkspaces() {
  return UiWorkspaces{};
}

bool getWindowOpen(const UiWorkspace& ws, std::string_view key, bool fallback) {
  const auto it = ws.windows.find(key);
  if (it == ws.windows.end()) return fallback;
  return it->second;
}

void setWindowOpen(UiWorkspace& ws, std::string key, bool open) {
  key = sanitizeKey(std::move(key));
  if (key.empty()) return;
  ws.windows[key] = open;
}

UiWorkspace* findWorkspace(UiWorkspaces& s, std::string_view name) {
  for (auto& ws : s.items) {
    if (ws.name == name) return &ws;
  }
  return nullptr;
}

const UiWorkspace* findWorkspace(const UiWorkspaces& s, std::string_view name) {
  for (const auto& ws : s.items) {
    if (ws.name == name) return &ws;
  }
  return nullptr;
}

const UiWorkspace* activeWorkspace(const UiWorkspaces& s) {
  if (s.active.empty()) return nullptr;
  return findWorkspace(s, s.active);
}

UiWorkspace& addOrUpdateWorkspace(UiWorkspaces& s, UiWorkspace ws) {
  ws.name = sanitizeWorkspaceName(ws.name);
  if (UiWorkspace* existing = findWorkspace(s, ws.name)) {
    *existing = std::move(ws);
    return *existing;
  }
  s.items.push_back(std::move(ws));
  return s.items.back();
}

bool removeWorkspace(UiWorkspaces& s, std::string_view name) {
  auto it = std::remove_if(s.items.begin(), s.items.end(), [&](const UiWorkspace& ws) {
    return ws.name == name;
  });
  const bool removed = (it != s.items.end());
  s.items.erase(it, s.items.end());
  if (removed && s.active == name) {
    s.active.clear();
  }
  return removed;
}

bool saveUiWorkspacesToFile(const UiWorkspaces& s, const std::string& path) {
  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Warn, "UiWorkspaces: failed to open for writing: " + path);
    return false;
  }

  f << "StellarForgeUiWorkspaces " << s.version << "\n";
  f << "autoSaveOnExit " << (s.autoSaveOnExit ? 1 : 0) << "\n";
  if (!s.active.empty()) {
    f << "active " << sanitizeWorkspaceName(s.active) << "\n";
  }

  // Deterministic order for diff-friendly files.
  std::vector<std::size_t> order;
  order.reserve(s.items.size());
  for (std::size_t i = 0; i < s.items.size(); ++i) order.push_back(i);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return s.items[a].name < s.items[b].name;
  });

  for (std::size_t idx : order) {
    const UiWorkspace& ws = s.items[idx];
    if (ws.name.empty()) continue;
    f << "workspace " << sanitizeWorkspaceName(ws.name) << "\n";
    if (!ws.imguiIniFile.empty()) {
      f << "ini " << ws.imguiIniFile << "\n";
    }
    for (const auto& kv : ws.windows) {
      if (kv.first.empty()) continue;
      f << "win " << sanitizeKey(kv.first) << " " << (kv.second ? 1 : 0) << "\n";
    }
    f << "end\n";
  }
  return true;
}

static void dedupeKeepLast(UiWorkspaces& s) {
  // Keep the last workspace for each name.
  std::vector<UiWorkspace> out;
  out.reserve(s.items.size());

  for (auto& ws : s.items) {
    ws.name = sanitizeWorkspaceName(ws.name);
    if (ws.name.empty()) continue;

    bool replaced = false;
    for (auto& o : out) {
      if (o.name == ws.name) {
        o = std::move(ws);
        replaced = true;
        break;
      }
    }
    if (!replaced) out.push_back(std::move(ws));
  }

  s.items = std::move(out);
}

bool loadUiWorkspacesFromFile(const std::string& path, UiWorkspaces& out) {
  std::ifstream f(path);
  if (!f) {
    stellar::core::log(stellar::core::LogLevel::Debug, "UiWorkspaces: file not found: " + path);
    return false;
  }

  std::string header;
  int version = 0;
  if (!(f >> header >> version)) {
    stellar::core::log(stellar::core::LogLevel::Warn, "UiWorkspaces: failed to read header");
    return false;
  }
  if (header != "StellarForgeUiWorkspaces") {
    stellar::core::log(stellar::core::LogLevel::Warn, "UiWorkspaces: bad header");
    return false;
  }

  UiWorkspaces s = makeDefaultUiWorkspaces();
  s.version = version;

  std::string line;
  std::getline(f, line); // consume remainder

  UiWorkspace* cur = nullptr;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream ss(line);
    std::string key;
    if (!(ss >> key)) continue;
    key = lowerAscii(key);

    if (key == "autosaveonexit") {
      int v = 1;
      ss >> v;
      s.autoSaveOnExit = (v != 0);
      continue;
    }

    if (key == "active") {
      std::string rest;
      std::getline(ss, rest);
      s.active = sanitizeWorkspaceName(rest);
      continue;
    }

    if (key == "workspace") {
      std::string rest;
      std::getline(ss, rest);
      UiWorkspace ws;
      ws.name = sanitizeWorkspaceName(rest);
      s.items.push_back(std::move(ws));
      cur = &s.items.back();
      continue;
    }

    if (key == "end") {
      cur = nullptr;
      continue;
    }

    if (!cur) {
      // Unknown line outside a workspace block: ignore.
      continue;
    }

    if (key == "ini" || key == "inifile" || key == "imguiinifile") {
      ss >> cur->imguiIniFile;
      continue;
    }

    if (key == "win" || key == "window") {
      std::string winKey;
      int v = 1;
      if (!(ss >> winKey >> v)) continue;
      winKey = sanitizeKey(std::move(winKey));
      if (winKey.empty()) continue;
      cur->windows[winKey] = (v != 0);
      continue;
    }

    // Unknown key inside workspace: ignore for forward-compat.
  }

  dedupeKeepLast(s);

  // If the active workspace doesn't exist, clear it.
  if (!s.active.empty() && !findWorkspace(s, s.active)) {
    s.active.clear();
  }

  out = std::move(s);
  return true;
}

} // namespace stellar::ui
