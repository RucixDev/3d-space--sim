#pragma once

#include <functional>
#include <string>

namespace stellar::game {

// Minimal in-game CVar browser/editor.
//
// This lives in the game app (not the core library) to avoid pulling ImGui into
// headless builds. The underlying CVar system is in stellar::core.

struct CVarWindowState {
  bool open{false};
  bool focusFilter{false};

  // Search/filter
  char filter[128]{};
  bool sortByRelevance{true};
  bool showChangedOnly{false};
  bool showArchivedOnly{false};
  bool showCheats{false};

  // Quick set line: "name = value" (press Enter)
  char setLine[256]{};

  // Persistence helpers (load/save file)
  char configPath[128]{"cvars.cfg"};

  // String edit modal
  bool requestEditPopup{false};
  std::string editName;
  char editBuf[512]{};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

// Open the window (optionally focusing the filter box).
void openCVarWindow(CVarWindowState& st, bool focusFilter = true);

// Draw the window. (No-op if st.open == false.)
void drawCVarWindow(CVarWindowState& st, const ToastFn& toast);

} // namespace stellar::game
