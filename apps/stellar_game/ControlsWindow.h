#pragma once

#include "ControlsConfig.h"

#include <SDL.h>

#include <functional>
#include <string>

namespace stellar::game {

// A richer in-game keybind editor with:
//  - fuzzy search
//  - conflict detection
//  - profile management (multiple controls files)
//
// This stays in the game app because it depends on SDL scancodes.
struct ControlsWindowState {
  bool open{false};

  // Search/filter
  char filter[128]{};
  bool focusFilter{false};
  bool showChangedOnly{false};
  bool showConflictsOnly{false};
  bool sortByRelevance{true};

  // Profile manager UI state
  int selectedProfile{0};
  char newProfileName[64]{"profile"};

  // Rebind capture state.
  enum class RebindKind : int {
    None = 0,
    Action,
    Hold,
    AxisPositive,
    AxisNegative,
  };

  struct RebindRequest {
    RebindKind kind{RebindKind::None};
    std::string label;
    KeyChord* chord{nullptr};
    SDL_Scancode* hold{nullptr};
    AxisPair* axis{nullptr};
  };

  RebindRequest rebind{};

  // If set, the controls window will attempt to scroll the matching row into view.
  std::string scrollToLabel;
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

// Open the window (and optionally focus the search box).
void openControlsWindow(ControlsWindowState& st, bool focusFilter = true);

// Handle a KEYDOWN while a rebind capture is active.
// Returns true if the event was consumed by the capture.
bool handleControlsRebindKeydown(ControlsWindowState& st,
                                const SDL_KeyboardEvent& ev,
                                ControlsConfig& controls,
                                bool& controlsDirty,
                                const ToastFn& toast);

// Draw the window.
void drawControlsWindow(ControlsWindowState& st,
                        ControlsConfig& controls,
                        const ControlsConfig& defaults,
                        bool& controlsDirty,
                        bool& autoSaveOnExit,
                        std::string& controlsPath,
                        const ToastFn& toast);

} // namespace stellar::game
