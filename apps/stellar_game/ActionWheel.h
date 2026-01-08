#pragma once

#include <functional>
#include <string>
#include <vector>

#include <imgui.h>

namespace stellar::game {

// A quick-access radial menu overlay (Dear ImGui).
//
// Intended usage:
//  - openActionWheel(state) on keydown
//  - while state.open, build a vector<ActionWheelItem> for the current context
//  - call drawActionWheelOverlay(state, items, pageName, pageIndex, pageCount)
//  - close/confirm from the app (e.g. on key release) using state.hovered
//
// NOTE: This is app-local (not in the core library) so headless builds don't depend on ImGui.

struct ActionWheelItem {
  std::string label;
  std::string detail;    // Optional, shown in the center while highlighted
  std::string shortcut;  // Optional key hint (e.g. "F1")
  bool enabled{true};
  std::function<void()> action;
};

struct ActionWheelState {
  bool open{false};
  bool justOpened{false};
  bool requestClose{false};
  int hovered{-1};

  ImVec2 center{0.0f, 0.0f};

  // Layout
  float innerRadius{48.0f};
  float outerRadius{170.0f};
  float deadzoneRadius{22.0f};
  float labelRadius{118.0f};

  // Optional: the overlay can display a page indicator.
  int page{0};

  // Animation state
  float openT{0.0f};
};

void openActionWheel(ActionWheelState& st);
void closeActionWheel(ActionWheelState& st);

// Draw the overlay and update st.hovered.
// Returns the index of the confirmed item when the user clicks (left mouse release), or -1.
int drawActionWheelOverlay(ActionWheelState& st,
                          const std::vector<ActionWheelItem>& items,
                          const char* pageName = nullptr,
                          int pageIndex = 0,
                          int pageCount = 0);

} // namespace stellar::game
