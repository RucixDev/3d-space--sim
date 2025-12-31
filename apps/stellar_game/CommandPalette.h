#pragma once

#include <functional>
#include <string>
#include <string_view>
#include <vector>

// Command Palette UI component (ImGui).
//
// Lives in the game app (not the core library) to avoid pulling ImGui into
// headless builds.

namespace stellar::game {

struct PaletteItem {
  std::string label;
  std::string detail;   // optional secondary text
  std::string shortcut; // optional key hint
  int priority{0};       // higher = shown first when query is empty
  std::function<void()> action;
};

struct CommandPaletteState {
  bool open{false};
  bool justOpened{false};
  bool focusInput{false};

  // Fixed buffer to keep ImGui integration simple (no imgui_stdlib dependency).
  char query[256]{};
  int selected{0};
};

// Open the palette. If resetQuery is true, clears the filter text.
void openCommandPalette(CommandPaletteState& st, bool resetQuery = true);

// Draw the palette. Returns true if a command was executed this frame.
//
// The palette does not own items; they can be rebuilt each frame.
bool drawCommandPalette(CommandPaletteState& st, const std::vector<PaletteItem>& items);

// Utility: a small fuzzy matcher used by the palette (exposed for potential
// reuse in other app UI).
int fuzzyMatchScore(std::string_view query, std::string_view text);

} // namespace stellar::game
