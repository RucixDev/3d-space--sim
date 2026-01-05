#pragma once

#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::ui {

struct UiWorkspace; // from UiWorkspaces.h

// -----------------------------------------------------------------------------
// Window Registry (headless)
// -----------------------------------------------------------------------------
//
// The game prototype historically stored window open/close state as a loose set
// of booleans in the main loop, which made it easy to forget to wire new windows
// into:
//   - workspace capture/apply
//   - command palette
//   - bulk "open all" / "close all" actions
//
// This registry is an opt-in bridge:
//   - it is still headless (no ImGui dependency)
//   - it can store pointers to runtime bools
//   - it can fire small callbacks on open/close (focus input, etc.)

struct WindowDesc {
  // Stable, file-friendly key (used by UiWorkspaces).
  std::string key;

  // Human-readable name (what the player sees).
  std::string label;

  // Optional grouping hint used by the app (Window Manager UI, etc.).
  std::string group;

  // Higher = appears earlier when query is empty (command palette).
  int palettePriority{0};

  // Default open state used by resetToDefaults().
  bool defaultOpen{false};

  // If false, this window won't be stored in UiWorkspace snapshots.
  bool persistInWorkspace{true};
};

struct WindowBinding {
  WindowDesc desc;

  // Pointer to the runtime open/closed flag.
  bool* open{nullptr};

  // Optional hooks fired when the open state changes.
  std::function<void()> onOpened{};
  std::function<void()> onClosed{};

  // Optional key-hint provider (e.g. "Ctrl+P").
  // The provider is used only by UI that wants to show the hint; it does not
  // participate in matching.
  std::function<std::string()> shortcutLabel{};
};

class WindowRegistry {
 public:
  WindowRegistry() = default;

  // Add or replace a binding by key. Returns the stored binding.
  WindowBinding& add(WindowBinding b);

  // Lookup (linear scan; the registry is small).
  WindowBinding* find(std::string_view key);
  const WindowBinding* find(std::string_view key) const;

  const std::vector<WindowBinding>& items() const { return items_; }

  // Set a window open/closed. Returns true if the state changed.
  bool setOpen(std::string_view key, bool open, bool fireCallbacks = true);

  // Bulk actions.
  void setAll(bool open, bool fireCallbacks = true);
  void resetToDefaults(bool fireCallbacks = true);

 private:
  std::vector<WindowBinding> items_;
};

// Capture/apply helpers for UiWorkspaces.
// These helpers only touch windows where desc.persistInWorkspace is true.
void captureWorkspaceWindows(const WindowRegistry& reg, UiWorkspace& ws);
void applyWorkspaceWindows(WindowRegistry& reg, const UiWorkspace& ws, bool fireCallbacks = true);

} // namespace stellar::ui
