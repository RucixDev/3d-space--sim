#pragma once

#include <functional>
#include <string>

namespace stellar::game {

// A lightweight "photo mode" / screenshot capture panel.
//
// The actual screenshot capture happens inside the main loop so we can
// capture either:
//  - the world render before any UI is drawn, or
//  - the final frame including UI overlays.
struct PhotoModeWindowState {
  bool open{false};

  // If true, the sim is paused while the window is open (handy for framing).
  bool pauseWhileOpen{true};

  // If true, capture after UI rendering (includes HUD/windows).
  // If false, capture before UI is drawn (world-only).
  bool captureIncludeUi{false};

  // If true, append a timestamp to the output filename.
  bool timestamp{true};

  // If true, copy the saved file path to clipboard on success.
  bool copyPathToClipboard{true};

  // Output directory and base filename.
  char outDir[160]{"screenshots"};
  char baseName[64]{"stellarforge"};

  // Last saved file path (for quick copy).
  std::string lastSavedPath;

  // One-shot focus helpers.
  bool focusDir{false};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

struct PhotoModeContext {
  int framebufferW{0};
  int framebufferH{0};

  // Called when the user presses "Take screenshot".
  // includeUi determines whether the capture should happen after UI rendering.
  std::function<void(bool includeUi,
                     const std::string& outDir,
                     const std::string& baseName,
                     bool timestamp,
                     bool copyPathToClipboard)> requestScreenshot;

  ToastFn toast;

  // Current pause state (to display).
  bool paused{false};
};

void drawPhotoModeWindow(PhotoModeWindowState& st, const PhotoModeContext& ctx);

} // namespace stellar::game
