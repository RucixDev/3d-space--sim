#pragma once

#include "stellar/ui/LogBuffer.h"

#include <cstdint>
#include <functional>
#include <string>

namespace stellar::game {

// In-game log viewer.
//
// The core library already writes log lines to stderr, but an in-game view is
// useful for:
//  - debugging control bindings and simulation behavior
//  - diagnosing problems on systems where stderr is hard to access (Windows)
//  - quickly copying log snippets for bug reports
//
// This window consumes a thread-safe ui::LogBuffer that is fed via a
// core::LogSink.

struct LogWindowState {
  bool open{false};

  // Live filter (fuzzy). Applied to the message string.
  char filter[128]{};
  bool sortByRelevance{false};
  bool autoScroll{true};

  // Per-level toggles.
  bool showTrace{false};
  bool showDebug{false};
  bool showInfo{true};
  bool showWarn{true};
  bool showError{true};

  // Current row selection.
  std::uint64_t selectedSeq{0};

  // Export path UI.
  char exportPath[200]{"logs/session.log"};

  // One-shot focus helpers.
  bool focusFilter{false};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawLogWindow(LogWindowState& st,
                   ui::LogBuffer& buffer,
                   const ToastFn& toast);

} // namespace stellar::game
