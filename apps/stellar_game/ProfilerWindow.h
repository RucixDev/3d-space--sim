#pragma once

#include "stellar/core/Profiler.h"

#include <functional>
#include <string>

namespace stellar::game {

// In-game CPU profiler viewer.
//
// Notes:
//  - The profiler capture itself lives in stellar::core (headless).
//  - This window is only compiled for the SDL/OpenGL + ImGui app.

struct ProfilerWindowState {
  bool open{false};

  bool pause{false};
  bool showFlameGraph{true};
  bool showAggregates{true};
  bool showPlot{true};

  // 0 = newest, 1 = previous, ... (only used when pause=true)
  int selectedFrameOffset{0};

  // Simple substring filter applied to aggregate rows.
  char filter[96]{};

  // ---- Trace export (chrome://tracing / Perfetto legacy JSON) ----
  // Default export path is relative to the working directory.
  char exportPath[256]{"profiler_trace.json"};
  bool exportAllFrames{true};
  bool exportIncludeFrameEvents{true};
  bool exportPretty{false};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawProfilerWindow(ProfilerWindowState& st,
                        core::Profiler& profiler,
                        const ToastFn& toast);

} // namespace stellar::game
