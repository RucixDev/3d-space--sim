#pragma once

#include <functional>
#include <string>

namespace stellar::game {

// Simple "About / Build Info" window.
//
// Purpose:
//  - Make bug reports actionable by exposing versions + feature flags.
//  - Provide a one-click "copy to clipboard" summary.
//
// This is intentionally a pure ImGui window (only compiled in the SDL/OpenGL app).

struct BuildInfoWindowState {
  bool open{false};

  // UI toggles
  bool showRuntimeGlInfo{true};
  bool showCompileTimeInfo{true};
  bool showFeatureFlags{true};
  bool showDependencyVersions{true};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawBuildInfoWindow(BuildInfoWindowState& st, const ToastFn& toast);

} // namespace stellar::game
