#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/ProceduralGraph.h"

#include <functional>
#include <string>

namespace stellar::game {

struct ProceduralLabWindowState {
  bool open{false};

  // Quickstart: start from a preset, then tweak nodes/palette.
  render::ProcGraphPreset preset{render::ProcGraphPreset::Nebula};
  bool lockToPreset{true};

  // Bake settings
  core::u64 seed{0xC0FFEEULL};
  int resolution{512};
  bool autoBake{true};

  // Bake quality
  bool bakeGenerateMips{true};
  float bakeDitherStrength{1.0f};
  bool bakePackHeightInAlpha{false};

  // Output options
  bool usePalette{true};
  bool showShaderSource{false};

  // Export
  char exportPath[256]{"screenshots/procedural_lab.png"};
  bool exportFlipY{true};

  // Graph I/O
  // Saved graphs are small, versioned text files for sharing/reuse.
  char graphPath[256]{"proc_graphs/nebula.procgraph"};

  // Graph data
  render::ProcGraph graph{render::ProcGraph::makeDefault()};
  bool dirty{true};

  // GPU baker + last error
  render::ProceduralGraphBaker baker{};
  std::string lastError{};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawProceduralLabWindow(ProceduralLabWindowState& st, float timeSec, const ToastFn& toast);

} // namespace stellar::game
