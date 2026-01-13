#pragma once

#include "stellar/proc/MiePhase.h"
#include "stellar/proc/PhaseMultipleScattering.h"
#include "stellar/render/Texture.h"

#include <functional>
#include <future>
#include <string>
#include <vector>

namespace stellar::game {

struct SpectralMieLabJobOutput {
  proc::SpectralMiePhaseResult result{};
  std::vector<std::uint8_t> rgba{};
  double ms{0.0};
  std::string error{};
};

struct SpectralMieLabWindowState {
  bool open{false};

  // Generation settings.
  proc::SpectralMiePhaseSettings settings{};
  bool autoRegenerate{true};

  // When enabled, the main renderer can opt into sampling the generated LUT
  // (hooked up in main.cpp / World Visuals).
  bool applyToAtmospheres{false};
  float atmosphereMieStrength{1.0f};

  // Analytic multiple scattering phase LUT (derived from the single-scatter
  // LUT via a Legendre moment / geometric-series model).
  bool enableMultipleScattering{true};
  proc::MultipleScatteringPhaseSettings msSettings{};

  // Optional runtime hookup for volumetric atmospheres.
  bool applyMsToVolumetric{false};
  float atmosphereMsPhaseStrength{1.0f};

  // Async job management.
  bool dirty{true};
  bool jobRunning{false};
  bool jobQueued{false};
  proc::SpectralMiePhaseSettings queuedSettings{};
  std::future<SpectralMieLabJobOutput> jobFuture{};

  // Latest results.
  proc::SpectralMiePhaseResult result{};
  std::vector<std::uint8_t> rgba{};
  double lastGenMs{0.0};
  std::string lastError{};

  // Multiple scattering derived results.
  proc::SpectralMiePhaseResult msResult{};
  std::vector<std::uint8_t> msRgba{};
  double lastMsGenMs{0.0};
  std::string lastMsError{};

  // GPU texture (height=1 RGBA8) holding phase curves.
  render::Texture2D lutTex{};
  render::Texture2D msLutTex{};

  // Export controls.
  char exportDir[256]{"screenshots"};
  char exportBaseName[256]{"mie_phase_lut"};
  bool exportTimestamp{true};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawSpectralMieLabWindow(SpectralMieLabWindowState& state,
                             float timeSec,
                             const ToastFn& toast);

} // namespace stellar::game
