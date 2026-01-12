#pragma once

#include "stellar/core/Types.h"

#include <string>

namespace stellar::render {

// Procedural "compositor shader" generator.
//
// This module generates deterministic PostFX composite shader *variants* from a
// seed by emitting a small set of `#define` lines. The defines are injected into
// PostFX's built-in composite fragment shader (see PostFX::buildCompositeVariant).
//
// Why defines instead of uniforms?
// - This is intentionally "engine-ish": it allows the game to build and cache
//   custom shader permutations without needing a full material editor.
// - The generated looks are stable per-seed and can do things that are awkward
//   to express purely via runtime parameters.

enum class CompositorKind : int {
  Cinematic = 0,
  Noir,
  Vaporwave,
  Analog,
  Random
};

enum class TonemapMode : int {
  Aces = 0,
  Reinhard,
  Uncharted2,
  Clamp
};

struct CompositorRecipe {
  TonemapMode tonemap{TonemapMode::Aces};

  // Basic grading controls (applied after tonemap and before output gamma).
  float saturation{1.0f};
  float contrast{1.0f};
  float hueShiftRad{0.0f};

  // Lift/Gamma/Gain (per-channel). Defaults are identity.
  float liftR{0.0f}, liftG{0.0f}, liftB{0.0f};
  float gammaR{1.0f}, gammaG{1.0f}, gammaB{1.0f};
  float gainR{1.0f}, gainG{1.0f}, gainB{1.0f};

  // Split-toning (tints shadows/highlights). If disabled, strength is ignored.
  bool splitTone{false};
  float splitStrength{0.0f};
  float shadowTintR{1.0f}, shadowTintG{1.0f}, shadowTintB{1.0f};
  float highlightTintR{1.0f}, highlightTintG{1.0f}, highlightTintB{1.0f};

  // Multipliers applied to existing PostFX uniforms (vignette/grain).
  float vignetteMul{1.0f};
  float grainMul{1.0f};
};

// Generate a deterministic recipe from (seed, kind).
CompositorRecipe makeProceduralCompositorRecipe(core::u64 seed, CompositorKind kind);

// Emit `#define ...` lines (no #version). This string is meant to be passed to:
//   PostFX::buildCompositeVariant(defines)
//
// The output is stable for a given recipe (deterministic).
std::string buildCompositorShaderDefines(const CompositorRecipe& r);

const char* compositorKindName(CompositorKind kind);

} // namespace stellar::render
