#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/PostFX.h"

namespace stellar::render {

// Seedable PostFX presets for quickly generating distinct visual looks.
//
// These are CPU-only helpers (safe for headless builds); the resulting settings are
// applied by PostFX::present at render time.

enum class PostFXPresetKind : core::u8 {
  Cinematic = 0,
  Retro = 1,
  Random = 2,
};

// Generate a deterministic PostFXSettings preset from a seed.
//
// - Cinematic: bloom/film-ish tuning, no retro quantization.
// - Retro: enables the Retro compositor and tunes it (pixelation, dithering, scanlines, ...).
// - Random: chooses a preset family deterministically from the seed.
PostFXSettings makeProceduralPostFXSettings(core::u64 seed, PostFXPresetKind kind = PostFXPresetKind::Random);

} // namespace stellar::render
