#include "stellar/render/ProceduralPostFX.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>

namespace stellar::render {

PostFXSettings makeProceduralPostFXSettings(core::u64 seed, PostFXPresetKind kind) {
  // Mix the seed with a stable tag to avoid accidental correlations if the caller
  // uses the same seed for other subsystems.
  core::SplitMix64 rng(core::hashCombine(seed, core::fnv1a64("postfx_preset")));

  PostFXPresetKind chosen = kind;
  if (chosen == PostFXPresetKind::Random) {
    // Deterministically pick a family.
    chosen = rng.chance(0.35) ? PostFXPresetKind::Retro : PostFXPresetKind::Cinematic;
  }

  PostFXSettings s{};
  s.enabled = true;

  // We intentionally keep warp/hyperspace at 0 so gameplay can drive those without
  // being stomped by a "generate look" button.
  s.warp = 0.0f;
  s.hyperspace = 0.0f;

  if (chosen == PostFXPresetKind::Cinematic) {
    s.retroEnabled = false;

    s.bloomEnabled = rng.chance(0.85);
    s.bloomThreshold = (float)rng.range(0.85, 1.35);
    s.bloomKnee = (float)rng.range(0.20, 1.10);
    s.bloomBoost = (float)rng.range(0.85, 1.55);
    s.bloomIntensity = (float)rng.range(0.55, 1.25);
    s.bloomPasses = rng.range(4, 12);

    s.exposure = (float)rng.range(0.75, 1.65);
    s.gamma = (float)rng.range(1.80, 2.40);

    // Auto-exposure is a natural fit for cinematic looks, but keep it optional.
    s.autoExposureEnabled = rng.chance(0.65);
    if (s.autoExposureEnabled) {
      s.autoExposureKey = (float)rng.range(0.14, 0.22);
      s.autoExposureMin = (float)rng.range(0.20, 0.65);
      s.autoExposureMax = (float)rng.range(1.80, 6.00);
      s.autoExposureSpeedUp = (float)rng.range(1.50, 5.50);
      s.autoExposureSpeedDown = (float)rng.range(0.60, 2.50);
      s.autoExposureCenterWeight = (float)rng.range(0.35, 0.90);
      s.autoExposureMaxSize = rng.range(96, 256);
    }

    s.vignette = (float)rng.range(0.05, 0.35);
    s.grain = (float)rng.range(0.0, 0.060);
    s.chromaticAberration = (float)rng.range(0.0000, 0.0030);

    // Retro fields unused, keep them in a sane neutral state.
    s.retroPixelSize = 1;
    s.retroColorSteps = 256;
    s.retroDitherStrength = 0.0f;
    s.retroScanlines = 0.0f;
    s.retroCurvature = 0.0f;
    s.retroJitter = 0.0f;
  } else {
    // --- Retro / CRT ---
    s.retroEnabled = true;

    s.bloomEnabled = rng.chance(0.55);
    s.bloomThreshold = (float)rng.range(0.65, 1.10);
    s.bloomKnee = (float)rng.range(0.30, 1.10);
    s.bloomBoost = (float)rng.range(0.90, 2.20);
    s.bloomIntensity = (float)rng.range(0.25, 1.20);
    s.bloomPasses = rng.range(2, 10);

    s.exposure = (float)rng.range(0.85, 1.55);
    s.gamma = (float)rng.range(2.00, 2.45);

    // Retro mode usually looks better with a stable exposure (avoid flicker).
    s.autoExposureEnabled = false;

    s.vignette = (float)rng.range(0.0, 0.35);
    s.grain = (float)rng.range(0.02, 0.12);
    s.chromaticAberration = (float)rng.range(0.0, 0.0040);

    s.retroPixelSize = rng.range(2, 6);

    // Pick from a small set of aesthetically pleasing step counts.
    constexpr int kStepOptions[] = {8, 12, 16, 24, 32, 48, 64};
    constexpr int kStepCount = (int)(sizeof(kStepOptions) / sizeof(kStepOptions[0]));
    s.retroColorSteps = kStepOptions[rng.range(0, kStepCount - 1)];

    s.retroDitherStrength = (float)rng.range(0.20, 1.00);
    s.retroScanlines = (float)rng.range(0.0, 0.85);
    s.retroCurvature = (float)rng.range(0.0, 0.70);
    s.retroJitter = (float)rng.range(0.0, 0.35);
  }

  // Clamp safety.
  s.bloomPasses = std::clamp(s.bloomPasses, 0, 16);
  s.autoExposureMaxSize = std::clamp(s.autoExposureMaxSize, 32, 512);
  s.retroPixelSize = std::clamp(s.retroPixelSize, 1, 64);
  s.retroColorSteps = std::clamp(s.retroColorSteps, 2, 256);
  s.retroDitherStrength = std::clamp(s.retroDitherStrength, 0.0f, 1.0f);
  s.retroScanlines = std::clamp(s.retroScanlines, 0.0f, 1.0f);
  s.retroCurvature = std::clamp(s.retroCurvature, 0.0f, 1.0f);
  s.retroJitter = std::clamp(s.retroJitter, 0.0f, 1.0f);

  return s;
}

} // namespace stellar::render
