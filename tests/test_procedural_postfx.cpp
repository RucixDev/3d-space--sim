#include "stellar/render/ProceduralPostFX.h"

#include <iostream>

int test_procedural_postfx() {
  int fails = 0;

  const stellar::core::u64 seed = 123456789ull;

  // Determinism: same seed + kind should generate identical settings.
  {
    const auto a = stellar::render::makeProceduralPostFXSettings(seed, stellar::render::PostFXPresetKind::Cinematic);
    const auto b = stellar::render::makeProceduralPostFXSettings(seed, stellar::render::PostFXPresetKind::Cinematic);

    if (a.bloomThreshold != b.bloomThreshold ||
        a.bloomKnee != b.bloomKnee ||
        a.bloomBoost != b.bloomBoost ||
        a.bloomIntensity != b.bloomIntensity ||
        a.bloomPasses != b.bloomPasses ||
        a.exposure != b.exposure ||
        a.gamma != b.gamma ||
        a.autoExposureEnabled != b.autoExposureEnabled ||
        a.autoExposureKey != b.autoExposureKey ||
        a.autoExposureMin != b.autoExposureMin ||
        a.autoExposureMax != b.autoExposureMax ||
        a.autoExposureSpeedUp != b.autoExposureSpeedUp ||
        a.autoExposureSpeedDown != b.autoExposureSpeedDown ||
        a.autoExposureCenterWeight != b.autoExposureCenterWeight ||
        a.autoExposureMaxSize != b.autoExposureMaxSize ||
        a.vignette != b.vignette ||
        a.grain != b.grain ||
        a.chromaticAberration != b.chromaticAberration ||
        a.retroEnabled != b.retroEnabled) {
      std::cerr << "[test_procedural_postfx] not deterministic for Cinematic preset\n";
      ++fails;
    }
  }

  // Sanity / range checks.
  {
    const auto s = stellar::render::makeProceduralPostFXSettings(seed, stellar::render::PostFXPresetKind::Retro);

    if (!s.enabled || !s.retroEnabled) {
      std::cerr << "[test_procedural_postfx] Retro preset should be enabled + retroEnabled\n";
      ++fails;
    }

    if (s.retroPixelSize < 1 || s.retroPixelSize > 64) {
      std::cerr << "[test_procedural_postfx] retroPixelSize out of range: " << s.retroPixelSize << "\n";
      ++fails;
    }

    if (s.retroColorSteps < 2 || s.retroColorSteps > 256) {
      std::cerr << "[test_procedural_postfx] retroColorSteps out of range: " << s.retroColorSteps << "\n";
      ++fails;
    }

    if (s.retroDitherStrength < 0.0f || s.retroDitherStrength > 1.0f) {
      std::cerr << "[test_procedural_postfx] retroDitherStrength out of range\n";
      ++fails;
    }

    if (s.bloomPasses < 0 || s.bloomPasses > 16) {
      std::cerr << "[test_procedural_postfx] bloomPasses out of range: " << s.bloomPasses << "\n";
      ++fails;
    }

    if (s.warp != 0.0f || s.hyperspace != 0.0f) {
      std::cerr << "[test_procedural_postfx] presets should not force warp/hyperspace\n";
      ++fails;
    }
  }

  // Different seeds should usually change something obvious.
  {
    const auto a = stellar::render::makeProceduralPostFXSettings(seed, stellar::render::PostFXPresetKind::Cinematic);
    const auto b = stellar::render::makeProceduralPostFXSettings(seed + 1, stellar::render::PostFXPresetKind::Cinematic);
    if (a.bloomThreshold == b.bloomThreshold && a.exposure == b.exposure && a.vignette == b.vignette) {
      std::cerr << "[test_procedural_postfx] adjacent seeds produced identical core params (unexpected)\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_procedural_postfx] pass\n";
  return fails;
}
