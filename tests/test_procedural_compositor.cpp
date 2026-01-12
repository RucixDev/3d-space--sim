#include "stellar/render/ProceduralCompositor.h"

#include <iostream>
#include <string>

int test_procedural_compositor() {
  using namespace stellar;

  int fails = 0;

  auto require = [&](bool cond, const char* msg) {
    if (!cond) {
      std::cerr << "[test_procedural_compositor] " << msg << "\n";
      ++fails;
    }
  };

  auto requireContains = [&](const std::string& s, const char* needle) {
    if (s.find(needle) == std::string::npos) {
      std::cerr << "[test_procedural_compositor] missing token: " << needle << "\n";
      ++fails;
    }
  };

  // Determinism: same seed + kind should generate identical define strings.
  {
    const core::u64 seed = 0x1234ULL;
    for (int k = 0; k < 5; ++k) {
      const render::CompositorKind kind = (render::CompositorKind)k;

      const render::CompositorRecipe r1 = render::makeProceduralCompositorRecipe(seed, kind);
      const render::CompositorRecipe r2 = render::makeProceduralCompositorRecipe(seed, kind);

      const std::string d1 = render::buildCompositorShaderDefines(r1);
      const std::string d2 = render::buildCompositorShaderDefines(r2);

      require(d1 == d2, "non-deterministic compositor defines for identical seed/kind");

      // Basic sanity: generated defines should include the core knobs we expect.
      requireContains(d1, "STELLAR_TONEMAP_MODE");
      requireContains(d1, "STELLAR_GRADE_ENABLED");
      requireContains(d1, "STELLAR_VIGNETTE_MUL");
      requireContains(d1, "STELLAR_GRAIN_MUL");
    }
  }

  // Sanity/range checks: recipes should stay within reasonable, identity-like bounds.
  {
    const core::u64 seed = 0xBADC0FFEEULL;
    const auto r = render::makeProceduralCompositorRecipe(seed, render::CompositorKind::Cinematic);

    require(r.saturation >= 0.0f, "saturation < 0");
    require(r.saturation <= 2.0f, "saturation > 2");
    require(r.contrast >= 0.0f, "contrast < 0");
    require(r.contrast <= 2.0f, "contrast > 2");
    require(r.vignetteMul >= 0.0f, "vignetteMul < 0");
    require(r.vignetteMul <= 3.0f, "vignetteMul > 3");
    require(r.grainMul >= 0.0f, "grainMul < 0");
    require(r.grainMul <= 3.0f, "grainMul > 3");
  }

  if (fails == 0) std::cout << "[test_procedural_compositor] pass\n";
  return fails;
}
