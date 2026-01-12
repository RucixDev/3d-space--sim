#include "stellar/render/ProceduralCompositor.h"

#include <catch2/catch_test_macros.hpp>

#include <string>

TEST_CASE("procedural compositor generates deterministic defines", "[proc][postfx]") {
  using namespace stellar;
  const core::u64 seed = 0x1234ULL;

  for (int k = 0; k < 5; ++k) {
    const render::CompositorKind kind = (render::CompositorKind)k;
    const render::CompositorRecipe r1 = render::makeProceduralCompositorRecipe(seed, kind);
    const render::CompositorRecipe r2 = render::makeProceduralCompositorRecipe(seed, kind);

    const std::string d1 = render::buildCompositorShaderDefines(r1);
    const std::string d2 = render::buildCompositorShaderDefines(r2);

    REQUIRE(d1 == d2);
    REQUIRE(d1.find("STELLAR_TONEMAP_MODE") != std::string::npos);
    REQUIRE(d1.find("STELLAR_GRADE_ENABLED") != std::string::npos);
    REQUIRE(d1.find("STELLAR_VIGNETTE_MUL") != std::string::npos);
    REQUIRE(d1.find("STELLAR_GRAIN_MUL") != std::string::npos);
  }
}

TEST_CASE("procedural compositor produces stable identity-like ranges", "[proc][postfx]") {
  using namespace stellar;
  const core::u64 seed = 0xBADC0FFEEULL;
  const auto r = render::makeProceduralCompositorRecipe(seed, render::CompositorKind::Cinematic);

  // These are sanity checks so recipes don't explode.
  REQUIRE(r.saturation >= 0.0f);
  REQUIRE(r.saturation <= 2.0f);
  REQUIRE(r.contrast >= 0.0f);
  REQUIRE(r.contrast <= 2.0f);
  REQUIRE(r.vignetteMul >= 0.0f);
  REQUIRE(r.vignetteMul <= 3.0f);
  REQUIRE(r.grainMul >= 0.0f);
  REQUIRE(r.grainMul <= 3.0f);
}
