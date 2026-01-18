#include "test_harness.h"

#include "stellar/ui/TextFx.h"

#include <cmath>
#include <string_view>

using namespace stellar;

static int test_textfx_compile_basic() {
  int failures = 0;

  const ui::textfx::Program p = ui::textfx::compile(
    "[color #ff0000]Hi[/color] [wave amp=2 freq=1 speed=0]X[/wave]"
  );

  CHECK(p.plain == "Hi X");
  CHECK(p.spans.size() == 2);

  // Expect first span to cover "Hi" (2 glyphs)
  CHECK((int)p.spans[0].startGlyph == 0);
  CHECK((int)p.spans[0].endGlyph == 2);

  // Second span covers "X" (after space)
  CHECK((int)p.spans[1].startGlyph == 3);
  CHECK((int)p.spans[1].endGlyph == 4);

  return failures;
}

static int test_textfx_unknown_tags_preserved() {
  int failures = 0;

  const ui::textfx::Program p = ui::textfx::compile(
    "Docking: [CONTRABAND] lanes detected"
  );

  CHECK(p.plain == "Docking: [CONTRABAND] lanes detected");
  CHECK(p.spans.empty());

  return failures;
}


static int test_textfx_program_glyph_count() {
  int failures = 0;

  // Mix ASCII + UTF-8 to ensure we count codepoints, not bytes.
  const ui::textfx::Program p = ui::textfx::compile(
    "[color #ff0000]Hi[/color] — [wave amp=2 freq=1 speed=0]π[/wave]"
  );

  CHECK(p.glyphCount == ui::textfx::glyphCountUtf8(p.plain));
  CHECK(p.plain == "Hi — π");
  CHECK(p.glyphCount == 6); // H i (space) — (space) π

  return failures;
}

static int test_textfx_program_cache() {
  int failures = 0;

  ui::textfx::ProgramCache cache(2);

  const auto& p1 = cache.get("[color #00ff00]OK[/color]");
  const auto& p2 = cache.get("[color #00ff00]OK[/color]");

  CHECK(p1.plain == "OK");
  CHECK(&p1 == &p2);
  CHECK(cache.size() == 1);

  // Fill and evict.
  cache.get("A");
  cache.get("B");
  CHECK(cache.size() == 2);

  // Ensure the original string still compiles correctly after possible eviction.
  const auto& p3 = cache.get("[color #00ff00]OK[/color]");
  CHECK(p3.plain == "OK");

  return failures;
}

static int test_textfx_eval_wave_and_pulse() {
  int failures = 0;

  // wave with speed=0 and phase=0 should be deterministic per glyphIndex.
  // pulse with speed=0 gives constant alpha at sin(phase) => 0.
  const ui::textfx::Program p = ui::textfx::compile(
    "[wave amp=10 freq=0.25 speed=0][pulse min=0.2 max=1 speed=0]ABCD[/pulse][/wave]"
  );

  ui::textfx::EvalStyle st;
  st.baseColor = ui::Color4f{1, 1, 1, 1};
  st.seed = 123;

  // glyphIndex 0, sin(0) => dy=0
  {
    const auto g0 = ui::textfx::evalGlyph(p, 0, 0.0f, 0.0f, 'A', 5.0f, st);
    CHECK(std::fabs(g0.dy - 0.0f) < 1e-4f);
    // pulse: speed=0 => ph=0 => sin(0)=0 => s01=0.5 => alpha=(0.2+1)/2=0.6
    CHECK(std::fabs(g0.color.a - 0.6f) < 1e-4f);
  }

  // glyphIndex 1, sin(2pi*0.25)=sin(pi/2)=1 => dy=+10
  {
    const auto g1 = ui::textfx::evalGlyph(p, 1, 0.0f, 0.0f, 'B', 5.0f, st);
    CHECK(std::fabs(g1.dy - 10.0f) < 1e-3f);
  }

  return failures;
}

int test_text_fx() {
  int failures = 0;

  failures += test_textfx_compile_basic();
  failures += test_textfx_unknown_tags_preserved();
  failures += test_textfx_program_glyph_count();
  failures += test_textfx_program_cache();
  failures += test_textfx_eval_wave_and_pulse();
  return failures;
}
