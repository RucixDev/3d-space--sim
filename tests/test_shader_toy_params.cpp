#include "stellar/render/ShaderToyParams.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

int test_shader_toy_params() {
  int failures = 0;

  render::ShaderToyParamSet ps{};

  // A clamped float param.
  {
    render::ShaderToyParamDef d{};
    d.type = render::ShaderToyParamType::Float;
    d.name = "a";
    d.defaultValue = {0.5f, 0, 0, 0};
    d.minValue = {0.0f, 0, 0, 0};
    d.maxValue = {1.0f, 0, 0, 0};
    ps.defs.push_back(d);
  }

  // A vec2 param.
  {
    render::ShaderToyParamDef d{};
    d.type = render::ShaderToyParamType::Vec2;
    d.name = "uv";
    d.defaultValue = {0.0f, 0.0f, 0.0f, 0.0f};
    d.minValue = {-1.0f, -1.0f, 0.0f, 0.0f};
    d.maxValue = {1.0f, 1.0f, 0.0f, 0.0f};
    ps.defs.push_back(d);
  }

  // A bool param.
  {
    render::ShaderToyParamDef d{};
    d.type = render::ShaderToyParamType::Bool;
    d.name = "flag";
    d.defaultValue = {0.0f, 0, 0, 0};
    d.minValue = {0.0f, 0, 0, 0};
    d.maxValue = {1.0f, 0, 0, 0};
    ps.defs.push_back(d);
  }

  ps.resetToDefaults();
  ps.rebuildIndex();

  const int idxA = ps.findIndex("a");
  const int idxUv = ps.findIndex("uv");
  const int idxFlag = ps.findIndex("flag");

  CHECK(idxA == 0);
  CHECK(idxUv == 1);
  CHECK(idxFlag == 2);

  // --- setValue clamps ---
  {
    CHECK(ps.setValue("a", {2.0f, 9.0f, 9.0f, 9.0f}));
    CHECK(std::abs(ps.values[idxA][0] - 1.0f) < 1.0e-6f);
    CHECK(ps.values[idxA][1] == 0.0f);
    CHECK(ps.values[idxA][2] == 0.0f);
    CHECK(ps.values[idxA][3] == 0.0f);

    CHECK(ps.setValue(idxA, {-1.0f, 1.0f, 1.0f, 1.0f}));
    CHECK(std::abs(ps.values[idxA][0] - 0.0f) < 1.0e-6f);
    CHECK(ps.values[idxA][1] == 0.0f);
    CHECK(ps.values[idxA][2] == 0.0f);
    CHECK(ps.values[idxA][3] == 0.0f);
  }

  // --- setRawValue does NOT clamp but still sanitizes unused lanes ---
  {
    CHECK(ps.setRawValue("a", {123.0f, 5.0f, 6.0f, 7.0f}));
    CHECK(std::abs(ps.values[idxA][0] - 123.0f) < 1.0e-6f);
    CHECK(ps.values[idxA][1] == 0.0f);
    CHECK(ps.values[idxA][2] == 0.0f);
    CHECK(ps.values[idxA][3] == 0.0f);

    CHECK(ps.setRawValue(idxUv, {2.0f, -2.0f, 9.0f, 9.0f}));
    CHECK(std::abs(ps.values[idxUv][0] - 2.0f) < 1.0e-6f);
    CHECK(std::abs(ps.values[idxUv][1] - -2.0f) < 1.0e-6f);
    CHECK(ps.values[idxUv][2] == 0.0f);
    CHECK(ps.values[idxUv][3] == 0.0f);

    CHECK(ps.setRawValue(idxFlag, {0.75f, 1.0f, 2.0f, 3.0f}));
    CHECK(std::abs(ps.values[idxFlag][0] - 0.75f) < 1.0e-6f);
    CHECK(ps.values[idxFlag][1] == 0.0f);
    CHECK(ps.values[idxFlag][2] == 0.0f);
    CHECK(ps.values[idxFlag][3] == 0.0f);
  }

  // Invalid index should fail.
  CHECK(!ps.setValue(-1, {0, 0, 0, 0}));
  CHECK(!ps.setRawValue(99, {0, 0, 0, 0}));

  return failures ? 1 : 0;
}
