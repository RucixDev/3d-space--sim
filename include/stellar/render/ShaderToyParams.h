#pragma once

#include <array>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace stellar::render {

class ShaderProgram;

// -----------------------------
// ShaderToy-style user parameters
// -----------------------------
//
// The Procedural Shader Lab supports a small, comment-driven "parameter" system
// that compiles into uniforms and exposes live sliders/controls in the UI.
//
// Syntax (inside any pass' GLSL snippet):
//
//   // @group Simulation
//   // @param float  feed      0.0367  0.0  0.08  0.0001
//   // @param float  kill      0.0649  0.0  0.08  0.0001
//   // @param vec2   offset    0.0 0.0  -1.0 -1.0   1.0  1.0   0.01
//   // @param color  tint      0.9 0.7 0.2
//   // @endgroup
//
// Notes:
//  - Parameter names must be valid GLSL identifiers.
//  - The wrapper will auto-inject `uniform` declarations for parsed parameters,
//    so you should *not* redeclare them yourself.
//  - Parameters are global to the ShaderToyGraph: a param defined in Buffer A
//    is visible in the Image pass and vice versa.
//  - Parameter values are saved in `.stoy` graph files.

enum class ShaderToyParamType : int {
  Float = 0,
  Int   = 1,
  Bool  = 2,
  Vec2  = 3,
  Vec3  = 4,
  Color3 = 5,
};

struct ShaderToyParamDef {
  ShaderToyParamType type{ShaderToyParamType::Float};

  // Uniform identifier used in GLSL.
  std::string name{};

  // UI-facing label. If empty, defaults to name.
  std::string label{};

  // Optional UI group (set via `// @group`).
  std::string group{};

  // Default, min and max are stored as up to 4 floats.
  // Scalars use x only; vec2 uses x/y; vec3 & color use x/y/z.
  std::array<float, 4> defaultValue{0, 0, 0, 0};
  std::array<float, 4> minValue{0, 0, 0, 0};
  std::array<float, 4> maxValue{1, 1, 1, 1};

  // Suggested UI drag speed.
  float step{0.01f};
};


// Transparent hashing helpers so we can look up parameters by std::string_view
// without allocating temporary std::string objects (C++20 heterogeneous lookup).
struct ShaderToyParamNameHash {
  using is_transparent = void;
  std::size_t operator()(std::string_view sv) const noexcept {
    return std::hash<std::string_view>{}(sv);
  }
};

struct ShaderToyParamNameEq {
  using is_transparent = void;
  bool operator()(std::string_view a, std::string_view b) const noexcept { return a == b; }
};

struct ShaderToyParamSet {
  std::vector<ShaderToyParamDef> defs{};
  std::vector<std::array<float, 4>> values{}; // 1:1 with defs

  // Fast lookup by name.
  std::unordered_map<std::string, int, ShaderToyParamNameHash, ShaderToyParamNameEq> indexByName{};

  bool empty() const { return defs.empty(); }
  void clear();

  void rebuildIndex();

  // Resets all parameter values to their defaults.
  void resetToDefaults();

  // Returns -1 if not present.
  int findIndex(std::string_view name) const;
  const ShaderToyParamDef* findDef(std::string_view name) const;

  // Set by name (clamped to min/max if present). Returns false if not found.
  bool setValue(std::string_view name, const std::array<float, 4>& v);

  // Set by index (faster; avoids name lookup).
  bool setValue(int index, const std::array<float, 4>& v);


  // Set by name **without** clamping to min/max. Returns false if not found.
  //
  // This is intended for internal/engine-driven uniforms that exceed any UI range
  // (e.g. world-space coordinates, radii in light-years, etc.).
  // Unused components are still zeroed according to the parameter type.
  bool setRawValue(std::string_view name, const std::array<float, 4>& v);

  // Set by index (faster; avoids name lookup).
  bool setRawValue(int index, const std::array<float, 4>& v);


  // Builds a uniform declaration block (GLSL) for these parameters.
  // This is inserted by ShaderToy's wrapper *before* `#line 1` so snippet
  // line numbers remain stable in compiler logs.
  std::string buildUniformDecls() const;

  // Applies the current parameter values as uniforms to a linked program.
  // (Safe to call even if some uniforms are optimized out.)
  void applyToShader(const ShaderProgram& shader) const;

  // Schema equality ignores UI-only data and ordering, and only checks that
  // `name -> type` matches. This is useful to detect when *recompiling all
  // passes* is required.
  bool schemaEquals(const ShaderToyParamSet& other) const;
};

// Parse parameters from a set of user snippets.
// If preserveValuesFrom is provided, values for matching (name,type) params
// are copied over (handy when recompiling after edits).
ShaderToyParamSet parseShaderToyParamsFromSources(std::span<const std::string_view> sources,
                                                  const ShaderToyParamSet* preserveValuesFrom = nullptr);

// Returns a stable GLSL type name for the given parameter type.
const char* shaderToyParamGlslType(ShaderToyParamType t);

} // namespace stellar::render
