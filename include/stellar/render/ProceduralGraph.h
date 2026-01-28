#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/RenderTarget.h"
#include "stellar/render/Shader.h"
#include "stellar/render/Texture.h"

#include <array>
#include <optional>
#include <string>
#include <vector>

namespace stellar::render {

// A tiny "procedural generation graphics engine" built around a scalar node graph
// that compiles into a GLSL fragment shader and bakes into an offscreen Texture2D.
//
// Design goals:
//  - Deterministic: (graph structure + params + seed) -> same result.
//  - Fast iteration: structure changes rebuild shader; param changes are uniforms.
//  - Self-contained: only depends on existing RenderTarget/Shader utilities.

constexpr int kProcGraphMaxNodes = 64;
constexpr int kProcGraphMaxPaletteStops = 8;

enum class ProcNodeOp : core::u8 {
  Constant = 0,

  // Coordinate sources
  UvX,
  UvY,

  // Math
  Add,
  Sub,
  Mul,
  Div,
  Min,
  Max,
  Abs,
  Invert,
  Fract,
  Clamp01,
  Smoothstep,
  Pow,
  Sine,

  // Noise / patterns
  Noise2D,
  Fbm2D,
  Perlin2D,
  FbmPerlin2D,
  RidgedFbmPerlin2D,
  Voronoi2D,

  // Domain ops (re-evaluates its input node at a modified UV)
  Warp,
  Pan,
};

const char* procNodeOpName(ProcNodeOp op);

// Parses either the canonical enum name (case-insensitive) or common aliases.
// Returns std::nullopt on failure.
std::optional<ProcNodeOp> procNodeOpFromString(const std::string& s);

struct ProcNode {
  ProcNodeOp op{ProcNodeOp::Constant};

  // Input node indices (use -1 for "none").
  // Canonical: in0/in1/in2. Legacy aliases: inA/inB/inC.
  union {
    int in0{-1};
    int inA;
  };
  union {
    int in1{-1};
    int inB;
  };
  union {
    int in2{-1};
    int inC;
  };

  // Freeform parameters; meaning depends on op.
  // Canonical: p0..p3. Legacy alias: p[4].
  union {
    struct {
      float p0{0.0f};
      float p1{0.0f};
      float p2{0.0f};
      float p3{0.0f};
    };
    float p[4];
  };

  // Optional node-local seed tweak.
  core::u64 seed{0};
};

struct ProcPaletteStop {
  // Canonical: pos, r/g/b. Legacy aliases: t, rgb[3].
  union {
    float pos{0.0f}; // [0,1]
    float t;
  };
  union {
    struct {
      float r{0.0f};
      float g{0.0f};
      float b{0.0f};
    };
    float rgb[3];
  };
};

struct ProcGraph {
  core::u64 seed{0xC0FFEEULL};

  std::vector<ProcNode> nodes;
  int output{-1};

  bool usePalette{true};
  int paletteCount{4};
  std::array<ProcPaletteStop, kProcGraphMaxPaletteStops> palette{};

  static ProcGraph makeDefault();
};

enum class ProcGraphPreset : core::u8 {
  Nebula = 0,
  Marble,
  Lava,
  AlienCircuit,
  Rocky,
};

const char* procGraphPresetName(ProcGraphPreset preset);
ProcGraph makeProceduralGraphPreset(ProcGraphPreset preset, core::u64 seed);

// ---- Graph file I/O (human-readable, versioned text) ----
//
// These helpers exist so the procedural engine can be used as an *asset pipeline*:
// users can save/share graphs, and other tools/windows can load them.
//
// File format is stable-ish but versioned so it can evolve.

bool saveProcGraphToFile(const ProcGraph& g, const std::string& path, std::string* outError = nullptr);
bool loadProcGraphFromFile(const std::string& path, ProcGraph& out, std::string* outError = nullptr);

struct ProcBakeStats {
  bool shaderRebuilt{false};
  double shaderBuildMs{0.0};
  double drawMs{0.0};

  bool mipsGenerated{false};
  double mipGenMs{0.0};
};

class ProceduralGraphBaker {
public:
  ProceduralGraphBaker() = default;
  ~ProceduralGraphBaker();

  ProceduralGraphBaker(const ProceduralGraphBaker&) = delete;
  ProceduralGraphBaker& operator=(const ProceduralGraphBaker&) = delete;

  // Quality knobs.
  //
  // - generateMips: allocates a mip chain for the baked texture, and generates it after each bake.
  //   This reduces aliasing when the baked texture is applied to 3D previews.
  // - ditherStrength: small ordered noise (in 8-bit space) applied before writing to RGBA8.
  //   Helps reduce color banding in smooth procedural gradients.
  void setGenerateMips(bool v) { generateMips_ = v; }
  bool generateMips() const { return generateMips_; }

  void setDitherStrength(float v) {
    // Keep within a sane range; values above ~2.0 start to look intentionally noisy.
    if (v < 0.0f) v = 0.0f;
    if (v > 2.0f) v = 2.0f;
    ditherStrength_ = v;
  }
  float ditherStrength() const { return ditherStrength_; }

  // When enabled, the scalar graph output `t` (clamped to [0,1]) is also written to
  // the alpha channel of the baked RGBA8 texture.
  //
  // This is intentionally optional because alpha is commonly used for transparency
  // in downstream pipelines. In tooling (raymarch preview, material prototyping)
  // the packed height is useful for deriving micro-normals.
  void setPackHeightInAlpha(bool v) { packHeightInAlpha_ = v; }
  bool packHeightInAlpha() const { return packHeightInAlpha_; }

  // Bake graph into an offscreen texture.
  // Returns false on shader compilation failure or FBO init failure; outError describes the problem.
  bool bake(const ProcGraph& g, int width, int height, float timeSec, std::string* outError = nullptr);

  const Texture2D& texture() const { return target_.color(); }
  bool isReady() const { return target_.isInited() && shader_.handle() != 0; }

  const ProcBakeStats& stats() const { return stats_; }
  const std::string& lastFragmentSource() const { return lastFragSrc_; }

private:
  static core::u64 structureKey(const ProcGraph& g);
  static std::string buildFragmentShader(const ProcGraph& g);

  bool ensureShader(const ProcGraph& g, std::string* outError);

  // Bake target (no depth needed; optional mip chain).
  RenderTarget2D target_{};

  ShaderProgram shader_{};
  core::u64 shaderKey_{0};
  unsigned int vao_{0};

  ProcBakeStats stats_{};
  std::string lastFragSrc_{};

  bool generateMips_{true};
  float ditherStrength_{1.0f};
  bool packHeightInAlpha_{false};

  // Packed uniform data (fixed max sizes).
  std::array<float, kProcGraphMaxNodes * 4> packedParams_{};
  std::array<int, kProcGraphMaxNodes> packedSeeds_{};
  std::array<float, kProcGraphMaxPaletteStops * 4> packedPalette_{};
};

} // namespace stellar::render
