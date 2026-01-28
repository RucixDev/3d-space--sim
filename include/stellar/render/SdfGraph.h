#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/SdfMesher.h" // ScalarField3D

#include <optional>
#include <string>
#include <vector>

namespace stellar::render {

// A tiny, CPU-side procedural geometry “engine”:
// a node graph that evaluates into a signed distance function (SDF).
//
// The output can be fed directly into SdfMesher (Marching Tetrahedra) to
// generate triangle meshes.
//
// NOTE: Unlike the GPU ProceduralGraph, this graph is evaluated on the CPU and
// is intended to be used for geometry (meshing), not 2D texture baking.

constexpr int kSdfGraphMaxNodes = 64;

enum class SdfNodeOp : core::u8 {
  Constant = 0,

  // Primitives (SDF, iso=0).
  Sphere,
  Box,
  Capsule,
  TorusY,
  Plane,

  // CSG / blending.
  Union,
  SmoothUnion,
  Intersect,
  Subtract,

  // Modifiers.
  NoiseDisplace,
  NoiseDisplacePerlin,
  Shell,

  // Space transforms (domain operations).
  Translate,
  RotateX,
  RotateY,
  RotateZ,
  Scale,
  Repeat,
  Mirror,
  TwistY,
};

const char* sdfNodeOpName(SdfNodeOp op);

// Canonical, file-format-friendly id for the op (no spaces or punctuation).
// Intended for graph serialization.
const char* sdfNodeOpId(SdfNodeOp op);

// Parses either sdfNodeOpId(), sdfNodeOpName(), or common aliases (case-insensitive).
// Returns std::nullopt on failure.
std::optional<SdfNodeOp> sdfNodeOpFromString(const std::string& s);

struct SdfNode {
  SdfNodeOp op{SdfNodeOp::Sphere};

  // Input node indices (use -1 for none).
  // Canonical names: in0 / in1. Legacy aliases: inA / inB.
  union {
    int in0{-1};
    int inA;
  };
  union {
    int in1{-1};
    int inB;
  };

  // Freeform parameters; interpretation depends on op.
  // Canonical names: p0..p7. Legacy alias: p[8].
  // Keep this POD-ish so it copies nicely for async jobs.
  union {
    struct {
      float p0{0.0f};
      float p1{0.0f};
      float p2{0.0f};
      float p3{0.0f};
      float p4{0.0f};
      float p5{0.0f};
      float p6{0.0f};
      float p7{0.0f};
    };
    float p[8];
  };

  // Optional node-local seed tweak (used by noise modifiers).
  core::u64 seed{0};
};


struct SdfGraph {
  core::u64 seed{0xBADC0FFEEULL};

  std::vector<SdfNode> nodes;
  int output{-1};

  static SdfGraph makeDefault();
};

enum class SdfGraphPreset : core::u8 {
  Asteroid = 0,
  Crystal,
  Torus,
  HollowBox,
  BooleanDemo,
};

const char* sdfGraphPresetName(SdfGraphPreset preset);
SdfGraph makeSdfGraphPreset(SdfGraphPreset preset, core::u64 seed);

// ---- Graph file I/O (human-readable, versioned text) ----
bool saveSdfGraphToFile(const SdfGraph& g, const std::string& path, std::string* outError = nullptr);
bool loadSdfGraphFromFile(const std::string& path, SdfGraph& out, std::string* outError = nullptr);

// Evaluate signed distance at a world-space point.
float evalSdfGraph(const SdfGraph& g, float x, float y, float z);

// Convenience: convert graph into a ScalarField3D compatible with SdfMesher.
// The returned field captures a *copy* of the graph (safe for async meshing).
ScalarField3D makeSdfField(const SdfGraph& g);

} // namespace stellar::render
