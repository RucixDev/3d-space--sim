#pragma once

#include "stellar/render/SdfMesher.h" // SdfMeshData, VertexPNUT

#include <cstddef>
#include <string>

namespace stellar::render {

// CPU-side triangle mesh simplification using Garland/Heckbert style
// Quadric Error Metrics (QEM) edge collapses.
//
// The project already has multiple ways to *generate* meshes procedurally
// (SDF meshing, procedural asteroid spheres, etc.), but there was no way to
// cheaply reduce polygon count for export or for authoring LOD chains.
//
// This module intentionally stays CPU-only and works on the existing
// SdfMeshData representation (VertexPNUT + index buffer). The output mesh is
// a clean reindexed triangle list with optional recomputed normals/UVs.

struct MeshSimplifyParams {
  // Target fraction of triangles to keep, in (0,1].
  // If >= 1.0 the input mesh is returned unchanged.
  float targetTriangleRatio{1.0f};

  // Optional explicit target triangle count.
  // If > 0 this overrides targetTriangleRatio.
  int targetTriangleCount{-1};

  // Safety cap: maximum number of successful edge collapses attempted.
  int maxIterations{250000};

  // Optional error threshold: if > 0, edge collapses with cost above this
  // value are rejected and simplification stops early.
  double maxError{0.0};

  // When true, vertex normals are recomputed from triangle geometry.
  bool recomputeNormals{true};

  // When true, UVs are recomputed from the vertex position using the same
  // spherical mapping convention used by the SDF meshers.
  bool recomputeSphericalUVs{true};
};

struct MeshSimplifyStats {
  std::size_t inputVertices{0};
  std::size_t inputTriangles{0};
  std::size_t outputVertices{0};
  std::size_t outputTriangles{0};

  // Time spent inside simplifyMeshQEM (milliseconds).
  double ms{0.0};

  // Largest collapse cost actually accepted.
  double maxAcceptedError{0.0};
};

// Simplify a triangle mesh with QEM edge collapses.
//
// On failure, returns the input mesh unchanged and writes a human-readable
// error string to outError (if provided).
SdfMeshData simplifyMeshQEM(const SdfMeshData& mesh,
                           const MeshSimplifyParams& params = {},
                           MeshSimplifyStats* outStats = nullptr,
                           std::string* outError = nullptr);

} // namespace stellar::render
