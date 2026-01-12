#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/Mesh.h"

#include <cstdint>
#include <vector>

namespace stellar::render {

// A tiny "procedural graphics engine" building block:
// deterministic mesh generation for rocky asteroids.
//
// This generator outputs a UV-sphere whose radius is displaced by 3D fBm noise
// and crater carving (spherical caps). The output is pure CPU data and can be
// uploaded to a GPU mesh via Mesh::upload().

struct AsteroidParams {
  // Base UV-sphere topology.
  int slices{24};
  int stacks{14};

  // Radius field (unit sphere -> displaced sphere).
  float baseRadius{1.0f};
  float noiseFrequency{2.60f};
  float noiseAmplitude{0.22f};
  int noiseOctaves{5};
  float noiseLacunarity{2.0f};
  float noiseGain{0.5f};

  // Crater carving. (All values are in "unit sphere" scale.)
  int craterCount{16};
  float craterRadiusMinDeg{6.0f};
  float craterRadiusMaxDeg{22.0f};
  float craterDepth{0.10f}; // relative radius delta at crater center
  float craterRim{0.25f};   // 0..1 small ridge near the crater edge

  // Safety clamp so extreme parameter combos never invert the mesh.
  float minRadius{0.65f};
  float maxRadius{1.35f};
};

struct AsteroidMeshData {
  std::vector<VertexPNUT> vertices;
  std::vector<std::uint32_t> indices;
};

struct AsteroidMeshStats {
  std::size_t vertexCount{0};
  std::size_t indexCount{0};
  float minRadius{0.0f};
  float maxRadius{0.0f};
};

// Generate a deterministic asteroid mesh.
//
// Notes:
// - The mesh is centered at the origin.
// - The average radius is ~params.baseRadius.
// - Normals are recomputed from displaced triangles (smoother lighting).
AsteroidMeshData generateAsteroidMesh(core::u64 seed, const AsteroidParams& params = {});

AsteroidMeshStats measureAsteroidMesh(const AsteroidMeshData& mesh);

} // namespace stellar::render
