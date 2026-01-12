#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/SdfMesher.h"

namespace stellar::render {

// A showcase generator built on top of the SdfMesher.
//
// The resulting mesh is intended to look like a "mysterious alien artifact"...
// It's essentially an SDF: sphere base + smooth metaball unions + symmetric noise
// displacement + planar cuts (faceting) + subtle grooves.
//
// Because it uses the generic mesher, it doubles as a reusable graphics-engine
// subsystem for future procedural shapes.

struct ArtifactParams {
  // Meshing domain and quality.
  int resolution{48};
  float bounds{1.35f};
  float iso{0.0f};

  // Base SDF (sphere).
  float baseRadius{0.90f};

  // Noise displacement (fBm).
  float noiseFrequency{2.1f};
  float noiseAmplitude{0.30f};
  int noiseOctaves{5};
  float noiseLacunarity{2.05f};
  float noiseGain{0.50f};

  // Smooth-union metaball "lumps".
  int blobCount{6};
  float blobRadiusMin{0.18f};
  float blobRadiusMax{0.55f};
  float blobCenterRadius{0.75f};
  float blobSmoothK{0.20f};

  // Random planar cuts (intersection with half-spaces) to create faceting.
  int cutCount{4};
  float cutOffsetMin{0.15f};
  float cutOffsetMax{0.55f};

  // Subtle trigonometric grooves for an "engineered" feel.
  float grooveStrength{0.10f};
  float grooveFrequency{8.0f};

  // Normal sampling for shading.
  float normalEps{0.0040f};
};

SdfMeshData generateArtifactMesh(core::u64 seed, const ArtifactParams& params = {});

} // namespace stellar::render
