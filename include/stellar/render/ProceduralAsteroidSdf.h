#pragma once

#include "stellar/core/Types.h"
#include "stellar/render/ProceduralAsteroid.h" // AsteroidMeshStats
#include "stellar/render/SdfMesher.h"

namespace stellar::render {

// Higher-fidelity asteroid generator based on SDF isosurface meshing.
//
// This complements the classic displaced-UV-sphere asteroid generator by:
//  - carving true craters (boolean subtraction)
//  - supporting random fracture planes (faceting)
//  - producing watertight, non-spherical topology
//
// It is intended for tooling/preview and for small-ish on-demand variant
// libraries (as used by the VFX Lab). The result is CPU mesh data that can
// be uploaded with Mesh::upload().

struct AsteroidSdfParams {
  // Meshing quality.
  int resolution{56};
  float bounds{1.55f};
  float iso{0.0f};

  // Base shape.
  float baseRadius{0.95f};
  // Axis scale applied to the domain before evaluating the SDF.
  // Values != 1.0 produce elongated/flattened asteroids.
  float axisScaleX{1.0f};
  float axisScaleY{1.0f};
  float axisScaleZ{1.0f};

  // Two-stage fBm displacement layered on top of the boolean cratered base.
  float noise1Frequency{2.6f};
  float noise1Amplitude{0.22f};
  int noise1Octaves{5};
  float noise1Lacunarity{2.0f};
  float noise1Gain{0.52f};

  float noise2Frequency{8.0f};
  float noise2Amplitude{0.06f};
  int noise2Octaves{3};
  float noise2Lacunarity{2.2f};
  float noise2Gain{0.55f};

  // Craters are carved by subtracting spheres from a near-spherical base.
  int craterCount{18};
  float craterRadiusMinDeg{4.0f};
  float craterRadiusMaxDeg{30.0f};
  float craterDepth{0.12f};       // fraction of crater radius (0..1)
  float craterSmoothK{0.08f};     // smooth-subtract edge softness

  // Optional fracture planes to create faceting and sharper silhouettes.
  int cutCount{2};
  float cutOffsetMin{0.22f};
  float cutOffsetMax{0.60f};

  // Optional subtle grooves (trig pattern) to break up large flat areas.
  float grooveStrength{0.02f};
  float grooveFrequency{7.0f};

  // Normal sampling for shading.
  float normalEps{0.0040f};
};

// Generate a deterministic SDF-based asteroid mesh.
SdfMeshData generateAsteroidSdfMesh(core::u64 seed, const AsteroidSdfParams& params = {});

// Convenience: compute stats in the same format as the classic generator so
// higher-level systems can treat both paths uniformly.
AsteroidMeshStats measureAsteroidSdfMesh(const SdfMeshData& mesh);

} // namespace stellar::render
