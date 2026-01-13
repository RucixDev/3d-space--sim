#pragma once

#include "stellar/render/SdfMesher.h" // ScalarField3D, SdfMeshData

namespace stellar::render {

// Dual Contouring is an isosurface extraction method that places **one vertex per cell**
// (instead of per-edge like Marching Cubes/Tetrahedra).
//
// This implementation is designed for procedural geometry workflows:
//  - Deterministic sampling on a uniform grid
//  - Per-cell vertex position solved from Hermite samples using a tiny 3x3 QEF
//  - Watertight connectivity by emitting quads around grid edges with sign changes
//
// Convention: by default, inside is f <= iso and outside is f > iso.
// For signed distance fields (SDFs), iso = 0 and f < 0 is "inside".

struct SdfDualContourParams {
  // Number of cells per axis in the cubic domain.
  // Total sampled grid points is (resolution+1)^3.
  int resolution{32};

  // Half-extent of the cubic domain: positions are sampled in [-bounds, +bounds].
  float bounds{1.25f};

  // Iso value threshold.
  float iso{0.0f};

  // Compute per-vertex normals from the scalar field gradient (central differences).
  bool computeNormalsFromField{true};

  // Step size for gradient sampling (in world units).
  float normalEps{0.0025f};

  // Flip triangle winding so generated geometry is front-facing for standard backface
  // culling when the field is an SDF (positive outside).
  bool fixWindingFromNormals{true};

  // ---- Dual Contouring specifics ----

  // QEF regularization (added to the diagonal of A^T A).
  // Larger values bias vertices toward the average intersection point and help
  // with near-singular cases.
  float qefRegularization{1.0e-6f};

  // Clamp the solved vertex position into the cell bounds. Helps prevent spikes.
  bool clampToCell{true};

  // After solving the QEF, optionally project the vertex onto the iso-surface by
  // taking a few Newton steps along the field gradient.
  bool projectToIso{true};

  // Number of projection iterations when projectToIso is enabled.
  int projectIterations{2};
};

SdfMeshData meshIsosurfaceDualContouring(const ScalarField3D& field,
                                        const SdfDualContourParams& params = {});

} // namespace stellar::render
