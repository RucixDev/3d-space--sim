#pragma once

#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"
#include "stellar/render/Shader.h"

#include <cstdint>
#include <string>
#include <vector>

namespace stellar::render {

// A minimal triangle descriptor used by the surfel generator, the visibility
// evaluator, and the next-best-path planner.
struct TriangleGeom {
  math::Vec3d p0{0, 0, 0};
  math::Vec3d p1{0, 0, 0};
  math::Vec3d p2{0, 0, 0};

  // Unit-length triangle normal (right-hand rule).
  math::Vec3d normal{0, 0, 1};

  // Convenience values for planning/preview.
  math::Vec3d centroid{0, 0, 0};
  double area{0.0};
};

TriangleGeom makeTriangleGeom(const math::Vec3d& p0,
                              const math::Vec3d& p1,
                              const math::Vec3d& p2);

// One geometry-aware surfel (a 2D Gaussian “disk” in 3D, aligned to the
// triangle's tangent plane).
struct GaussianSurfel {
  // Position (world)
  float px{0}, py{0}, pz{0};

  // Orthonormal basis
  float nx{0}, ny{0}, nz{1};
  float tx{1}, ty{0}, tz{0};
  float bx{0}, by{1}, bz{0};

  // World-space radii for the surfel quad. Interpreted as ~3σ extents when used
  // with the default renderer sharpness (≈ 9).
  float ru{0.05f}, rv{0.05f};

  // Debug / preview shading
  float r{1}, g{1}, b{1};
  float opacity{1.0f};

  // Optional linkage back to the source triangle / mapping.
  std::int32_t triIndex{-1};
  std::uint16_t obsCount{0};
};

// Build a geometry-aware surfel from a triangle (principal-axis ellipse on the
// triangle plane). `requiredObs` is used to map obsCount -> opacity.
GaussianSurfel makeGaussianSurfelFromTriangle(const TriangleGeom& tri,
                                              std::uint16_t obsCount,
                                              std::uint16_t requiredObs,
                                              float alphaMul = 0.85f);

// A simple camera basis used by the visibility evaluator.
struct Viewpoint {
  math::Vec3d pos{0, 0, 0};
  math::Vec3d forward{0, 0, -1};
  math::Vec3d right{1, 0, 0};
  math::Vec3d up{0, 1, 0};
};

Viewpoint makeLookAtView(const math::Vec3d& eye,
                         const math::Vec3d& target = math::Vec3d{0, 0, 0},
                         const math::Vec3d& worldUp = math::Vec3d{0, 1, 0});

// Visibility evaluation settings.
//
// If occlusionAware=true, we rasterize the mesh into a tiny Z-buffer from the
// candidate view (fast enough for interactive NBV emphasizers) and return the
// set of triangles that win at least one pixel.
struct VisibilitySettings {
  double fovYRad{math::degToRad(60.0)};
  int rasterRes{96};
  bool occlusionAware{true};
  double nearPlane{0.02};
  double farPlane{50.0};
};

struct ViewVisibility {
  std::vector<int> visibleTriangles{};
};

ViewVisibility computeViewVisibility(const std::vector<TriangleGeom>& tris,
                                     const Viewpoint& view,
                                     const VisibilitySettings& s);

// Beam-search “next-best-path” (multi-step lookahead NBV) planner settings.
struct PathPlannerSettings {
  int horizon{3};
  int beamWidth{8};

  // Penalty per radian of angular movement between successive views.
  double moveCostWeight{0.15};

  // Extra penalty when reusing a view already in the path or in the provided
  // history (captured views).
  double revisitPenalty{0.35};

  // Per-state branching factor (keeps CPU cost bounded). If <=0, defaults to
  // beamWidth.
  int localBranching{8};
};

struct PlannedViewPath {
  std::vector<int> views{};
  double score{0.0};
};

PlannedViewPath planNextBestPathBeam(const std::vector<ViewVisibility>& vis,
                                     const std::vector<math::Vec3d>& viewPositions,
                                     const std::vector<TriangleGeom>& tris,
                                     const std::vector<std::uint8_t>& obsCount,
                                     std::uint8_t requiredObs,
                                     int startViewIdx,
                                     const PathPlannerSettings& ps,
                                     const std::vector<int>* history = nullptr);

// -----------------------------------------------------------------------------
// Simple instanced Gaussian-surfels renderer (OpenGL 3.3)
// -----------------------------------------------------------------------------
class GaussianSurfelRenderer {
public:
  ~GaussianSurfelRenderer();

  bool init(std::string* outError = nullptr);

  void setViewProj(const float* view, const float* proj);

  // Controls the gaussian falloff. If you interpret ru/rv as ~3σ extents, use
  // ~9 for a very soft edge (exp(-4.5) at the quad boundary).
  void setSharpness(float v) { sharpness_ = v; }

  void draw(const std::vector<GaussianSurfel>& surfels);

private:
  ShaderProgram shader_{};

  std::uint32_t vao_{0};
  std::uint32_t vbo_{0};
  std::uint32_t instanceVbo_{0};

  float view_[16]{};
  float proj_[16]{};

  float sharpness_{9.0f};
};

} // namespace stellar::render
