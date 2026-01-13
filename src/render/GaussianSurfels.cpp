#include "stellar/render/GaussianSurfels.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

namespace stellar::render {

TriangleGeom makeTriangleGeom(const math::Vec3d& p0,
                              const math::Vec3d& p1,
                              const math::Vec3d& p2) {
  TriangleGeom t{};
  t.p0 = p0;
  t.p1 = p1;
  t.p2 = p2;

  const math::Vec3d e0 = p1 - p0;
  const math::Vec3d e1 = p2 - p0;
  const math::Vec3d cr = math::cross(e0, e1);
  const double len = cr.length();
  t.area = 0.5 * len;
  t.centroid = (p0 + p1 + p2) / 3.0;

  if (len > 1e-12) {
    t.normal = cr / len;
  } else {
    t.normal = math::Vec3d{0, 0, 1};
  }
  return t;
}

static math::Vec3d safeNormalize(const math::Vec3d& v, const math::Vec3d& fallback) {
  const double lsq = v.lengthSq();
  if (lsq < 1e-18) return fallback;
  return v / std::sqrt(lsq);
}

GaussianSurfel makeGaussianSurfelFromTriangle(const TriangleGeom& tri,
                                              std::uint16_t obsCount,
                                              std::uint16_t requiredObs,
                                              float alphaMul) {
  GaussianSurfel s{};
  s.px = (float)tri.centroid.x;
  s.py = (float)tri.centroid.y;
  s.pz = (float)tri.centroid.z;

  math::Vec3d n = tri.normal;
  n = safeNormalize(n, math::Vec3d{0, 0, 1});

  // Start from an arbitrary but stable tangent (edge-based).
  math::Vec3d t0 = tri.p1 - tri.p0;
  t0 = safeNormalize(t0, math::Vec3d{1, 0, 0});
  math::Vec3d b0 = math::cross(n, t0);
  b0 = safeNormalize(b0, math::Vec3d{0, 1, 0});
  // Re-orthonormalize.
  t0 = safeNormalize(math::cross(b0, n), math::Vec3d{1, 0, 0});
  b0 = safeNormalize(math::cross(n, t0), math::Vec3d{0, 1, 0});

  // Project triangle vertices into the local tangent frame about the centroid.
  const math::Vec3d c = tri.centroid;
  const math::Vec3d p[3] = {tri.p0, tri.p1, tri.p2};

  double u[3]{}, v[3]{};
  for (int i = 0; i < 3; ++i) {
    const math::Vec3d d = p[i] - c;
    u[i] = math::dot(d, t0);
    v[i] = math::dot(d, b0);
  }

  // 2x2 covariance in (t0,b0) basis. Mean is ~0 due to centroid.
  const double invN = 1.0 / 3.0;
  const double covXX = invN * (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  const double covYY = invN * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  const double covXY = invN * (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);

  // Eigen-decomposition of symmetric 2x2.
  const double trace = covXX + covYY;
  const double det = covXX * covYY - covXY * covXY;
  const double halfTrace = 0.5 * trace;
  const double disc = std::sqrt(std::max(halfTrace * halfTrace - det, 0.0));
  const double lambda1 = halfTrace + disc;
  const double lambda2 = halfTrace - disc;

  // Principal eigenvector in uv-plane.
  double vx = 1.0;
  double vy = 0.0;
  if (std::abs(covXY) > 1e-14 || std::abs(lambda1 - covXX) > 1e-14) {
    vx = covXY;
    vy = (lambda1 - covXX);
  }
  const double vlen = std::sqrt(vx*vx + vy*vy);
  if (vlen > 1e-18) { vx /= vlen; vy /= vlen; }
  else { vx = 1.0; vy = 0.0; }

  // Rotate the frame into principal-axis space.
  math::Vec3d t = safeNormalize(t0 * vx + b0 * vy, t0);
  math::Vec3d b = safeNormalize(math::cross(n, t), b0);
  t = safeNormalize(math::cross(b, n), t0);

  // Project the vertices into the eigen frame to estimate extents.
  double maxAbsU = 0.0;
  double maxAbsV = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double ue = u[i] * vx + v[i] * vy;
    const double ve = -u[i] * vy + v[i] * vx;
    maxAbsU = std::max(maxAbsU, std::abs(ue));
    maxAbsV = std::max(maxAbsV, std::abs(ve));
  }

  // Convert eigenvalues (variances) to a conservative ~3σ extent, but also
  // ensure the quad encloses the triangle's projected vertices.
  constexpr double kEps = 1e-12;
  const double sigmaU = std::sqrt(std::max(lambda1, kEps));
  const double sigmaV = std::sqrt(std::max(lambda2, kEps));

  double ru = std::max(3.0 * sigmaU, maxAbsU * 1.05);
  double rv = std::max(3.0 * sigmaV, maxAbsV * 1.05);

  // Guard against degenerate / extremely tiny triangles (keeps them visible).
  const double minR = 0.0025;
  ru = std::max(ru, minR);
  rv = std::max(rv, minR);

  s.nx = (float)n.x; s.ny = (float)n.y; s.nz = (float)n.z;
  s.tx = (float)t.x; s.ty = (float)t.y; s.tz = (float)t.z;
  s.bx = (float)b.x; s.by = (float)b.y; s.bz = (float)b.z;
  s.ru = (float)ru;
  s.rv = (float)rv;

  // Debug color from normal (handy for checking reconstruction coverage).
  s.r = (float)math::clamp(0.5 * (n.x + 1.0), 0.0, 1.0);
  s.g = (float)math::clamp(0.5 * (n.y + 1.0), 0.0, 1.0);
  s.b = (float)math::clamp(0.5 * (n.z + 1.0), 0.0, 1.0);

  if (requiredObs == 0) requiredObs = 1;
  const float conf = (float)math::clamp((double)obsCount / (double)requiredObs, 0.0, 1.0);
  s.opacity = math::clamp(alphaMul * conf, 0.0f, 1.0f);

  s.triIndex = -1;
  s.obsCount = obsCount;

  return s;
}

Viewpoint makeLookAtView(const math::Vec3d& eye,
                         const math::Vec3d& target,
                         const math::Vec3d& worldUp) {
  Viewpoint v{};
  v.pos = eye;

  v.forward = safeNormalize(target - eye, math::Vec3d{0, 0, -1});
  v.right = safeNormalize(math::cross(v.forward, worldUp), math::Vec3d{1, 0, 0});
  v.up = safeNormalize(math::cross(v.right, v.forward), math::Vec3d{0, 1, 0});
  return v;
}

static bool projectToNdc(const Viewpoint& view,
                         double tanHalfFovY,
                         double nearPlane,
                         const math::Vec3d& p,
                         double& outX,
                         double& outY,
                         double& outZ) {
  const math::Vec3d d = p - view.pos;
  const double x = math::dot(d, view.right);
  const double y = math::dot(d, view.up);
  const double z = math::dot(d, view.forward);

  if (z <= nearPlane) return false;

  outX = x / (z * tanHalfFovY);
  outY = y / (z * tanHalfFovY);
  outZ = z;
  return true;
}

ViewVisibility computeViewVisibility(const std::vector<TriangleGeom>& tris,
                                     const Viewpoint& view,
                                     const VisibilitySettings& s) {
  ViewVisibility out{};

  if (tris.empty()) return out;

  const double tanHalfFovY = std::tan(std::max(1e-6, s.fovYRad * 0.5));
  const double nearPlane = std::max(1e-6, s.nearPlane);

  // Fast path: no occlusion, just a front-face + frustum check.
  if (!s.occlusionAware) {
    out.visibleTriangles.reserve(tris.size() / 2);

    for (int i = 0; i < (int)tris.size(); ++i) {
      const TriangleGeom& t = tris[(std::size_t)i];

      // Front-face test.
      const math::Vec3d toCam = safeNormalize(view.pos - t.centroid, math::Vec3d{0, 0, 1});
      const double cosAng = math::dot(t.normal, toCam);
      if (cosAng <= 0.0) continue;

      // Cheap frustum test: project the centroid.
      double xn, yn, zn;
      if (!projectToNdc(view, tanHalfFovY, nearPlane, t.centroid, xn, yn, zn)) continue;
      if (xn < -1.1 || xn > 1.1 || yn < -1.1 || yn > 1.1) continue;

      out.visibleTriangles.push_back(i);
    }
    return out;
  }

  // Occlusion-aware path: tiny CPU rasterizer into a depth buffer.
  const int res = std::max(8, s.rasterRes);
  const int w = res;
  const int h = res;
  std::vector<float> depth((std::size_t)w * (std::size_t)h, 1e30f);
  std::vector<int> triId((std::size_t)w * (std::size_t)h, -1);

  auto toScreen = [&](double xNdc, double yNdc, double& sx, double& sy) {
    // xNdc,yNdc in [-1,1] => pixel center space.
    sx = (xNdc * 0.5 + 0.5) * (double)w;
    sy = (0.5 - yNdc * 0.5) * (double)h; // flip Y so +Y is up.
  };

  for (int ti = 0; ti < (int)tris.size(); ++ti) {
    const TriangleGeom& t = tris[(std::size_t)ti];

    // Skip back-faces early.
    const math::Vec3d toCam = safeNormalize(view.pos - t.centroid, math::Vec3d{0, 0, 1});
    const double cosAng = math::dot(t.normal, toCam);
    if (cosAng <= 0.0) continue;

    double x0, y0, z0;
    double x1, y1, z1;
    double x2, y2, z2;
    if (!projectToNdc(view, tanHalfFovY, nearPlane, t.p0, x0, y0, z0)) continue;
    if (!projectToNdc(view, tanHalfFovY, nearPlane, t.p1, x1, y1, z1)) continue;
    if (!projectToNdc(view, tanHalfFovY, nearPlane, t.p2, x2, y2, z2)) continue;

    // Clip to a slightly expanded NDC square.
    if ((x0 < -2.0 && x1 < -2.0 && x2 < -2.0) ||
        (x0 >  2.0 && x1 >  2.0 && x2 >  2.0) ||
        (y0 < -2.0 && y1 < -2.0 && y2 < -2.0) ||
        (y0 >  2.0 && y1 >  2.0 && y2 >  2.0)) {
      continue;
    }

    double sx0, sy0;
    double sx1, sy1;
    double sx2, sy2;
    toScreen(x0, y0, sx0, sy0);
    toScreen(x1, y1, sx1, sy1);
    toScreen(x2, y2, sx2, sy2);

    const double minX = std::floor(std::min({sx0, sx1, sx2}));
    const double maxX = std::ceil (std::max({sx0, sx1, sx2}));
    const double minY = std::floor(std::min({sy0, sy1, sy2}));
    const double maxY = std::ceil (std::max({sy0, sy1, sy2}));

    int xMin = (int)std::max(0.0, minX);
    int yMin = (int)std::max(0.0, minY);
    int xMax = (int)std::min((double)(w - 1), maxX);
    int yMax = (int)std::min((double)(h - 1), maxY);
    if (xMin > xMax || yMin > yMax) continue;

    const double denom = (sy1 - sy2) * (sx0 - sx2) + (sx2 - sx1) * (sy0 - sy2);
    if (std::abs(denom) < 1e-12) continue;

    // Rasterize.
    for (int py = yMin; py <= yMax; ++py) {
      for (int px = xMin; px <= xMax; ++px) {
        const double fx = (double)px + 0.5;
        const double fy = (double)py + 0.5;

        const double w0 = ((sy1 - sy2) * (fx - sx2) + (sx2 - sx1) * (fy - sy2)) / denom;
        const double w1 = ((sy2 - sy0) * (fx - sx2) + (sx0 - sx2) * (fy - sy2)) / denom;
        const double w2 = 1.0 - w0 - w1;

        if (w0 < -1e-6 || w1 < -1e-6 || w2 < -1e-6) continue;

        const double z = w0 * z0 + w1 * z1 + w2 * z2;
        const std::size_t idx = (std::size_t)py * (std::size_t)w + (std::size_t)px;
        if (z < (double)depth[idx]) {
          depth[idx] = (float)z;
          triId[idx] = ti;
        }
      }
    }
  }

  std::vector<std::uint8_t> seen(tris.size(), 0u);
  for (int id : triId) {
    if (id >= 0) seen[(std::size_t)id] = 1u;
  }

  out.visibleTriangles.reserve(tris.size() / 3);
  for (int i = 0; i < (int)seen.size(); ++i) {
    if (seen[(std::size_t)i]) out.visibleTriangles.push_back(i);
  }

  return out;
}

// ----------------------------------------------------------------------------
// Next-best-path planner (beam search with multi-step lookahead)
// ----------------------------------------------------------------------------

static bool containsInt(const std::vector<int>& v, int x) {
  for (int a : v) if (a == x) return true;
  return false;
}

static double angularDistanceRad(const math::Vec3d& a, const math::Vec3d& b) {
  const math::Vec3d na = safeNormalize(a, math::Vec3d{1, 0, 0});
  const math::Vec3d nb = safeNormalize(b, math::Vec3d{1, 0, 0});
  const double d = math::clamp(math::dot(na, nb), -1.0, 1.0);
  return std::acos(d);
}

PlannedViewPath planNextBestPathBeam(const std::vector<ViewVisibility>& vis,
                                     const std::vector<math::Vec3d>& viewPositions,
                                     const std::vector<TriangleGeom>& tris,
                                     const std::vector<std::uint8_t>& obsCount,
                                     std::uint8_t requiredObs,
                                     int startViewIdx,
                                     const PathPlannerSettings& ps,
                                     const std::vector<int>* history) {
  PlannedViewPath out{};

  const int viewCount = (int)vis.size();
  if (viewCount == 0) return out;
  if ((int)viewPositions.size() != viewCount) return out;
  if (tris.empty() || obsCount.size() != tris.size()) return out;

  const int horizon = std::max(0, ps.horizon);
  if (horizon <= 0) return out;

  const int beamWidth = std::max(1, ps.beamWidth);
  int branching = ps.localBranching;
  if (branching <= 0) branching = beamWidth;
  branching = std::max(1, std::min(branching, viewCount));

  if (requiredObs == 0) requiredObs = 1;

  auto gainForView = [&](int viewIdx, const std::vector<std::uint8_t>& obs) -> double {
    const auto& vv = vis[(std::size_t)viewIdx].visibleTriangles;
    double gain = 0.0;
    const math::Vec3d eye = viewPositions[(std::size_t)viewIdx];

    for (int triIdx : vv) {
      if (triIdx < 0 || triIdx >= (int)tris.size()) continue;

      const std::uint8_t c = obs[(std::size_t)triIdx];
      if (c >= requiredObs) continue;

      const std::uint8_t remaining = (std::uint8_t)(requiredObs - c);
      const TriangleGeom& t = tris[(std::size_t)triIdx];

      const math::Vec3d toCam = safeNormalize(eye - t.centroid, math::Vec3d{0, 0, 1});
      const double cosAng = math::dot(t.normal, toCam);
      if (cosAng <= 0.0) continue;

      const double w = (double)remaining / (double)requiredObs;
      gain += t.area * w * cosAng;
    }
    return gain;
  };

  struct State {
    std::vector<int> path{};
    int last{-1};
    double score{0.0};
    std::vector<std::uint8_t> obs{};
  };

  State root{};
  root.last = startViewIdx;
  root.score = 0.0;
  root.obs = obsCount;

  std::vector<State> beam{};
  beam.push_back(std::move(root));

  for (int depth = 0; depth < horizon; ++depth) {
    std::vector<State> next{};
    next.reserve((std::size_t)beam.size() * (std::size_t)branching);

    for (const State& st : beam) {
      // Find local top-k candidates for this state.
      std::vector<std::pair<double, int>> best{};
      best.reserve((std::size_t)branching);

      for (int v = 0; v < viewCount; ++v) {
        if (v == st.last) continue;

        double stepScore = gainForView(v, st.obs);

        // Motion penalty.
        if (st.last >= 0 && st.last < viewCount) {
          const double ang = angularDistanceRad(viewPositions[(std::size_t)st.last],
                                                viewPositions[(std::size_t)v]);
          stepScore -= ps.moveCostWeight * ang;
        }

        // Revisit penalty.
        bool revisit = containsInt(st.path, v);
        if (!revisit && history) revisit = containsInt(*history, v);
        if (revisit) stepScore -= ps.revisitPenalty;

        // Keep only top `branching` scores.
        if ((int)best.size() < branching) {
          best.emplace_back(stepScore, v);
          std::push_heap(best.begin(), best.end(),
                         [](const auto& a, const auto& b) { return a.first > b.first; });
        } else {
          // best is a min-heap by score at the front when using inverted compare.
          auto worst = best.front();
          if (stepScore > worst.first) {
            std::pop_heap(best.begin(), best.end(),
                          [](const auto& a, const auto& b) { return a.first > b.first; });
            best.back() = {stepScore, v};
            std::push_heap(best.begin(), best.end(),
                           [](const auto& a, const auto& b) { return a.first > b.first; });
          }
        }
      }

      // Heap -> sorted descending.
      std::sort(best.begin(), best.end(),
                [](const auto& a, const auto& b) { return a.first > b.first; });

      // Expand.
      for (const auto& kv : best) {
        const double stepScore = kv.first;
        const int vIdx = kv.second;

        State child{};
        child.path = st.path;
        child.path.push_back(vIdx);
        child.last = vIdx;
        child.score = st.score + stepScore;
        child.obs = st.obs;

        // Apply observation update.
        const auto& vv = vis[(std::size_t)vIdx].visibleTriangles;
        for (int triIdx : vv) {
          if (triIdx < 0 || triIdx >= (int)child.obs.size()) continue;
          if (child.obs[(std::size_t)triIdx] < requiredObs) child.obs[(std::size_t)triIdx] += 1;
        }

        next.push_back(std::move(child));
      }
    }

    if (next.empty()) break;

    // Keep global top beamWidth states.
    std::sort(next.begin(), next.end(),
              [](const State& a, const State& b) { return a.score > b.score; });
    if ((int)next.size() > beamWidth) next.resize((std::size_t)beamWidth);

    beam = std::move(next);
  }

  if (beam.empty()) return out;
  out.views = beam.front().path;
  out.score = beam.front().score;
  return out;
}

// ----------------------------------------------------------------------------
// OpenGL Gaussian surfel renderer
// ----------------------------------------------------------------------------

GaussianSurfelRenderer::~GaussianSurfelRenderer() {
  if (instanceVbo_) {
    gl::DeleteBuffers(1, &instanceVbo_);
    instanceVbo_ = 0;
  }
  if (vbo_) {
    gl::DeleteBuffers(1, &vbo_);
    vbo_ = 0;
  }
  if (vao_) {
    gl::DeleteVertexArrays(1, &vao_);
    vao_ = 0;
  }
}

static const char* kSurfelsVS = R"GLSL(
#version 330 core
layout(location=0) in vec2 aCorner;

layout(location=1) in vec3 iPos;
layout(location=2) in vec3 iTangent;
layout(location=3) in vec3 iBitangent;
layout(location=4) in vec2 iRadii;
layout(location=5) in vec4 iColorOpacity;

uniform mat4 uView;
uniform mat4 uProj;

out vec2 vCorner;
out vec4 vColorOpacity;

void main() {
  vec3 world = iPos
             + aCorner.x * iRadii.x * iTangent
             + aCorner.y * iRadii.y * iBitangent;

  gl_Position = uProj * uView * vec4(world, 1.0);
  vCorner = aCorner;
  vColorOpacity = iColorOpacity;
}
)GLSL";

static const char* kSurfelsFS = R"GLSL(
#version 330 core
in vec2 vCorner;
in vec4 vColorOpacity;

uniform float uSharpness;

out vec4 FragColor;

void main() {
  float r2 = dot(vCorner, vCorner);
  float w = exp(-0.5 * uSharpness * r2);

  float a = clamp(vColorOpacity.a * w, 0.0, 1.0);

  // Premultiplied alpha output (recommended blend: ONE, ONE_MINUS_SRC_ALPHA)
  FragColor = vec4(vColorOpacity.rgb * a, a);
}
)GLSL";

bool GaussianSurfelRenderer::init(std::string* outError) {
  if (!shader_.build(kSurfelsVS, kSurfelsFS, outError)) return false;

  gl::GenVertexArrays(1, &vao_);
  gl::GenBuffers(1, &vbo_);
  gl::GenBuffers(1, &instanceVbo_);

  // 2-triangle quad in clip-local space.
  const float quad[12] = {
    -1.0f, -1.0f,
     1.0f, -1.0f,
     1.0f,  1.0f,
    -1.0f, -1.0f,
     1.0f,  1.0f,
    -1.0f,  1.0f
  };

  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, vbo_);
  gl::BufferData(GL_ARRAY_BUFFER, (GLsizeiptr)sizeof(quad), quad, GL_STATIC_DRAW);

  gl::EnableVertexAttribArray(0);
  gl::VertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 2, (void*)0);

  // Instance attributes will be configured in draw() after uploading instance buffer.
  gl::BindBuffer(GL_ARRAY_BUFFER, instanceVbo_);
  gl::BufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);

  struct Instance {
    float px, py, pz;
    float tx, ty, tz;
    float bx, by, bz;
    float ru, rv;
    float r, g, b, a;
  };

  const GLsizei stride = (GLsizei)sizeof(Instance);

  gl::EnableVertexAttribArray(1);
  gl::VertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Instance, px));
  gl::VertexAttribDivisor(1, 1);

  gl::EnableVertexAttribArray(2);
  gl::VertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Instance, tx));
  gl::VertexAttribDivisor(2, 1);

  gl::EnableVertexAttribArray(3);
  gl::VertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Instance, bx));
  gl::VertexAttribDivisor(3, 1);

  gl::EnableVertexAttribArray(4);
  gl::VertexAttribPointer(4, 2, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Instance, ru));
  gl::VertexAttribDivisor(4, 1);

  gl::EnableVertexAttribArray(5);
  gl::VertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, stride, (void*)offsetof(Instance, r));
  gl::VertexAttribDivisor(5, 1);

  // Default to identity matrices.
  for (int i = 0; i < 16; ++i) {
    view_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
    proj_[i] = (i % 5 == 0) ? 1.0f : 0.0f;
  }

  return true;
}

void GaussianSurfelRenderer::setViewProj(const float* view, const float* proj) {
  std::memcpy(view_, view, sizeof(float) * 16);
  std::memcpy(proj_, proj, sizeof(float) * 16);
}

void GaussianSurfelRenderer::draw(const std::vector<GaussianSurfel>& surfels) {
  if (surfels.empty()) return;

  struct Instance {
    float px, py, pz;
    float tx, ty, tz;
    float bx, by, bz;
    float ru, rv;
    float r, g, b, a;
  };

  std::vector<Instance> instances;
  instances.resize(surfels.size());

  for (std::size_t i = 0; i < surfels.size(); ++i) {
    const GaussianSurfel& s = surfels[i];
    Instance inst{};
    inst.px = s.px; inst.py = s.py; inst.pz = s.pz;
    inst.tx = s.tx; inst.ty = s.ty; inst.tz = s.tz;
    inst.bx = s.bx; inst.by = s.by; inst.bz = s.bz;
    inst.ru = s.ru; inst.rv = s.rv;
    inst.r = s.r; inst.g = s.g; inst.b = s.b;
    inst.a = s.opacity;
    instances[i] = inst;
  }

  shader_.bind();
  shader_.setUniformMat4("uView", view_);
  shader_.setUniformMat4("uProj", proj_);
  shader_.setUniform1f("uSharpness", sharpness_);

  gl::BindVertexArray(vao_);
  gl::BindBuffer(GL_ARRAY_BUFFER, instanceVbo_);
  gl::BufferData(GL_ARRAY_BUFFER,
                 (GLsizeiptr)(instances.size() * sizeof(Instance)),
                 instances.data(),
                 GL_DYNAMIC_DRAW);

  gl::DrawArraysInstanced(GL_TRIANGLES, 0, 6, (GLsizei)instances.size());
}

} // namespace stellar::render
