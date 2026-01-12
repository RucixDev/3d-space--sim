#include "stellar/render/SdfMesher.h"

#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::render {

namespace {

static inline math::Vec3d lerp3(const math::Vec3d& a, const math::Vec3d& b, double t) {
  return a + (b - a) * t;
}

static inline math::Vec3d safeNormalize(const math::Vec3d& v, double eps = 1e-18) {
  return v.normalized(eps);
}

static inline void sphericalUV(const math::Vec3d& p, float& outU, float& outV) {
  const double r2 = p.x * p.x + p.y * p.y + p.z * p.z;
  const double r = std::sqrt(std::max(0.0, r2));
  const double invR = 1.0 / (r + 1e-18);
  const double nx = p.x * invR;
  const double ny = p.y * invR;
  const double nz = p.z * invR;

  double u = std::atan2(nz, nx) / (2.0 * math::kPi); // [-0.5,0.5]
  if (u < 0.0) u += 1.0;

  // v=1 at +Y pole, v=0 at -Y pole (matches Mesh::makeUvSphere())
  const double vy = math::clamp(ny, -1.0, 1.0);
  const double v = 0.5 + std::asin(vy) / math::kPi;

  outU = static_cast<float>(u);
  outV = static_cast<float>(v);
}

struct TempVert {
  math::Vec3d p;
  math::Vec3d n;
  float u{0.0f};
  float v{0.0f};
};

static inline float sampleField(const ScalarField3D& field, const math::Vec3d& p) {
  return field(static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z));
}

static inline math::Vec3d fieldNormal(const ScalarField3D& field, const math::Vec3d& p, double eps) {
  const math::Vec3d ex{eps, 0.0, 0.0};
  const math::Vec3d ey{0.0, eps, 0.0};
  const math::Vec3d ez{0.0, 0.0, eps};

  const double dx = (double)sampleField(field, p + ex) - (double)sampleField(field, p - ex);
  const double dy = (double)sampleField(field, p + ey) - (double)sampleField(field, p - ey);
  const double dz = (double)sampleField(field, p + ez) - (double)sampleField(field, p - ez);
  return safeNormalize(math::Vec3d{dx, dy, dz});
}

// Marching Tetrahedra tables (edge flags + triangulation), using the standard
// edge ordering for a tetrahedron with vertices 0..3.
// Edges:
//  0: 0-1
//  1: 1-2
//  2: 2-0
//  3: 0-3
//  4: 1-3
//  5: 2-3
static constexpr int kTetEdgeConnection[6][2] = {
    {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3},
};

static constexpr int kTetEdgeFlags[16] = {
    0x00, 0x0d, 0x13, 0x1e, 0x26, 0x2b, 0x35, 0x38,
    0x38, 0x35, 0x2b, 0x26, 0x1e, 0x13, 0x0d, 0x00,
};

// For each case, list 0-2 triangles as edge triples, terminated by -1.
static constexpr int kTetTriangles[16][7] = {
    {-1, -1, -1, -1, -1, -1, -1},
    {0, 3, 2, -1, -1, -1, -1},
    {0, 1, 4, -1, -1, -1, -1},
    {1, 4, 2, 2, 4, 3, -1},

    {1, 2, 5, -1, -1, -1, -1},
    {0, 3, 5, 0, 5, 1, -1},
    {0, 2, 5, 0, 5, 4, -1},
    {5, 4, 3, -1, -1, -1, -1},

    {3, 4, 5, -1, -1, -1, -1},
    {4, 5, 0, 5, 2, 0, -1},
    {1, 5, 0, 5, 3, 0, -1},
    {5, 2, 1, -1, -1, -1, -1},

    {3, 4, 2, 2, 4, 1, -1},
    {4, 1, 0, -1, -1, -1, -1},
    {2, 3, 0, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1},
};

// Cube split into six tetrahedra that share the main diagonal 0-6.
// Corner order for the cube:
//  0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0)
//  4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
static constexpr int kTetsInCube[6][4] = {
    {0, 1, 2, 6},
    {0, 2, 3, 6},
    {0, 3, 7, 6},
    {0, 7, 4, 6},
    {0, 4, 5, 6},
    {0, 5, 1, 6},
};

} // namespace

SdfMeshData meshIsosurfaceMarchingTetrahedra(const ScalarField3D& field,
                                             const SdfMesherParams& params) {
  SdfMeshData out;

  if (!field) {
    return out;
  }

  const int N = std::clamp(params.resolution, 4, 256);
  const float bounds = std::max(1e-4f, params.bounds);
  const float iso = params.iso;

  const float step = 2.0f * bounds / static_cast<float>(N);
  const int n1 = N + 1;

  auto gridIndex = [n1](int i, int j, int k) -> std::size_t {
    return static_cast<std::size_t>((k * n1 + j) * n1 + i);
  };

  // Pre-sample field at grid points so edge intersections are consistent.
  std::vector<float> values(static_cast<std::size_t>(n1) * n1 * n1);
  for (int k = 0; k <= N; ++k) {
    const float z = -bounds + step * static_cast<float>(k);
    for (int j = 0; j <= N; ++j) {
      const float y = -bounds + step * static_cast<float>(j);
      for (int i = 0; i <= N; ++i) {
        const float x = -bounds + step * static_cast<float>(i);
        values[gridIndex(i, j, k)] = field(x, y, z);
      }
    }
  }

  // Conservative reserve: if the field produces very dense geometry this may still grow.
  out.vertices.reserve(static_cast<std::size_t>(N) * N * 3);
  out.indices.reserve(static_cast<std::size_t>(N) * N * 9);

  auto gridSample = [&](int i, int j, int k) -> float {
    return values[gridIndex(i, j, k)];
  };

  // For each cube cell, run 6 tetrahedra.
  for (int k = 0; k < N; ++k) {
    const double z0 = -bounds + step * static_cast<double>(k);
    const double z1 = z0 + step;

    for (int j = 0; j < N; ++j) {
      const double y0 = -bounds + step * static_cast<double>(j);
      const double y1 = y0 + step;

      for (int i = 0; i < N; ++i) {
        const double x0 = -bounds + step * static_cast<double>(i);
        const double x1 = x0 + step;

        // Cube corners.
        math::Vec3d cPos[8] = {
            {x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},
            {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1},
        };

        float cVal[8] = {
            gridSample(i, j, k),
            gridSample(i + 1, j, k),
            gridSample(i + 1, j + 1, k),
            gridSample(i, j + 1, k),
            gridSample(i, j, k + 1),
            gridSample(i + 1, j, k + 1),
            gridSample(i + 1, j + 1, k + 1),
            gridSample(i, j + 1, k + 1),
        };

        for (int t = 0; t < 6; ++t) {
          math::Vec3d p[4];
          float v[4];

          for (int tv = 0; tv < 4; ++tv) {
            const int ci = kTetsInCube[t][tv];
            p[tv] = cPos[ci];
            v[tv] = cVal[ci];
          }

          int caseIndex = 0;
          for (int vi = 0; vi < 4; ++vi) {
            if (v[vi] <= iso) {
              caseIndex |= (1 << vi);
            }
          }

          const int edgeFlags = kTetEdgeFlags[caseIndex];
          if (edgeFlags == 0) {
            continue;
          }

          TempVert edgeVert[6];
          bool edgeValid[6] = {false, false, false, false, false, false};

          // Compute vertices on intersected edges.
          for (int e = 0; e < 6; ++e) {
            if (!(edgeFlags & (1 << e))) continue;

            const int a = kTetEdgeConnection[e][0];
            const int b = kTetEdgeConnection[e][1];
            const float va = v[a];
            const float vb = v[b];

            double t01 = 0.5;
            const double denom = static_cast<double>(vb) - static_cast<double>(va);
            if (std::abs(denom) > 1e-20) {
              t01 = (static_cast<double>(iso) - static_cast<double>(va)) / denom;
            }
            t01 = std::clamp(t01, 0.0, 1.0);

            const math::Vec3d pp = lerp3(p[a], p[b], t01);

            math::Vec3d nn{0.0, 1.0, 0.0};
            if (params.computeNormalsFromField) {
              nn = fieldNormal(field, pp, std::max(1e-8, (double)params.normalEps));
            }

            float uu = 0.0f, vv = 0.0f;
            sphericalUV(pp, uu, vv);

            edgeVert[e] = TempVert{pp, nn, uu, vv};
            edgeValid[e] = true;
          }

          // Emit triangles.
          for (int tri = 0; tri < 2; ++tri) {
            const int e0 = kTetTriangles[caseIndex][tri * 3 + 0];
            if (e0 < 0) break;

            const int e1 = kTetTriangles[caseIndex][tri * 3 + 1];
            const int e2 = kTetTriangles[caseIndex][tri * 3 + 2];
            if (e1 < 0 || e2 < 0) break;

            if (!edgeValid[e0] || !edgeValid[e1] || !edgeValid[e2]) {
              continue;
            }

            TempVert a = edgeVert[e0];
            TempVert b = edgeVert[e1];
            TempVert c = edgeVert[e2];

            // If we have field normals, fix triangle winding so it faces the same way.
            if (params.computeNormalsFromField && params.fixWindingFromNormals) {
              const math::Vec3d triN = math::cross(b.p - a.p, c.p - a.p);
              const math::Vec3d avgN = safeNormalize((a.n + b.n + c.n) * (1.0 / 3.0));
              if (math::dot(triN, avgN) < 0.0) {
                std::swap(b, c);
              }
            }

            const std::uint32_t base = static_cast<std::uint32_t>(out.vertices.size());

            auto push = [&](const TempVert& tv) {
              VertexPNUT vtx{};
              vtx.px = static_cast<float>(tv.p.x);
              vtx.py = static_cast<float>(tv.p.y);
              vtx.pz = static_cast<float>(tv.p.z);
              vtx.nx = static_cast<float>(tv.n.x);
              vtx.ny = static_cast<float>(tv.n.y);
              vtx.nz = static_cast<float>(tv.n.z);
              vtx.u = tv.u;
              vtx.v = tv.v;
              out.vertices.push_back(vtx);
            };

            push(a);
            push(b);
            push(c);

            out.indices.push_back(base + 0);
            out.indices.push_back(base + 1);
            out.indices.push_back(base + 2);
          }
        }
      }
    }
  }

  return out;
}

SdfMeshStats measureSdfMesh(const SdfMeshData& mesh) {
  SdfMeshStats s{};
  s.vertexCount = mesh.vertices.size();
  s.indexCount = mesh.indices.size();

  if (mesh.vertices.empty()) {
    return s;
  }

  float minX = mesh.vertices[0].px;
  float minY = mesh.vertices[0].py;
  float minZ = mesh.vertices[0].pz;
  float maxX = mesh.vertices[0].px;
  float maxY = mesh.vertices[0].py;
  float maxZ = mesh.vertices[0].pz;

  for (const auto& v : mesh.vertices) {
    minX = std::min(minX, v.px);
    minY = std::min(minY, v.py);
    minZ = std::min(minZ, v.pz);
    maxX = std::max(maxX, v.px);
    maxY = std::max(maxY, v.py);
    maxZ = std::max(maxZ, v.pz);
  }

  s.minX = minX;
  s.minY = minY;
  s.minZ = minZ;
  s.maxX = maxX;
  s.maxY = maxY;
  s.maxZ = maxZ;
  return s;
}

} // namespace stellar::render
