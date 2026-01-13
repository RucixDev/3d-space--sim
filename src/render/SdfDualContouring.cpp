#include "stellar/render/SdfDualContouring.h"

#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace stellar::render {

namespace {

static inline math::Vec3d lerp3(const math::Vec3d& a, const math::Vec3d& b, double t) {
  return a + (b - a) * t;
}

static inline math::Vec3d safeNormalize(const math::Vec3d& v, double eps = 1e-18) {
  return v.normalized(eps);
}

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

static bool solve3x3(double a[3][3], double b[3], double outX[3]) {
  // Gaussian elimination with partial pivoting on an augmented 3x4 matrix.
  double m[3][4] = {
      {a[0][0], a[0][1], a[0][2], b[0]},
      {a[1][0], a[1][1], a[1][2], b[1]},
      {a[2][0], a[2][1], a[2][2], b[2]},
  };

  for (int col = 0; col < 3; ++col) {
    // Pivot
    int piv = col;
    double best = std::fabs(m[col][col]);
    for (int r = col + 1; r < 3; ++r) {
      const double v = std::fabs(m[r][col]);
      if (v > best) {
        best = v;
        piv = r;
      }
    }

    if (best < 1e-16) {
      return false;
    }

    if (piv != col) {
      for (int c = col; c < 4; ++c) {
        std::swap(m[col][c], m[piv][c]);
      }
    }

    // Normalize pivot row
    const double inv = 1.0 / m[col][col];
    for (int c = col; c < 4; ++c) {
      m[col][c] *= inv;
    }

    // Eliminate other rows
    for (int r = 0; r < 3; ++r) {
      if (r == col) continue;
      const double f = m[r][col];
      if (std::fabs(f) < 1e-20) continue;
      for (int c = col; c < 4; ++c) {
        m[r][c] -= f * m[col][c];
      }
    }
  }

  outX[0] = m[0][3];
  outX[1] = m[1][3];
  outX[2] = m[2][3];
  return true;
}

static inline void qefAddPlane(double ata[3][3], double atb[3], const math::Vec3d& n, const math::Vec3d& p) {
  const double nx = n.x;
  const double ny = n.y;
  const double nz = n.z;
  const double d = nx * p.x + ny * p.y + nz * p.z;

  // A^T A += n n^T
  ata[0][0] += nx * nx;
  ata[0][1] += nx * ny;
  ata[0][2] += nx * nz;

  ata[1][0] += ny * nx;
  ata[1][1] += ny * ny;
  ata[1][2] += ny * nz;

  ata[2][0] += nz * nx;
  ata[2][1] += nz * ny;
  ata[2][2] += nz * nz;

  // A^T b += n * (n·p)
  atb[0] += nx * d;
  atb[1] += ny * d;
  atb[2] += nz * d;
}

static inline double clampd(double v, double lo, double hi) {
  return std::max(lo, std::min(v, hi));
}

} // namespace

SdfMeshData meshIsosurfaceDualContouring(const ScalarField3D& field,
                                        const SdfDualContourParams& params) {
  SdfMeshData out;

  if (!field) {
    return out;
  }

  const int N = std::clamp(params.resolution, 4, 256);
  const float bounds = std::max(1.0e-4f, params.bounds);
  const float iso = params.iso;
  const double step = 2.0 * (double)bounds / (double)N;
  const int n1 = N + 1;

  auto gridIndex = [n1](int i, int j, int k) -> std::size_t {
    return static_cast<std::size_t>((k * n1 + j) * n1 + i);
  };

  // Pre-sample field at grid points so topology is consistent.
  std::vector<float> values(static_cast<std::size_t>(n1) * n1 * n1);
  for (int k = 0; k <= N; ++k) {
    const float z = -bounds + (float)(step * (double)k);
    for (int j = 0; j <= N; ++j) {
      const float y = -bounds + (float)(step * (double)j);
      for (int i = 0; i <= N; ++i) {
        const float x = -bounds + (float)(step * (double)i);
        values[gridIndex(i, j, k)] = field(x, y, z);
      }
    }
  }

  auto gridSample = [&](int i, int j, int k) -> float {
    return values[gridIndex(i, j, k)];
  };

  auto cellIndex = [N](int i, int j, int k) -> std::size_t {
    return static_cast<std::size_t>((k * N + j) * N + i);
  };

  // Map each cell -> vertex index (or -1 if cell has no surface crossing).
  std::vector<int> cellVert(static_cast<std::size_t>(N) * N * N, -1);

  // Conservative reserves.
  out.vertices.reserve(static_cast<std::size_t>(N) * N * 2);
  out.indices.reserve(static_cast<std::size_t>(N) * N * 6);

  // Cube edge list (corner pairs).
  static constexpr int kEdgePairs[12][2] = {
      {0, 1}, {1, 2}, {2, 3}, {3, 0},
      {4, 5}, {5, 6}, {6, 7}, {7, 4},
      {0, 4}, {1, 5}, {2, 6}, {3, 7},
  };

  // Build per-cell vertices.
  for (int k = 0; k < N; ++k) {
    const double z0 = -(double)bounds + step * (double)k;
    const double z1 = z0 + step;

    for (int j = 0; j < N; ++j) {
      const double y0 = -(double)bounds + step * (double)j;
      const double y1 = y0 + step;

      for (int i = 0; i < N; ++i) {
        const double x0 = -(double)bounds + step * (double)i;
        const double x1 = x0 + step;

        // Corner positions.
        const math::Vec3d cPos[8] = {
            {x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0},
            {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1},
        };

        // Corner values.
        const float cVal[8] = {
            gridSample(i, j, k),
            gridSample(i + 1, j, k),
            gridSample(i + 1, j + 1, k),
            gridSample(i, j + 1, k),
            gridSample(i, j, k + 1),
            gridSample(i + 1, j, k + 1),
            gridSample(i + 1, j + 1, k + 1),
            gridSample(i, j + 1, k + 1),
        };

        bool anyIn = false;
        bool anyOut = false;
        for (int c = 0; c < 8; ++c) {
          const bool inside = (cVal[c] <= iso);
          anyIn |= inside;
          anyOut |= !inside;
        }

        if (!(anyIn && anyOut)) {
          continue; // no crossing in this cell
        }

        // QEF accumulation from Hermite samples.
        double ata[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        double atb[3] = {0, 0, 0};
        math::Vec3d pSum{0, 0, 0};
        int pCount = 0;

        const double eps = std::max(1.0e-6, (double)params.normalEps);

        for (int e = 0; e < 12; ++e) {
          const int a = kEdgePairs[e][0];
          const int b = kEdgePairs[e][1];

          const float va = cVal[a];
          const float vb = cVal[b];
          const bool sa = (va <= iso);
          const bool sb = (vb <= iso);
          if (sa == sb) continue;

          const double denom = (double)vb - (double)va;
          double t = 0.5;
          if (std::fabs(denom) > 1.0e-12) {
            t = ((double)iso - (double)va) / denom;
          }
          t = std::clamp(t, 0.0, 1.0);

          const math::Vec3d p = lerp3(cPos[a], cPos[b], t);
          const math::Vec3d n = fieldNormal(field, p, eps);
          if (n.lengthSq() < 1.0e-18) continue;

          qefAddPlane(ata, atb, n, p);
          pSum += p;
          ++pCount;
        }

        if (pCount == 0) {
          continue;
        }

        const math::Vec3d pAvg = pSum / (double)pCount;

        // Regularization.
        const double reg = std::max(0.0, (double)params.qefRegularization);
        ata[0][0] += reg;
        ata[1][1] += reg;
        ata[2][2] += reg;

        math::Vec3d x = pAvg;

        // Solve (A^T A) x = (A^T b)
        double sol[3] = {0, 0, 0};
        double aCopy[3][3] = {
            {ata[0][0], ata[0][1], ata[0][2]},
            {ata[1][0], ata[1][1], ata[1][2]},
            {ata[2][0], ata[2][1], ata[2][2]},
        };
        double bCopy[3] = {atb[0], atb[1], atb[2]};
        if (solve3x3(aCopy, bCopy, sol)) {
          x = {sol[0], sol[1], sol[2]};
        }

        // Clamp to cell bounds to avoid spikes.
        if (params.clampToCell) {
          x.x = clampd(x.x, x0, x1);
          x.y = clampd(x.y, y0, y1);
          x.z = clampd(x.z, z0, z1);
        }

        // Optional projection onto iso surface.
        if (params.projectToIso && params.projectIterations > 0) {
          for (int it = 0; it < params.projectIterations; ++it) {
            const float f = sampleField(field, x);
            const double d = (double)f - (double)iso;
            if (std::fabs(d) < 1.0e-6) break;

            const math::Vec3d n = fieldNormal(field, x, eps);
            if (n.lengthSq() < 1.0e-18) break;

            // Newton step for signed distance-ish fields: x -= n * (f-iso)
            x -= n * d;

            if (params.clampToCell) {
              x.x = clampd(x.x, x0, x1);
              x.y = clampd(x.y, y0, y1);
              x.z = clampd(x.z, z0, z1);
            }
          }
        }

        VertexPNUT v{};
        v.px = (float)x.x;
        v.py = (float)x.y;
        v.pz = (float)x.z;

        if (params.computeNormalsFromField) {
          const math::Vec3d n = fieldNormal(field, x, eps);
          v.nx = (float)n.x;
          v.ny = (float)n.y;
          v.nz = (float)n.z;
        } else {
          v.nx = 0.0f;
          v.ny = 1.0f;
          v.nz = 0.0f;
        }

        sphericalUV(x, v.u, v.v);

        const int vertIndex = (int)out.vertices.size();
        out.vertices.push_back(v);
        cellVert[cellIndex(i, j, k)] = vertIndex;
      }
    }
  }

  auto getCellVert = [&](int i, int j, int k) -> int {
    if (i < 0 || j < 0 || k < 0 || i >= N || j >= N || k >= N) return -1;
    return cellVert[cellIndex(i, j, k)];
  };

  auto edgeHasCrossing = [&](float a, float b) -> bool {
    const bool sa = (a <= iso);
    const bool sb = (b <= iso);
    return sa != sb;
  };

  // Emit quads around grid edges with sign changes.
  // Each quad is triangulated as (a,b,c) (a,c,d). Winding will be fixed later.

  // X edges: i in [0..N-1], j in [1..N-1], k in [1..N-1]
  for (int k = 1; k < N; ++k) {
    for (int j = 1; j < N; ++j) {
      for (int i = 0; i < N; ++i) {
        const float v0 = gridSample(i, j, k);
        const float v1 = gridSample(i + 1, j, k);
        if (!edgeHasCrossing(v0, v1)) continue;

        const int a = getCellVert(i, j - 1, k - 1);
        const int b = getCellVert(i, j,     k - 1);
        const int c = getCellVert(i, j,     k);
        const int d = getCellVert(i, j - 1, k);
        if (a < 0 || b < 0 || c < 0 || d < 0) continue;

        out.indices.push_back((std::uint32_t)a);
        out.indices.push_back((std::uint32_t)b);
        out.indices.push_back((std::uint32_t)c);

        out.indices.push_back((std::uint32_t)a);
        out.indices.push_back((std::uint32_t)c);
        out.indices.push_back((std::uint32_t)d);
      }
    }
  }

  // Y edges: i in [1..N-1], j in [0..N-1], k in [1..N-1]
  for (int k = 1; k < N; ++k) {
    for (int j = 0; j < N; ++j) {
      for (int i = 1; i < N; ++i) {
        const float v0 = gridSample(i, j, k);
        const float v1 = gridSample(i, j + 1, k);
        if (!edgeHasCrossing(v0, v1)) continue;

        const int a = getCellVert(i - 1, j, k - 1);
        const int b = getCellVert(i,     j, k - 1);
        const int c = getCellVert(i,     j, k);
        const int d = getCellVert(i - 1, j, k);
        if (a < 0 || b < 0 || c < 0 || d < 0) continue;

        out.indices.push_back((std::uint32_t)a);
        out.indices.push_back((std::uint32_t)b);
        out.indices.push_back((std::uint32_t)c);

        out.indices.push_back((std::uint32_t)a);
        out.indices.push_back((std::uint32_t)c);
        out.indices.push_back((std::uint32_t)d);
      }
    }
  }

  // Z edges: i in [1..N-1], j in [1..N-1], k in [0..N-1]
  for (int k = 0; k < N; ++k) {
    for (int j = 1; j < N; ++j) {
      for (int i = 1; i < N; ++i) {
        const float v0 = gridSample(i, j, k);
        const float v1 = gridSample(i, j, k + 1);
        if (!edgeHasCrossing(v0, v1)) continue;

        const int a = getCellVert(i - 1, j - 1, k);
        const int b = getCellVert(i,     j - 1, k);
        const int c = getCellVert(i,     j,     k);
        const int d = getCellVert(i - 1, j,     k);
        if (a < 0 || b < 0 || c < 0 || d < 0) continue;

        out.indices.push_back((std::uint32_t)a);
        out.indices.push_back((std::uint32_t)b);
        out.indices.push_back((std::uint32_t)c);

        out.indices.push_back((std::uint32_t)a);
        out.indices.push_back((std::uint32_t)c);
        out.indices.push_back((std::uint32_t)d);
      }
    }
  }

  // If requested, compute normals even when we didn't during vertex build.
  if (params.computeNormalsFromField) {
    const double eps = std::max(1.0e-6, (double)params.normalEps);
    for (auto& v : out.vertices) {
      const math::Vec3d p{(double)v.px, (double)v.py, (double)v.pz};
      const math::Vec3d n = fieldNormal(field, p, eps);
      v.nx = (float)n.x;
      v.ny = (float)n.y;
      v.nz = (float)n.z;
    }
  }

  // Fix winding by comparing each triangle's geometric normal to the field normal.
  if (params.fixWindingFromNormals && !out.indices.empty()) {
    const double eps = std::max(1.0e-6, (double)params.normalEps);
    for (std::size_t t = 0; t + 2 < out.indices.size(); t += 3) {
      const std::uint32_t ia = out.indices[t + 0];
      const std::uint32_t ib = out.indices[t + 1];
      const std::uint32_t ic = out.indices[t + 2];
      if (ia >= out.vertices.size() || ib >= out.vertices.size() || ic >= out.vertices.size()) continue;

      const VertexPNUT& va = out.vertices[ia];
      const VertexPNUT& vb = out.vertices[ib];
      const VertexPNUT& vc = out.vertices[ic];

      const math::Vec3d pa{(double)va.px, (double)va.py, (double)va.pz};
      const math::Vec3d pb{(double)vb.px, (double)vb.py, (double)vb.pz};
      const math::Vec3d pc{(double)vc.px, (double)vc.py, (double)vc.pz};

      const math::Vec3d triN = safeNormalize(math::cross(pb - pa, pc - pa));
      if (triN.lengthSq() < 1.0e-18) continue;

      // Prefer vertex normals if available; otherwise sample the field.
      math::Vec3d hintN{(double)(va.nx + vb.nx + vc.nx),
                        (double)(va.ny + vb.ny + vc.ny),
                        (double)(va.nz + vb.nz + vc.nz)};
      if (hintN.lengthSq() < 1.0e-12) {
        const math::Vec3d centroid = (pa + pb + pc) * (1.0 / 3.0);
        hintN = fieldNormal(field, centroid, eps);
      } else {
        hintN = safeNormalize(hintN);
      }

      if (math::dot(triN, hintN) < 0.0) {
        std::swap(out.indices[t + 1], out.indices[t + 2]);
      }
    }
  }

  return out;
}

} // namespace stellar::render
