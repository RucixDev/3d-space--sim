#include "stellar/render/MeshSimplify.h"

#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <unordered_set>
#include <vector>

namespace stellar::render {

namespace {

// --- Helpers ---------------------------------------------------------------

static inline math::Vec3d toVec3d(const VertexPNUT& v) {
  return { (double)v.px, (double)v.py, (double)v.pz };
}

static inline void setPos(VertexPNUT& v, const math::Vec3d& p) {
  v.px = (float)p.x;
  v.py = (float)p.y;
  v.pz = (float)p.z;
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

  // v=1 at +Y pole, v=0 at -Y pole.
  const double vy = math::clamp(ny, -1.0, 1.0);
  const double v = 0.5 + std::asin(vy) / math::kPi;

  outU = (float)u;
  outV = (float)v;
}

// --- Quadric Error Metric --------------------------------------------------

// Symmetric 4x4 quadric stored in 10 unique coefficients.
//
// Ordering:
//  m0  = q00
//  m1  = q01
//  m2  = q02
//  m3  = q03
//  m4  = q11
//  m5  = q12
//  m6  = q13
//  m7  = q22
//  m8  = q23
//  m9  = q33
struct Quadric {
  std::array<double, 10> m{};

  Quadric& operator+=(const Quadric& o) {
    for (int i = 0; i < 10; ++i) m[(std::size_t)i] += o.m[(std::size_t)i];
    return *this;
  }

  static Quadric fromPlane(double a, double b, double c, double d, double w = 1.0) {
    // K = w * [a b c d]^T [a b c d]
    Quadric q;
    const double aa = a * a * w;
    const double ab = a * b * w;
    const double ac = a * c * w;
    const double ad = a * d * w;
    const double bb = b * b * w;
    const double bc = b * c * w;
    const double bd = b * d * w;
    const double cc = c * c * w;
    const double cd = c * d * w;
    const double dd = d * d * w;

    q.m = {aa, ab, ac, ad, bb, bc, bd, cc, cd, dd};
    return q;
  }

  double evaluate(const math::Vec3d& p) const {
    const double x = p.x;
    const double y = p.y;
    const double z = p.z;

    const double q00 = m[0];
    const double q01 = m[1];
    const double q02 = m[2];
    const double q03 = m[3];
    const double q11 = m[4];
    const double q12 = m[5];
    const double q13 = m[6];
    const double q22 = m[7];
    const double q23 = m[8];
    const double q33 = m[9];

    // v^T Q v with v=[x y z 1]
    const double v0 = q00 * x + q01 * y + q02 * z + q03;
    const double v1 = q01 * x + q11 * y + q12 * z + q13;
    const double v2 = q02 * x + q12 * y + q22 * z + q23;
    const double v3 = q03 * x + q13 * y + q23 * z + q33;

    return x * v0 + y * v1 + z * v2 + 1.0 * v3;
  }
};

static bool solve3x3(double a[3][3], double b[3], double outX[3]) {
  // Gaussian elimination with partial pivoting on an augmented 3x4 matrix.
  double m[3][4] = {
      {a[0][0], a[0][1], a[0][2], b[0]},
      {a[1][0], a[1][1], a[1][2], b[1]},
      {a[2][0], a[2][1], a[2][2], b[2]},
  };

  for (int col = 0; col < 3; ++col) {
    int piv = col;
    double best = std::fabs(m[col][col]);
    for (int r = col + 1; r < 3; ++r) {
      const double v = std::fabs(m[r][col]);
      if (v > best) {
        best = v;
        piv = r;
      }
    }
    if (best < 1e-18) return false;

    if (piv != col) {
      for (int c = col; c < 4; ++c) std::swap(m[col][c], m[piv][c]);
    }

    const double inv = 1.0 / m[col][col];
    for (int c = col; c < 4; ++c) m[col][c] *= inv;

    for (int r = 0; r < 3; ++r) {
      if (r == col) continue;
      const double f = m[r][col];
      if (std::fabs(f) < 1e-20) continue;
      for (int c = col; c < 4; ++c) m[r][c] -= f * m[col][c];
    }
  }

  outX[0] = m[0][3];
  outX[1] = m[1][3];
  outX[2] = m[2][3];
  return true;
}

struct VertexW {
  math::Vec3d p;
  Quadric q;
  bool alive{true};
  std::uint32_t version{0};
};

struct Tri {
  int a{0}, b{0}, c{0};
  bool alive{true};
};

struct EdgeCand {
  int v0{0};
  int v1{0};
  double cost{0.0};
  math::Vec3d pos;
  std::uint32_t ver0{0};
  std::uint32_t ver1{0};
};

struct EdgeCompare {
  bool operator()(const EdgeCand& a, const EdgeCand& b) const {
    return a.cost > b.cost; // min-heap
  }
};

static inline std::uint64_t edgeKey(int a, int b) {
  const std::uint32_t u = (std::uint32_t)std::min(a, b);
  const std::uint32_t v = (std::uint32_t)std::max(a, b);
  return (std::uint64_t(u) << 32) | std::uint64_t(v);
}

static bool computeEdgeCandidate(const std::vector<VertexW>& verts,
                                int a,
                                int b,
                                EdgeCand& out) {
  if (a == b) return false;
  if (a < 0 || b < 0 || a >= (int)verts.size() || b >= (int)verts.size()) return false;
  if (!verts[(std::size_t)a].alive || !verts[(std::size_t)b].alive) return false;

  const VertexW& va = verts[(std::size_t)a];
  const VertexW& vb = verts[(std::size_t)b];
  Quadric q = va.q;
  q += vb.q;

  // Solve for optimal position: A x = -b
  // A = upper-left 3x3, b = [q03,q13,q23]
  double A[3][3] = {
      {q.m[0], q.m[1], q.m[2]},
      {q.m[1], q.m[4], q.m[5]},
      {q.m[2], q.m[5], q.m[7]},
  };
  double bb[3] = {-q.m[3], -q.m[6], -q.m[8]};

  math::Vec3d bestPos;
  double bestCost = std::numeric_limits<double>::infinity();

  double x[3];
  if (solve3x3(A, bb, x)) {
    bestPos = math::Vec3d{x[0], x[1], x[2]};
    bestCost = q.evaluate(bestPos);
  }

  // Fallbacks: endpoints + midpoint
  const math::Vec3d pa = va.p;
  const math::Vec3d pb = vb.p;
  const math::Vec3d pm = (pa + pb) * 0.5;
  const double ca = q.evaluate(pa);
  const double cb = q.evaluate(pb);
  const double cm = q.evaluate(pm);
  if (ca < bestCost) { bestCost = ca; bestPos = pa; }
  if (cb < bestCost) { bestCost = cb; bestPos = pb; }
  if (cm < bestCost) { bestCost = cm; bestPos = pm; }

  out.v0 = std::min(a, b);
  out.v1 = std::max(a, b);
  out.pos = bestPos;
  out.cost = bestCost;
  out.ver0 = verts[(std::size_t)out.v0].version;
  out.ver1 = verts[(std::size_t)out.v1].version;
  return true;
}

static void gatherNeighbors(const std::vector<VertexW>& verts,
                            const std::vector<Tri>& tris,
                            const std::vector<std::vector<int>>& vertTris,
                            int v,
                            std::vector<int>& outNeighbors) {
  outNeighbors.clear();
  if (v < 0 || v >= (int)verts.size()) return;
  if (!verts[(std::size_t)v].alive) return;

  std::unordered_set<int> set;
  set.reserve(64);

  for (int ti : vertTris[(std::size_t)v]) {
    if (ti < 0 || ti >= (int)tris.size()) continue;
    const Tri& t = tris[(std::size_t)ti];
    if (!t.alive) continue;
    const int vs[3] = {t.a, t.b, t.c};
    for (int k = 0; k < 3; ++k) {
      const int u = vs[k];
      if (u == v) continue;
      if (u < 0 || u >= (int)verts.size()) continue;
      if (!verts[(std::size_t)u].alive) continue;
      set.insert(u);
    }
  }

  outNeighbors.reserve(set.size());
  for (int u : set) outNeighbors.push_back(u);
}

static void recomputeNormalsAndUVs(SdfMeshData& mesh,
                                  bool recomputeNormals,
                                  bool recomputeSphericalUVs) {
  if (mesh.vertices.empty() || mesh.indices.empty()) return;

  // Zero normals
  if (recomputeNormals) {
    for (auto& v : mesh.vertices) {
      v.nx = v.ny = v.nz = 0.0f;
    }

    for (std::size_t i = 0; i + 2 < mesh.indices.size(); i += 3) {
      const std::uint32_t ia = mesh.indices[i + 0];
      const std::uint32_t ib = mesh.indices[i + 1];
      const std::uint32_t ic = mesh.indices[i + 2];
      if (ia >= mesh.vertices.size() || ib >= mesh.vertices.size() || ic >= mesh.vertices.size()) continue;

      const math::Vec3d a = toVec3d(mesh.vertices[ia]);
      const math::Vec3d b = toVec3d(mesh.vertices[ib]);
      const math::Vec3d c = toVec3d(mesh.vertices[ic]);

      const math::Vec3d n = math::cross(b - a, c - a);
      const double len = n.length();
      if (len < 1e-18) continue;
      const math::Vec3d nn = n / len;

      mesh.vertices[ia].nx += (float)nn.x;
      mesh.vertices[ia].ny += (float)nn.y;
      mesh.vertices[ia].nz += (float)nn.z;
      mesh.vertices[ib].nx += (float)nn.x;
      mesh.vertices[ib].ny += (float)nn.y;
      mesh.vertices[ib].nz += (float)nn.z;
      mesh.vertices[ic].nx += (float)nn.x;
      mesh.vertices[ic].ny += (float)nn.y;
      mesh.vertices[ic].nz += (float)nn.z;
    }

    for (auto& v : mesh.vertices) {
      const math::Vec3d n{(double)v.nx, (double)v.ny, (double)v.nz};
      const math::Vec3d nn = n.normalized(1e-18);
      v.nx = (float)nn.x;
      v.ny = (float)nn.y;
      v.nz = (float)nn.z;
    }
  }

  if (recomputeSphericalUVs) {
    for (auto& v : mesh.vertices) {
      float u = 0.0f, w = 0.0f;
      sphericalUV(toVec3d(v), u, w);
      v.u = u;
      v.v = w;
    }
  }
}

} // namespace

SdfMeshData simplifyMeshQEM(const SdfMeshData& mesh,
                           const MeshSimplifyParams& params,
                           MeshSimplifyStats* outStats,
                           std::string* outError) {
  const auto t0 = std::chrono::high_resolution_clock::now();

  if (outError) outError->clear();
  if (outStats) *outStats = {};

  if (mesh.vertices.empty() || mesh.indices.empty()) {
    if (outError) *outError = "Empty mesh.";
    return mesh;
  }
  if (mesh.indices.size() % 3 != 0) {
    if (outError) *outError = "Index buffer is not a multiple of 3.";
    return mesh;
  }

  const std::size_t inVerts = mesh.vertices.size();
  const std::size_t inTris = mesh.indices.size() / 3u;

  if (outStats) {
    outStats->inputVertices = inVerts;
    outStats->inputTriangles = inTris;
  }

  // Determine target triangle count.
  std::size_t targetTris = inTris;
  if (params.targetTriangleCount > 0) {
    targetTris = (std::size_t)std::max(1, params.targetTriangleCount);
  } else {
    const double r = std::clamp((double)params.targetTriangleRatio, 0.0, 1.0);
    targetTris = (std::size_t)std::max<std::size_t>(1, (std::size_t)std::llround((double)inTris * r));
  }

  if (targetTris >= inTris || params.targetTriangleRatio >= 0.999f) {
    // Nothing to do.
    if (outStats) {
      outStats->outputVertices = inVerts;
      outStats->outputTriangles = inTris;
      const auto t1 = std::chrono::high_resolution_clock::now();
      outStats->ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    }
    return mesh;
  }

  // Build working vertices.
  std::vector<VertexW> verts(inVerts);
  for (std::size_t i = 0; i < inVerts; ++i) {
    verts[i].p = toVec3d(mesh.vertices[i]);
    verts[i].alive = true;
  }

  // Build triangles + adjacency.
  std::vector<Tri> tris(inTris);
  std::vector<std::vector<int>> vertTris(inVerts);
  vertTris.shrink_to_fit();

  for (std::size_t ti = 0; ti < inTris; ++ti) {
    const std::uint32_t ia = mesh.indices[ti * 3 + 0];
    const std::uint32_t ib = mesh.indices[ti * 3 + 1];
    const std::uint32_t ic = mesh.indices[ti * 3 + 2];
    if (ia >= inVerts || ib >= inVerts || ic >= inVerts) {
      if (outError) *outError = "Index out of range while building triangle list.";
      return mesh;
    }

    Tri t;
    t.a = (int)ia;
    t.b = (int)ib;
    t.c = (int)ic;
    t.alive = true;
    tris[ti] = t;

    vertTris[ia].push_back((int)ti);
    vertTris[ib].push_back((int)ti);
    vertTris[ic].push_back((int)ti);
  }

  // Accumulate quadrics from triangle planes.
  for (const Tri& t : tris) {
    const math::Vec3d a = verts[(std::size_t)t.a].p;
    const math::Vec3d b = verts[(std::size_t)t.b].p;
    const math::Vec3d c = verts[(std::size_t)t.c].p;
    const math::Vec3d n = math::cross(b - a, c - a);
    const double area2 = n.length();
    if (area2 < 1e-18) continue;
    const math::Vec3d nn = n / area2;
    const double d = -math::dot(nn, a);

    // Weight by triangle area so large faces matter more.
    const double w = area2;
    const Quadric q = Quadric::fromPlane(nn.x, nn.y, nn.z, d, w);
    verts[(std::size_t)t.a].q += q;
    verts[(std::size_t)t.b].q += q;
    verts[(std::size_t)t.c].q += q;
  }

  // Build initial edge heap.
  std::priority_queue<EdgeCand, std::vector<EdgeCand>, EdgeCompare> heap;
  std::unordered_set<std::uint64_t> seen;
  seen.reserve(inTris * 3);

  auto pushEdge = [&](int a, int b) {
    if (a == b) return;
    const std::uint64_t key = edgeKey(a, b);
    if (seen.find(key) != seen.end()) return;
    seen.insert(key);

    EdgeCand e;
    if (computeEdgeCandidate(verts, a, b, e)) {
      heap.push(e);
    }
  };

  for (const Tri& t : tris) {
    pushEdge(t.a, t.b);
    pushEdge(t.b, t.c);
    pushEdge(t.c, t.a);
  }

  std::size_t aliveTris = inTris;
  double maxAccepted = 0.0;
  int collapses = 0;

  std::vector<int> neighbors;
  neighbors.reserve(64);

  while (aliveTris > targetTris && !heap.empty() && collapses < params.maxIterations) {
    EdgeCand e = heap.top();
    heap.pop();

    if (e.v0 < 0 || e.v1 < 0 || e.v0 >= (int)verts.size() || e.v1 >= (int)verts.size()) continue;
    VertexW& v0 = verts[(std::size_t)e.v0];
    VertexW& v1 = verts[(std::size_t)e.v1];
    if (!v0.alive || !v1.alive) continue;

    // Stale candidate? recompute.
    if (v0.version != e.ver0 || v1.version != e.ver1) {
      EdgeCand ne;
      if (computeEdgeCandidate(verts, e.v0, e.v1, ne)) {
        heap.push(ne);
      }
      continue;
    }

    if (params.maxError > 0.0 && e.cost > params.maxError) {
      break; // No more acceptable collapses.
    }

    // Collapse v1 into v0 (keep lower index for determinism).
    v0.p = e.pos;
    v0.q += v1.q;
    v0.version++;
    v1.alive = false;
    v1.version++;

    // Rewire triangles that reference v1.
    for (int ti : vertTris[(std::size_t)e.v1]) {
      if (ti < 0 || ti >= (int)tris.size()) continue;
      Tri& t = tris[(std::size_t)ti];
      if (!t.alive) continue;

      if (t.a == e.v1) t.a = e.v0;
      if (t.b == e.v1) t.b = e.v0;
      if (t.c == e.v1) t.c = e.v0;

      if (t.a == t.b || t.b == t.c || t.c == t.a) {
        t.alive = false;
        if (aliveTris > 0) --aliveTris;
        continue;
      }

      // Add to v0's incident list (duplicates are acceptable; we filter later).
      vertTris[(std::size_t)e.v0].push_back(ti);
    }

    vertTris[(std::size_t)e.v1].clear();

    maxAccepted = std::max(maxAccepted, e.cost);
    ++collapses;

    // Add new edge candidates around v0.
    gatherNeighbors(verts, tris, vertTris, e.v0, neighbors);
    for (int nb : neighbors) {
      EdgeCand ne;
      if (computeEdgeCandidate(verts, e.v0, nb, ne)) {
        heap.push(ne);
      }
    }
  }

  // Build compacted output mesh.
  std::vector<bool> used(inVerts, false);
  used.shrink_to_fit();

  std::size_t usedCount = 0;
  for (const Tri& t : tris) {
    if (!t.alive) continue;
    used[(std::size_t)t.a] = true;
    used[(std::size_t)t.b] = true;
    used[(std::size_t)t.c] = true;
  }
  for (std::size_t i = 0; i < inVerts; ++i) {
    if (verts[i].alive && used[i]) ++usedCount;
  }

  std::vector<int> remap(inVerts, -1);
  remap.shrink_to_fit();

  SdfMeshData out;
  out.vertices.reserve(usedCount);
  out.indices.reserve(aliveTris * 3u);

  int next = 0;
  for (std::size_t i = 0; i < inVerts; ++i) {
    if (!verts[i].alive || !used[i]) continue;
    remap[i] = next++;

    VertexPNUT v{};
    setPos(v, verts[i].p);
    v.nx = 0.0f; v.ny = 1.0f; v.nz = 0.0f;
    v.u = 0.0f; v.v = 0.0f;
    out.vertices.push_back(v);
  }

  std::size_t outAliveTris = 0;
  for (const Tri& t : tris) {
    if (!t.alive) continue;
    const int ra = remap[(std::size_t)t.a];
    const int rb = remap[(std::size_t)t.b];
    const int rc = remap[(std::size_t)t.c];
    if (ra < 0 || rb < 0 || rc < 0) continue;
    if (ra == rb || rb == rc || rc == ra) continue;
    out.indices.push_back((std::uint32_t)ra);
    out.indices.push_back((std::uint32_t)rb);
    out.indices.push_back((std::uint32_t)rc);
    ++outAliveTris;
  }

  // Final attribute recompute.
  recomputeNormalsAndUVs(out, params.recomputeNormals, params.recomputeSphericalUVs);

  const auto t1 = std::chrono::high_resolution_clock::now();
  if (outStats) {
    outStats->outputVertices = out.vertices.size();
    outStats->outputTriangles = out.indices.size() / 3u;
    outStats->ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    outStats->maxAcceptedError = maxAccepted;
  }

  // If we somehow produced an empty mesh, fail gracefully.
  if (out.vertices.empty() || out.indices.empty()) {
    if (outError) *outError = "Simplifier produced an empty mesh.";
    return mesh;
  }

  return out;
}

} // namespace stellar::render
