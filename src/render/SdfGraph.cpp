#include "stellar/render/SdfGraph.h"

#include "stellar/core/Hash.h"
#include "stellar/math/Math.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace stellar::render {

namespace {

constexpr float kFar = 1e9f;
constexpr float kNegFar = -1e9f;

struct V3 {
  float x;
  float y;
  float z;
};

struct V2 {
  float x;
  float y;
};

static inline V3 v3(float x, float y, float z) { return {x, y, z}; }
static inline V2 v2(float x, float y) { return {x, y}; }

static inline V3 add(const V3& a, const V3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline V3 sub(const V3& a, const V3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline V3 mul(const V3& a, float s) { return {a.x * s, a.y * s, a.z * s}; }

static inline V3 abs3(const V3& v) { return {std::fabs(v.x), std::fabs(v.y), std::fabs(v.z)}; }
static inline V3 max3(const V3& v, float s) { return {std::max(v.x, s), std::max(v.y, s), std::max(v.z, s)}; }

static inline float dot3(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline float length3(const V3& v) { return std::sqrt(dot3(v, v)); }

static inline float dot2(const V2& a, const V2& b) { return a.x * b.x + a.y * b.y; }
static inline float length2(const V2& v) { return std::sqrt(dot2(v, v)); }

static inline float clampf(float v, float lo, float hi) { return std::max(lo, std::min(v, hi)); }

// --- SDF primitives ---

static inline float sdfSphere(const V3& p, float r) {
  return length3(p) - r;
}

static inline float sdfBox(const V3& p, const V3& b) {
  // b is half-size
  const V3 q = sub(abs3(p), b);
  const float outside = length3(max3(q, 0.0f));
  const float inside = std::min(std::max(q.x, std::max(q.y, q.z)), 0.0f);
  return outside + inside;
}

static inline float sdfCapsule(const V3& p, const V3& a, const V3& b, float r) {
  // distance to segment AB
  const V3 pa = sub(p, a);
  const V3 ba = sub(b, a);
  const float baba = dot3(ba, ba);
  if (baba <= 1e-12f) return length3(pa) - r;
  const float h = clampf(dot3(pa, ba) / baba, 0.0f, 1.0f);
  const V3 d = sub(pa, mul(ba, h));
  return length3(d) - r;
}

static inline float sdfTorusY(const V3& p, float majorR, float minorR) {
  const V2 q = v2(length2(v2(p.x, p.z)) - majorR, p.y);
  return length2(q) - minorR;
}

static inline float sdfPlane(const V3& p, const V3& n, float d) {
  // n should be normalized.
  return dot3(p, n) + d;
}

// Polynomial smooth union (Inigo Quilez style).
static inline float opSmoothUnion(float d1, float d2, float k) {
  k = std::max(k, 1e-6f);
  const float h = clampf(0.5f + 0.5f * (d2 - d1) / k, 0.0f, 1.0f);
  // mix(d2, d1, h) - k*h*(1-h)
  const float m = d2 + (d1 - d2) * h;
  return m - k * h * (1.0f - h);
}

static inline double fbmAmpSum(int octaves, double gain) {
  octaves = std::max(1, octaves);
  double amp = 0.5;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp;
    amp *= gain;
  }
  return sum;
}

static inline float signedFbm01(core::u64 seed,
                               double x, double y, double z,
                               int octaves,
                               double lacunarity,
                               double gain) {
  const double s = proc::fbm3D(seed, x, y, z, octaves, lacunarity, gain);
  const double aSum = fbmAmpSum(octaves, gain);
  const double n01 = (aSum > 0.0) ? (s / aSum) : 0.0;
  const double nClamped = std::max(0.0, std::min(1.0, n01));
  return static_cast<float>(nClamped * 2.0 - 1.0);
}

static inline float signedFbmPerlin01(core::u64 seed,
                                     double x, double y, double z,
                                     int octaves,
                                     double lacunarity,
                                     double gain) {
  const double s = proc::fbmPerlin3D(seed, x, y, z, octaves, lacunarity, gain);
  const double aSum = fbmAmpSum(octaves, gain);
  const double n01 = (aSum > 0.0) ? (s / aSum) : 0.0;
  const double nClamped = std::max(0.0, std::min(1.0, n01));
  return static_cast<float>(nClamped * 2.0 - 1.0);
}

constexpr float kPi = 3.14159265358979323846f;

static inline float degToRad(float deg) { return deg * (kPi / 180.0f); }

static inline V3 rotX(const V3& p, float a) {
  const float c = std::cos(a);
  const float s = std::sin(a);
  return {p.x, c * p.y - s * p.z, s * p.y + c * p.z};
}

static inline V3 rotY(const V3& p, float a) {
  const float c = std::cos(a);
  const float s = std::sin(a);
  return {c * p.x + s * p.z, p.y, -s * p.x + c * p.z};
}

static inline V3 rotZ(const V3& p, float a) {
  const float c = std::cos(a);
  const float s = std::sin(a);
  return {c * p.x - s * p.y, s * p.x + c * p.y, p.z};
}

// Wrap `x` into the range [-period/2, +period/2).
static inline float wrapPeriod(float x, float period) {
  period = std::max(period, 1e-6f);
  const float half = 0.5f * period;
  float y = std::fmod(x + half, period);
  if (y < 0.0f) y += period;
  return y - half;
}

static float evalNodeRec(const SdfGraph& g, int idx, const V3& p, int depth) {
  const int n = std::min((int)g.nodes.size(), kSdfGraphMaxNodes);
  if (idx < 0 || idx >= n) return kFar;

  // Cycle guard: well-formed graphs should be DAGs; if not, bail out safely.
  if (depth > n + 4) return kFar;

  const SdfNode& node = g.nodes[(std::size_t)idx];

  auto evalSame = [&](int child, float def) -> float {
    if (child < 0 || child >= n) return def;
    return evalNodeRec(g, child, p, depth + 1);
  };

  auto evalAt = [&](int child, const V3& pp, float def) -> float {
    if (child < 0 || child >= n) return def;
    return evalNodeRec(g, child, pp, depth + 1);
  };

  switch (node.op) {
    case SdfNodeOp::Constant: {
      return node.p0;
    }

    case SdfNodeOp::Sphere: {
      const V3 c = v3(node.p0, node.p1, node.p2);
      const float r = std::max(0.0f, node.p3);
      return sdfSphere(sub(p, c), r);
    }

    case SdfNodeOp::Box: {
      const V3 c = v3(node.p0, node.p1, node.p2);
      const V3 b = v3(std::max(0.0f, node.p3),
                      std::max(0.0f, node.p4),
                      std::max(0.0f, node.p5));
      return sdfBox(sub(p, c), b);
    }

    case SdfNodeOp::Capsule: {
      const V3 a = v3(node.p0, node.p1, node.p2);
      const V3 b = v3(node.p3, node.p4, node.p5);
      const float r = std::max(0.0f, node.p6);
      return sdfCapsule(p, a, b, r);
    }

    case SdfNodeOp::TorusY: {
      const V3 c = v3(node.p0, node.p1, node.p2);
      const float R = std::max(0.0f, node.p3);
      const float r = std::max(0.0f, node.p4);
      return sdfTorusY(sub(p, c), R, r);
    }

    case SdfNodeOp::Plane: {
      V3 nrm = v3(node.p0, node.p1, node.p2);
      const float len = length3(nrm);
      if (len > 1e-6f) nrm = mul(nrm, 1.0f / len);
      else nrm = v3(0.0f, 1.0f, 0.0f);
      return sdfPlane(p, nrm, node.p3);
    }

    case SdfNodeOp::Union: {
      const float a = evalSame(node.in0, kFar);
      const float b = evalSame(node.in1, kFar);
      return std::min(a, b);
    }

    case SdfNodeOp::SmoothUnion: {
      const float a = evalSame(node.in0, kFar);
      const float b = evalSame(node.in1, kFar);
      return opSmoothUnion(a, b, node.p0);
    }

    case SdfNodeOp::Intersect: {
      const float a = evalSame(node.in0, kNegFar);
      const float b = evalSame(node.in1, kNegFar);
      return std::max(a, b);
    }

    case SdfNodeOp::Subtract: {
      const float a = evalSame(node.in0, kFar);
      const float b = evalSame(node.in1, kFar);
      return std::max(a, -b);
    }

    case SdfNodeOp::NoiseDisplace: {
      const float d = evalSame(node.in0, kFar);
      const float freq = std::max(0.0f, node.p0);
      const float amp = node.p1;
      const int oct = std::max(1, std::min(12, (int)std::lround(node.p2)));
      const double lac = (node.p3 <= 0.0f) ? 2.0 : (double)node.p3;
      const double gain = (node.p4 <= 0.0f) ? 0.5 : (double)node.p4;

      const double ox = (double)node.p5;
      const double oy = (double)node.p6;
      const double oz = (double)node.p7;

      const core::u64 ns = core::hashCombine(g.seed, node.seed);

      const float nSigned = signedFbm01(ns,
                                       (double)p.x * (double)freq + ox,
                                       (double)p.y * (double)freq + oy,
                                       (double)p.z * (double)freq + oz,
                                       oct, lac, gain);

      return d + amp * nSigned;
    }

    case SdfNodeOp::NoiseDisplacePerlin: {
      const float d = evalSame(node.in0, kFar);
      const float freq = std::max(0.0f, node.p0);
      const float amp = node.p1;
      const int oct = std::max(1, std::min(12, (int)std::lround(node.p2)));
      const double lac = (node.p3 <= 0.0f) ? 2.0 : (double)node.p3;
      const double gain = (node.p4 <= 0.0f) ? 0.5 : (double)node.p4;

      const double ox = (double)node.p5;
      const double oy = (double)node.p6;
      const double oz = (double)node.p7;

      const core::u64 ns = core::hashCombine(g.seed, node.seed);

      const float nSigned = signedFbmPerlin01(ns,
                                             (double)p.x * (double)freq + ox,
                                             (double)p.y * (double)freq + oy,
                                             (double)p.z * (double)freq + oz,
                                             oct, lac, gain);

      return d + amp * nSigned;
    }

    case SdfNodeOp::Shell: {
      const float d = evalSame(node.in0, kFar);
      const float t = std::max(0.0f, node.p0);
      return std::fabs(d) - t;
    }

    // --- Space transforms / domain operations ---
    case SdfNodeOp::Translate: {
      const V3 t = v3(node.p0, node.p1, node.p2);
      return evalAt(node.in0, sub(p, t), kFar);
    }

    case SdfNodeOp::RotateX: {
      const float a = degToRad(node.p0);
      const V3 pivot = v3(node.p1, node.p2, node.p3);
      V3 q = sub(p, pivot);
      q = rotX(q, -a);
      q = add(q, pivot);
      return evalAt(node.in0, q, kFar);
    }

    case SdfNodeOp::RotateY: {
      const float a = degToRad(node.p0);
      const V3 pivot = v3(node.p1, node.p2, node.p3);
      V3 q = sub(p, pivot);
      q = rotY(q, -a);
      q = add(q, pivot);
      return evalAt(node.in0, q, kFar);
    }

    case SdfNodeOp::RotateZ: {
      const float a = degToRad(node.p0);
      const V3 pivot = v3(node.p1, node.p2, node.p3);
      V3 q = sub(p, pivot);
      q = rotZ(q, -a);
      q = add(q, pivot);
      return evalAt(node.in0, q, kFar);
    }

    case SdfNodeOp::Scale: {
      const float sc = std::max(1e-6f, std::fabs(node.p0));
      const V3 pivot = v3(node.p1, node.p2, node.p3);
      V3 q = add(mul(sub(p, pivot), 1.0f / sc), pivot);
      const float d = evalAt(node.in0, q, kFar);
      return d * sc;
    }

    case SdfNodeOp::Repeat: {
      V3 period = v3(std::max(1e-6f, node.p0),
                     std::max(1e-6f, node.p1),
                     std::max(1e-6f, node.p2));
      const V3 off = v3(node.p3, node.p4, node.p5);
      V3 q = add(p, off);
      q.x = wrapPeriod(q.x, period.x);
      q.y = wrapPeriod(q.y, period.y);
      q.z = wrapPeriod(q.z, period.z);
      return evalAt(node.in0, q, kFar);
    }

    case SdfNodeOp::Mirror: {
      const bool mx = node.p0 > 0.5f;
      const bool my = node.p1 > 0.5f;
      const bool mz = node.p2 > 0.5f;
      const V3 pivot = v3(node.p3, node.p4, node.p5);

      V3 q = sub(p, pivot);
      if (mx) q.x = std::fabs(q.x);
      if (my) q.y = std::fabs(q.y);
      if (mz) q.z = std::fabs(q.z);
      q = add(q, pivot);

      return evalAt(node.in0, q, kFar);
    }

    case SdfNodeOp::TwistY: {
      const float k = degToRad(node.p0); // deg/unit -> rad/unit
      const V3 pivot = v3(node.p1, node.p2, node.p3);

      V3 q = sub(p, pivot);
      const float a = k * q.y;
      q = rotY(q, -a);
      q = add(q, pivot);

      return evalAt(node.in0, q, kFar);
    }

    default: {
      return kFar;
    }
  }
}

} // namespace

const char* sdfNodeOpName(SdfNodeOp op) {
  switch (op) {
    case SdfNodeOp::Constant: return "Constant";
    case SdfNodeOp::Sphere: return "Sphere";
    case SdfNodeOp::Box: return "Box";
    case SdfNodeOp::Capsule: return "Capsule";
    case SdfNodeOp::TorusY: return "Torus (Y)";
    case SdfNodeOp::Plane: return "Plane";
    case SdfNodeOp::Union: return "Union";
    case SdfNodeOp::SmoothUnion: return "SmoothUnion";
    case SdfNodeOp::Intersect: return "Intersect";
    case SdfNodeOp::Subtract: return "Subtract";
    case SdfNodeOp::NoiseDisplace: return "NoiseDisplace";
    case SdfNodeOp::NoiseDisplacePerlin: return "NoiseDisplace (Perlin)";
    case SdfNodeOp::Shell: return "Shell";
    case SdfNodeOp::Translate: return "Translate";
    case SdfNodeOp::RotateX: return "Rotate (X)";
    case SdfNodeOp::RotateY: return "Rotate (Y)";
    case SdfNodeOp::RotateZ: return "Rotate (Z)";
    case SdfNodeOp::Scale: return "Scale";
    case SdfNodeOp::Repeat: return "Repeat";
    case SdfNodeOp::Mirror: return "Mirror";
    case SdfNodeOp::TwistY: return "Twist (Y)";
    default: return "(unknown)";
  }
}


const char* sdfNodeOpId(SdfNodeOp op) {
  switch (op) {
    case SdfNodeOp::Constant: return "Constant";
    case SdfNodeOp::Sphere: return "Sphere";
    case SdfNodeOp::Box: return "Box";
    case SdfNodeOp::Capsule: return "Capsule";
    case SdfNodeOp::TorusY: return "TorusY";
    case SdfNodeOp::Plane: return "Plane";
    case SdfNodeOp::Union: return "Union";
    case SdfNodeOp::SmoothUnion: return "SmoothUnion";
    case SdfNodeOp::Intersect: return "Intersect";
    case SdfNodeOp::Subtract: return "Subtract";
    case SdfNodeOp::NoiseDisplace: return "NoiseDisplace";
    case SdfNodeOp::NoiseDisplacePerlin: return "NoiseDisplacePerlin";
    case SdfNodeOp::Shell: return "Shell";
    case SdfNodeOp::Translate: return "Translate";
    case SdfNodeOp::RotateX: return "RotateX";
    case SdfNodeOp::RotateY: return "RotateY";
    case SdfNodeOp::RotateZ: return "RotateZ";
    case SdfNodeOp::Scale: return "Scale";
    case SdfNodeOp::Repeat: return "Repeat";
    case SdfNodeOp::Mirror: return "Mirror";
    case SdfNodeOp::TwistY: return "TwistY";
    default: return "Unknown";
  }
}


static std::string lowerAscii(std::string s) {
  for (char& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

std::optional<SdfNodeOp> sdfNodeOpFromString(const std::string& s) {
  const std::string k = lowerAscii(s);
  if (k == "constant") return SdfNodeOp::Constant;
  if (k == "sphere") return SdfNodeOp::Sphere;
  if (k == "box") return SdfNodeOp::Box;
  if (k == "capsule") return SdfNodeOp::Capsule;
  if (k == "torusy" || k == "torus_y" || k == "torus(y)" || k == "torus" || k == "torus (y)") return SdfNodeOp::TorusY;
  if (k == "plane") return SdfNodeOp::Plane;

  if (k == "union") return SdfNodeOp::Union;
  if (k == "smoothunion" || k == "smooth_union" || k == "smooth" || k == "smin") return SdfNodeOp::SmoothUnion;
  if (k == "intersect" || k == "intersection") return SdfNodeOp::Intersect;
  if (k == "subtract" || k == "difference" || k == "sub") return SdfNodeOp::Subtract;

  if (k == "noisedisplace" || k == "noise_displace" || k == "displace" || k == "warp") return SdfNodeOp::NoiseDisplace;
  if (k == "noisedisplaceperlin" || k == "noise_displace_perlin" || k == "perlindisplace" || k == "perlin_displace") return SdfNodeOp::NoiseDisplacePerlin;
  if (k == "shell") return SdfNodeOp::Shell;

  // Space transforms / domain ops.
  if (k == "translate" || k == "move" || k == "offset") return SdfNodeOp::Translate;

  if (k == "rotatex" || k == "rotate_x" || k == "rotx") return SdfNodeOp::RotateX;
  if (k == "rotatey" || k == "rotate_y" || k == "roty") return SdfNodeOp::RotateY;
  if (k == "rotatez" || k == "rotate_z" || k == "rotz") return SdfNodeOp::RotateZ;

  if (k == "scale" || k == "uniformscale" || k == "uniform_scale") return SdfNodeOp::Scale;

  if (k == "repeat" || k == "tile" || k == "tiling" || k == "replicate") return SdfNodeOp::Repeat;

  if (k == "mirror" || k == "symmetry" || k == "sym") return SdfNodeOp::Mirror;

  if (k == "twisty" || k == "twist_y" || k == "twist") return SdfNodeOp::TwistY;

  return std::nullopt;
}


SdfGraph SdfGraph::makeDefault() {
  SdfGraph g;
  g.seed = 0xBADC0FFEEULL;

  // Single sphere at origin.
  SdfNode s;
  s.op = SdfNodeOp::Sphere;
  s.p0 = 0.0f; // cx
  s.p1 = 0.0f; // cy
  s.p2 = 0.0f; // cz
  s.p3 = 1.0f; // r
  g.nodes.push_back(s);
  g.output = 0;
  return g;
}

const char* sdfGraphPresetName(SdfGraphPreset preset) {
  switch (preset) {
    case SdfGraphPreset::Asteroid: return "Asteroid";
    case SdfGraphPreset::Crystal: return "Crystal";
    case SdfGraphPreset::Torus: return "Torus";
    case SdfGraphPreset::HollowBox: return "Hollow Box";
    case SdfGraphPreset::BooleanDemo: return "Booleans Demo";
    default: return "(unknown)";
  }
}

SdfGraph makeSdfGraphPreset(SdfGraphPreset preset, core::u64 seed) {
  SdfGraph g;
  g.seed = seed;

  auto push = [&](const SdfNode& n) {
    if ((int)g.nodes.size() < kSdfGraphMaxNodes) g.nodes.push_back(n);
  };

  if (preset == SdfGraphPreset::Asteroid) {
    // Sphere -> 2-stage noise displacement.
    SdfNode base;
    base.op = SdfNodeOp::Sphere;
    base.p3 = 0.95f;
    push(base);

    SdfNode warp1;
    warp1.op = SdfNodeOp::NoiseDisplace;
    warp1.in0 = 0;
    warp1.p0 = 2.6f;  // freq
    warp1.p1 = 0.22f; // amp
    warp1.p2 = 5.0f;  // octaves
    warp1.p3 = 2.0f;  // lacunarity
    warp1.p4 = 0.52f; // gain
    warp1.seed = core::hashCombine(seed, core::fnv1a64("asteroid_warp1"));
    push(warp1);

    SdfNode warp2;
    warp2.op = SdfNodeOp::NoiseDisplace;
    warp2.in0 = 1;
    warp2.p0 = 8.0f;  // freq
    warp2.p1 = 0.06f; // amp
    warp2.p2 = 3.0f;  // octaves
    warp2.p3 = 2.2f;  // lacunarity
    warp2.p4 = 0.55f; // gain
    warp2.seed = core::hashCombine(seed, core::fnv1a64("asteroid_warp2"));
    push(warp2);

    g.output = 2;
    return g;
  }

  if (preset == SdfGraphPreset::Crystal) {
    // Box ∩ Sphere, then subtle noise to break perfect symmetry.
    SdfNode box;
    box.op = SdfNodeOp::Box;
    box.p0 = 0.0f; box.p1 = 0.0f; box.p2 = 0.0f;
    box.p3 = 0.35f; box.p4 = 0.85f; box.p5 = 0.35f;
    push(box);

    SdfNode sphere;
    sphere.op = SdfNodeOp::Sphere;
    sphere.p3 = 0.95f;
    push(sphere);

    SdfNode isect;
    isect.op = SdfNodeOp::Intersect;
    isect.in0 = 0;
    isect.in1 = 1;
    push(isect);

    SdfNode warp;
    // Perlin-style gradient noise tends to avoid "blocky" grid artifacts.
    warp.op = SdfNodeOp::NoiseDisplacePerlin;
    warp.in0 = 2;
    warp.p0 = 6.0f;
    warp.p1 = 0.045f;
    warp.p2 = 4.0f;
    warp.p3 = 2.0f;
    warp.p4 = 0.55f;
    warp.seed = core::hashCombine(seed, core::fnv1a64("crystal_warp"));
    push(warp);

    g.output = 3;
    return g;
  }

  if (preset == SdfGraphPreset::Torus) {
    SdfNode t;
    t.op = SdfNodeOp::TorusY;
    t.p3 = 0.80f; // major
    t.p4 = 0.24f; // minor
    push(t);

    SdfNode warp;
    warp.op = SdfNodeOp::NoiseDisplacePerlin;
    warp.in0 = 0;
    warp.p0 = 4.0f;
    warp.p1 = 0.04f;
    warp.p2 = 4.0f;
    warp.p3 = 2.0f;
    warp.p4 = 0.5f;
    warp.seed = core::hashCombine(seed, core::fnv1a64("torus_warp"));
    push(warp);

    g.output = 1;
    return g;
  }

  if (preset == SdfGraphPreset::HollowBox) {
    // Outer box minus inner box.
    SdfNode outer;
    outer.op = SdfNodeOp::Box;
    outer.p3 = 0.95f;
    outer.p4 = 0.95f;
    outer.p5 = 0.95f;
    push(outer);

    SdfNode inner;
    inner.op = SdfNodeOp::Box;
    inner.p3 = 0.70f;
    inner.p4 = 0.70f;
    inner.p5 = 0.70f;
    push(inner);

    SdfNode sub;
    sub.op = SdfNodeOp::Subtract;
    sub.in0 = 0;
    sub.in1 = 1;
    push(sub);

    g.output = 2;
    return g;
  }

  // Boolean demo: two overlapping spheres, with derived union/subtract/intersect.
  {
    SdfNode a;
    a.op = SdfNodeOp::Sphere;
    a.p0 = -0.35f;
    a.p3 = 0.75f;
    push(a);

    SdfNode b;
    b.op = SdfNodeOp::Sphere;
    b.p0 = 0.35f;
    b.p3 = 0.75f;
    push(b);

    SdfNode uni;
    uni.op = SdfNodeOp::Union;
    uni.in0 = 0;
    uni.in1 = 1;
    push(uni);

    SdfNode sub;
    sub.op = SdfNodeOp::Subtract;
    sub.in0 = 0;
    sub.in1 = 1;
    push(sub);

    SdfNode isect;
    isect.op = SdfNodeOp::Intersect;
    isect.in0 = 0;
    isect.in1 = 1;
    push(isect);

    g.output = 2; // union by default
    return g;
  }
}

float evalSdfGraph(const SdfGraph& g, float x, float y, float z) {
  const int n = std::min((int)g.nodes.size(), kSdfGraphMaxNodes);
  if (n <= 0) return kFar;
  if (g.output < 0 || g.output >= n) return kFar;

  const V3 p = v3(x, y, z);
  return evalNodeRec(g, g.output, p, 0);
}


bool saveSdfGraphToFile(const SdfGraph& g, const std::string& path, std::string* outError) {
  namespace fs = std::filesystem;

  if (path.empty()) {
    if (outError) *outError = "Empty graph path.";
    return false;
  }

  std::error_code ec;
  fs::path p(path);
  if (p.has_parent_path()) {
    fs::create_directories(p.parent_path(), ec);
  }

  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.is_open()) {
    if (outError) *outError = "Failed to open graph file for writing.";
    return false;
  }

  const int nodeCount = std::min<int>((int)g.nodes.size(), kSdfGraphMaxNodes);

  f.setf(std::ios::fixed);
  f.precision(6);

  f << "StellarForgeSdfGraph 1\n";
  f << "seed " << (unsigned long long)g.seed << "\n";
  f << "output " << g.output << "\n";
  f << "nodes " << nodeCount << "\n";
  for (int i = 0; i < nodeCount; ++i) {
    const SdfNode& n = g.nodes[(std::size_t)i];
    f << "node " << sdfNodeOpId(n.op)
      << " " << n.in0 << " " << n.in1
      << " " << n.p0 << " " << n.p1 << " " << n.p2 << " " << n.p3
      << " " << n.p4 << " " << n.p5 << " " << n.p6 << " " << n.p7
      << " " << (unsigned long long)n.seed << "\n";
  }

  return true;
}

static bool parseI32(const std::string& s, int& out) {
  try {
    std::size_t pos = 0;
    int v = std::stoi(s, &pos, 0);
    if (pos != s.size()) return false;
    out = v;
    return true;
  } catch (...) {
    return false;
  }
}

static bool parseU64(const std::string& s, core::u64& out) {
  try {
    std::size_t pos = 0;
    unsigned long long v = std::stoull(s, &pos, 0);
    if (pos != s.size()) return false;
    out = (core::u64)v;
    return true;
  } catch (...) {
    return false;
  }
}

bool loadSdfGraphFromFile(const std::string& path, SdfGraph& out, std::string* outError) {
  std::ifstream f(path);
  if (!f.is_open()) {
    if (outError) *outError = "Failed to open graph file for reading.";
    return false;
  }

  std::string header;
  int version = 0;
  if (!(f >> header >> version)) {
    if (outError) *outError = "Failed to read graph header.";
    return false;
  }
  if (header != "StellarForgeSdfGraph") {
    if (outError) *outError = "Not a StellarForgeSdfGraph file.";
    return false;
  }

  SdfGraph g{};
  g.seed = 0xBADC0FFEEULL;
  g.output = -1;
  g.nodes.clear();

  (void)version; // Currently only v1.

  std::string line;
  std::getline(f, line); // consume remainder of header line
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    std::istringstream ss(line);
    std::string key;
    if (!(ss >> key)) continue;
    key = lowerAscii(key);

    if (key == "seed") {
      std::string v;
      if (ss >> v) {
        core::u64 tmp = 0;
        if (parseU64(v, tmp)) g.seed = tmp;
      }
    } else if (key == "output") {
      ss >> g.output;
    } else if (key == "node") {
      SdfNode n{};
      std::string opTok;
      if (!(ss >> opTok)) continue;

      if (auto op = sdfNodeOpFromString(opTok)) {
        n.op = *op;
      } else {
        int opInt = 0;
        if (parseI32(opTok, opInt)) {
          n.op = (SdfNodeOp)std::clamp(opInt, 0, 255);
        } else {
          continue;
        }
      }

      ss >> n.in0 >> n.in1;
      ss >> n.p0 >> n.p1 >> n.p2 >> n.p3 >> n.p4 >> n.p5 >> n.p6 >> n.p7;
      std::string seedTok;
      if (ss >> seedTok) {
        core::u64 tmp = 0;
        if (parseU64(seedTok, tmp)) n.seed = tmp;
      }

      if ((int)g.nodes.size() < kSdfGraphMaxNodes) {
        g.nodes.push_back(n);
      }
    }
  }

  // Clamp output.
  if (g.nodes.empty()) {
    g.output = -1;
  } else {
    if (g.output < 0 || g.output >= (int)g.nodes.size()) {
      g.output = (int)g.nodes.size() - 1;
    }
  }

  out = std::move(g);
  return true;
}

ScalarField3D makeSdfField(const SdfGraph& g) {
  // Capture a copy so this can be safely used from a background thread.
  const SdfGraph gg = g;
  return [gg](float x, float y, float z) -> float {
    return evalSdfGraph(gg, x, y, z);
  };
}

} // namespace stellar::render
