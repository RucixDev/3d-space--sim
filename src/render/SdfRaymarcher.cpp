#include "stellar/render/SdfRaymarcher.h"

#include "stellar/core/Hash.h"
#include "stellar/render/Gl.h"

#include <algorithm>
#include <chrono>
#include <sstream>

namespace stellar::render {

namespace {

int clampNodeCount(int n) {
  return std::clamp(n, 0, kSdfGraphMaxNodes);
}

const char* kFullscreenTriVS = R"GLSL(
#version 330 core
out vec2 vUv;

const vec2 kPos[3] = vec2[](
  vec2(-1.0, -1.0),
  vec2( 3.0, -1.0),
  vec2(-1.0,  3.0)
);

void main() {
  gl_Position = vec4(kPos[gl_VertexID], 0.0, 1.0);
  vUv = 0.5 * (gl_Position.xy + 1.0);
}
)GLSL";

} // namespace

SdfRaymarcher::~SdfRaymarcher() {
  if (vao_ != 0) {
    gl::DeleteVertexArrays(1, &vao_);
    vao_ = 0;
  }
}

core::u64 SdfRaymarcher::structureKey(const SdfGraph& g) {
  // Shader structure key: depends only on topology (op + connections).
  // Parameters and seeds are uploaded via uniforms and do NOT require a recompile.
  core::u64 h = core::fnv1a64("SdfRaymarchShaderV1");
  const int nodeCount = clampNodeCount((int)g.nodes.size());
  h = core::hashCombine(h, (core::u64)nodeCount);
  for (int i = 0; i < nodeCount; ++i) {
    const auto& n = g.nodes[(std::size_t)i];
    const core::u64 packed =
        (core::u64)n.op |
        ((core::u64)(core::u32)(n.in0 + 2) << 8) |
        ((core::u64)(core::u32)(n.in1 + 2) << 20);
    h = core::hashCombine(h, packed);
  }
  return h;
}

std::string SdfRaymarcher::buildFragmentShader(const SdfGraph& g) {
  const int nodeCount = clampNodeCount((int)g.nodes.size());

  auto callNode = [&](int idx, const char* pVar, const char* def) -> std::string {
    if (idx >= 0 && idx < nodeCount) {
      std::ostringstream oss;
      oss << "n" << idx << "(" << pVar << ")";
      return oss.str();
    }
    return def;
  };

  std::ostringstream ss;
  ss << "#version 330 core\n";
  ss << "in vec2 vUv;\n";
  ss << "out vec4 oColor;\n\n";
  ss << "#define MAX_NODES " << kSdfGraphMaxNodes << "\n";
  ss << "#define FAR 1e9\n";
  ss << "#define NEG_FAR -1e9\n\n";

  ss << "uniform vec2  uResolution;\n";
  ss << "uniform float uTime;\n";
  ss << "uniform int   uNodeCount;\n";
  ss << "uniform int   uOutput;\n";
  ss << "uniform float uBounds;\n";
  ss << "uniform float uIso;\n";
  ss << "uniform vec4  uP[MAX_NODES*2];\n";
  ss << "uniform int   uS[MAX_NODES];\n\n";

  ss << "uniform vec3  uCamPos;\n";
  ss << "uniform vec3  uCamFwd;\n";
  ss << "uniform vec3  uCamRight;\n";
  ss << "uniform vec3  uCamUp;\n";
  ss << "uniform float uFovY;\n\n";

  ss << "uniform int   uMaxSteps;\n";
  ss << "uniform float uEps;\n";
  ss << "uniform float uMaxDist;\n";
  ss << "uniform int   uUseBoundsAabb;\n\n";

  ss << "uniform vec3  uLightDir;\n";
  ss << "uniform vec3  uBaseColor;\n";
  ss << "uniform int   uUseAlbedoTex;\n";
  ss << "uniform sampler2D uAlbedoTex;\n";
  ss << "uniform float uTexScale;\n";
  ss << "uniform float uTexRotateDeg;\n";
  ss << "uniform vec2  uTexOffset;\n";
  ss << "uniform vec3  uBgColor;\n";
  ss << "uniform float uAmbient;\n";
  ss << "uniform float uDiffuse;\n";
  ss << "uniform float uSpecular;\n";
  ss << "uniform float uShininess;\n\n";

  ss << "uniform int   uSoftShadow;\n";
  ss << "uniform int   uShadowSteps;\n";
  ss << "uniform float uShadowMaxDist;\n";
  ss << "uniform float uShadowK;\n\n";

  ss << "uniform int   uAO;\n";
  ss << "uniform int   uAoSteps;\n";
  ss << "uniform float uAoStepSize;\n";
  ss << "uniform float uAoStrength;\n\n";

  ss << "uniform int   uDebug;\n\n";

  // --- Hash + value noise + FBm (3D) ---
  ss << R"GLSL(
uint hash_u32(uint x) {
  x ^= x >> 16;
  x *= 0x7feb352du;
  x ^= x >> 15;
  x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}

uint hash3(uvec3 v) {
  uint h = v.x;
  h = hash_u32(h ^ hash_u32(v.y));
  h = hash_u32(h ^ hash_u32(v.z));
  return h;
}

float hash3f(ivec3 p, uint seed) {
  uvec3 v = uvec3(p) + uvec3(seed, seed * 1664525u + 1013904223u, seed * 22695477u + 1u);
  return float(hash3(v)) / 4294967295.0;
}

float fade(float t) {
  return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

float valueNoise3(vec3 p, uint seed) {
  ivec3 i = ivec3(floor(p));
  vec3 f = fract(p);
  vec3 u = vec3(fade(f.x), fade(f.y), fade(f.z));

  float n000 = hash3f(i + ivec3(0,0,0), seed);
  float n100 = hash3f(i + ivec3(1,0,0), seed);
  float n010 = hash3f(i + ivec3(0,1,0), seed);
  float n110 = hash3f(i + ivec3(1,1,0), seed);
  float n001 = hash3f(i + ivec3(0,0,1), seed);
  float n101 = hash3f(i + ivec3(1,0,1), seed);
  float n011 = hash3f(i + ivec3(0,1,1), seed);
  float n111 = hash3f(i + ivec3(1,1,1), seed);

  float nx00 = mix(n000, n100, u.x);
  float nx10 = mix(n010, n110, u.x);
  float nx01 = mix(n001, n101, u.x);
  float nx11 = mix(n011, n111, u.x);
  float nxy0 = mix(nx00, nx10, u.y);
  float nxy1 = mix(nx01, nx11, u.y);
  return mix(nxy0, nxy1, u.z);
}

vec3 grad3(ivec3 p, uint seed) {
  // 12 classic Perlin gradients, normalized (length 1).
  const float INV_SQRT2 = 0.70710678;
  uint h = hash3(uvec3(p) + uvec3(seed, seed * 1664525u + 1013904223u, seed * 22695477u + 1u));
  int idx = int(h % 12u);
  if (idx == 0) return vec3( 1, 1, 0) * INV_SQRT2;
  if (idx == 1) return vec3(-1, 1, 0) * INV_SQRT2;
  if (idx == 2) return vec3( 1,-1, 0) * INV_SQRT2;
  if (idx == 3) return vec3(-1,-1, 0) * INV_SQRT2;
  if (idx == 4) return vec3( 1, 0, 1) * INV_SQRT2;
  if (idx == 5) return vec3(-1, 0, 1) * INV_SQRT2;
  if (idx == 6) return vec3( 1, 0,-1) * INV_SQRT2;
  if (idx == 7) return vec3(-1, 0,-1) * INV_SQRT2;
  if (idx == 8) return vec3( 0, 1, 1) * INV_SQRT2;
  if (idx == 9) return vec3( 0,-1, 1) * INV_SQRT2;
  if (idx == 10) return vec3(0, 1,-1) * INV_SQRT2;
  return vec3(0,-1,-1) * INV_SQRT2;
}

float perlinNoise3(vec3 p, uint seed) {
  ivec3 i = ivec3(floor(p));
  vec3 f = fract(p);
  vec3 u = vec3(fade(f.x), fade(f.y), fade(f.z));

  vec3 g000 = grad3(i + ivec3(0,0,0), seed);
  vec3 g100 = grad3(i + ivec3(1,0,0), seed);
  vec3 g010 = grad3(i + ivec3(0,1,0), seed);
  vec3 g110 = grad3(i + ivec3(1,1,0), seed);
  vec3 g001 = grad3(i + ivec3(0,0,1), seed);
  vec3 g101 = grad3(i + ivec3(1,0,1), seed);
  vec3 g011 = grad3(i + ivec3(0,1,1), seed);
  vec3 g111 = grad3(i + ivec3(1,1,1), seed);

  float d000 = dot(g000, f - vec3(0,0,0));
  float d100 = dot(g100, f - vec3(1,0,0));
  float d010 = dot(g010, f - vec3(0,1,0));
  float d110 = dot(g110, f - vec3(1,1,0));
  float d001 = dot(g001, f - vec3(0,0,1));
  float d101 = dot(g101, f - vec3(1,0,1));
  float d011 = dot(g011, f - vec3(0,1,1));
  float d111 = dot(g111, f - vec3(1,1,1));

  float nx00 = mix(d000, d100, u.x);
  float nx10 = mix(d010, d110, u.x);
  float nx01 = mix(d001, d101, u.x);
  float nx11 = mix(d011, d111, u.x);
  float nxy0 = mix(nx00, nx10, u.y);
  float nxy1 = mix(nx01, nx11, u.y);
  float n = mix(nxy0, nxy1, u.z);

  // Typical normalization for 3D Perlin; keeps output mostly in [-1,1].
  n *= 0.577350269;
  return clamp(0.5 + 0.5 * n, 0.0, 1.0);
}

float fbmAmpSum(int octaves, float gain) {
  float amp = 0.5;
  float sum = 0.0;
  for (int i = 0; i < 12; ++i) {
    if (i >= octaves) break;
    sum += amp;
    amp *= gain;
  }
  return sum;
}

float signedFbm01(uint seed, vec3 p, int octaves, float lacunarity, float gain) {
  octaves = clamp(octaves, 1, 12);
  if (lacunarity <= 0.0) lacunarity = 2.0;
  if (gain <= 0.0) gain = 0.5;

  float sum = 0.0;
  float amp = 0.5;
  float freq = 1.0;
  for (int i = 0; i < 12; ++i) {
    if (i >= octaves) break;
    float n = valueNoise3(p * freq, seed + uint(i) * 101u);
    sum += n * amp;
    freq *= lacunarity;
    amp *= gain;
  }

  float aSum = fbmAmpSum(octaves, gain);
  float n01 = (aSum > 0.0) ? (sum / aSum) : 0.0;
  n01 = clamp(n01, 0.0, 1.0);
  return n01 * 2.0 - 1.0;
}

float signedFbmPerlin01(uint seed, vec3 p, int octaves, float lacunarity, float gain) {
  octaves = clamp(octaves, 1, 12);
  if (lacunarity <= 0.0) lacunarity = 2.0;
  if (gain <= 0.0) gain = 0.5;

  float sum = 0.0;
  float amp = 0.5;
  float freq = 1.0;
  for (int i = 0; i < 12; ++i) {
    if (i >= octaves) break;
    float n = perlinNoise3(p * freq, seed + uint(i) * 101u);
    sum += n * amp;
    freq *= lacunarity;
    amp *= gain;
  }

  float aSum = fbmAmpSum(octaves, gain);
  float n01 = (aSum > 0.0) ? (sum / aSum) : 0.0;
  n01 = clamp(n01, 0.0, 1.0);
  return n01 * 2.0 - 1.0;
}
)GLSL";

  // --- SDF primitives + ops ---
  ss << R"GLSL(
float sdfSphere(vec3 p, float r) {
  return length(p) - r;
}

float sdfBox(vec3 p, vec3 b) {
  vec3 q = abs(p) - b;
  float outside = length(max(q, vec3(0.0)));
  float inside = min(max(q.x, max(q.y, q.z)), 0.0);
  return outside + inside;
}

float sdfCapsule(vec3 p, vec3 a, vec3 b, float r) {
  vec3 pa = p - a;
  vec3 ba = b - a;
  float baba = dot(ba, ba);
  float h = (baba > 1e-12) ? clamp(dot(pa, ba) / baba, 0.0, 1.0) : 0.0;
  vec3 d = pa - ba * h;
  return length(d) - r;
}

float sdfTorusY(vec3 p, float majorR, float minorR) {
  vec2 q = vec2(length(p.xz) - majorR, p.y);
  return length(q) - minorR;
}

float sdfPlane(vec3 p, vec3 n, float d) {
  float lenN = length(n);
  n = (lenN > 1e-6) ? (n / lenN) : vec3(0.0, 1.0, 0.0);
  return dot(p, n) + d;
}

float opSmoothUnion(float d1, float d2, float k) {
  k = max(k, 1e-6);
  float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
  float m = mix(d2, d1, h);
  return m - k * h * (1.0 - h);
}
)GLSL";

  // --- Domain helpers (transforms / repetition) ---
  ss << R"GLSL(
vec3 rotX(vec3 p, float a) {
  float c = cos(a), s = sin(a);
  return vec3(p.x, c * p.y - s * p.z, s * p.y + c * p.z);
}

vec3 rotY(vec3 p, float a) {
  float c = cos(a), s = sin(a);
  return vec3(c * p.x + s * p.z, p.y, -s * p.x + c * p.z);
}

vec3 rotZ(vec3 p, float a) {
  float c = cos(a), s = sin(a);
  return vec3(c * p.x - s * p.y, s * p.x + c * p.y, p.z);
}

vec3 repeatP(vec3 p, vec3 period, vec3 off) {
  period = max(period, vec3(1e-6));
  p += off;
  p = mod(p + 0.5 * period, period) - 0.5 * period;
  return p;
}

vec3 mirrorP(vec3 p, vec3 mask, vec3 pivot) {
  vec3 q = p - pivot;
  if (mask.x > 0.5) q.x = abs(q.x);
  if (mask.y > 0.5) q.y = abs(q.y);
  if (mask.z > 0.5) q.z = abs(q.z);
  return q + pivot;
}
)GLSL";

  // --- Per-node SDF functions ---
  if (nodeCount > 0) {
    // Prototypes first (avoid ordering issues if loaded graphs reference forward nodes).
    for (int i = 0; i < nodeCount; ++i) {
      ss << "float n" << i << "(vec3 p);\n";
    }
    ss << "\n";

    for (int i = 0; i < nodeCount; ++i) {
      const SdfNode& node = g.nodes[(std::size_t)i];
      const int p0 = i * 2;
      const int p1 = i * 2 + 1;

      ss << "float n" << i << "(vec3 p) {\n";
      ss << "  vec4 a = uP[" << p0 << "];\n";
      ss << "  vec4 b = uP[" << p1 << "];\n";

      switch (node.op) {
        case SdfNodeOp::Constant: {
          ss << "  return a.x;\n";
        } break;

        case SdfNodeOp::Sphere: {
          ss << "  vec3 c = a.xyz;\n";
          ss << "  float r = max(0.0, a.w);\n";
          ss << "  return sdfSphere(p - c, r);\n";
        } break;

        case SdfNodeOp::Box: {
          ss << "  vec3 c = a.xyz;\n";
          ss << "  vec3 halfSize = vec3(max(0.0, a.w), max(0.0, b.x), max(0.0, b.y));\n";
          ss << "  return sdfBox(p - c, halfSize);\n";
        } break;

        case SdfNodeOp::Capsule: {
          ss << "  vec3 A = a.xyz;\n";
          ss << "  vec3 B = vec3(a.w, b.x, b.y);\n";
          ss << "  float r = max(0.0, b.z);\n";
          ss << "  return sdfCapsule(p, A, B, r);\n";
        } break;

        case SdfNodeOp::TorusY: {
          ss << "  vec3 c = a.xyz;\n";
          ss << "  float R = max(0.0, a.w);\n";
          ss << "  float r = max(0.0, b.x);\n";
          ss << "  return sdfTorusY(p - c, R, r);\n";
        } break;

        case SdfNodeOp::Plane: {
          ss << "  vec3 nrm = a.xyz;\n";
          ss << "  float d = a.w;\n";
          ss << "  return sdfPlane(p, nrm, d);\n";
        } break;

        case SdfNodeOp::Union: {
          ss << "  float da = " << callNode(node.in0, "p", "FAR") << ";\n";
          ss << "  float db = " << callNode(node.in1, "p", "FAR") << ";\n";
          ss << "  return min(da, db);\n";
        } break;

        case SdfNodeOp::SmoothUnion: {
          ss << "  float da = " << callNode(node.in0, "p", "FAR") << ";\n";
          ss << "  float db = " << callNode(node.in1, "p", "FAR") << ";\n";
          ss << "  return opSmoothUnion(da, db, a.x);\n";
        } break;

        case SdfNodeOp::Intersect: {
          ss << "  float da = " << callNode(node.in0, "p", "NEG_FAR") << ";\n";
          ss << "  float db = " << callNode(node.in1, "p", "NEG_FAR") << ";\n";
          ss << "  return max(da, db);\n";
        } break;

        case SdfNodeOp::Subtract: {
          ss << "  float da = " << callNode(node.in0, "p", "FAR") << ";\n";
          ss << "  float db = " << callNode(node.in1, "p", "FAR") << ";\n";
          ss << "  return max(da, -db);\n";
        } break;

        case SdfNodeOp::NoiseDisplace: {
          ss << "  float d = " << callNode(node.in0, "p", "FAR") << ";\n";
          ss << "  uint seed = uint(uS[" << i << "]);\n";
          ss << "  float freq = max(0.0, a.x);\n";
          ss << "  float amp = a.y;\n";
          ss << "  int oct = int(clamp(floor(a.z + 0.5), 1.0, 12.0));\n";
          ss << "  float lac = (a.w <= 0.0) ? 2.0 : a.w;\n";
          ss << "  float gain = (b.x <= 0.0) ? 0.5 : b.x;\n";
          ss << "  vec3 off = vec3(b.y, b.z, b.w);\n";
          ss << "  float n = signedFbm01(seed, p * freq + off, oct, lac, gain);\n";
          ss << "  return d + amp * n;\n";
        } break;

        case SdfNodeOp::NoiseDisplacePerlin: {
          ss << "  float d = " << callNode(node.in0, "p", "FAR") << ";\n";
          ss << "  uint seed = uint(uS[" << i << "]);\n";
          ss << "  float freq = max(0.0, a.x);\n";
          ss << "  float amp = a.y;\n";
          ss << "  int oct = int(clamp(floor(a.z + 0.5), 1.0, 12.0));\n";
          ss << "  float lac = (a.w <= 0.0) ? 2.0 : a.w;\n";
          ss << "  float gain = (b.x <= 0.0) ? 0.5 : b.x;\n";
          ss << "  vec3 off = vec3(b.y, b.z, b.w);\n";
          ss << "  float n = signedFbmPerlin01(seed, p * freq + off, oct, lac, gain);\n";
          ss << "  return d + amp * n;\n";
        } break;

        case SdfNodeOp::Shell: {
          ss << "  float d = " << callNode(node.in0, "p", "FAR") << ";\n";
          ss << "  float t = max(0.0, a.x);\n";
          ss << "  return abs(d) - t;\n";
        } break;

        // --- Space transforms / domain operations ---
        case SdfNodeOp::Translate: {
          ss << "  vec3 t = a.xyz;\n";
          ss << "  vec3 q = p - t;\n";
          ss << "  return " << callNode(node.in0, "q", "FAR") << ";\n";
        } break;

        case SdfNodeOp::RotateX: {
          ss << "  float ang = radians(a.x);\n";
          ss << "  vec3 pivot = vec3(a.y, a.z, a.w);\n";
          ss << "  vec3 q = p - pivot;\n";
          ss << "  q = rotX(q, -ang);\n";
          ss << "  q += pivot;\n";
          ss << "  return " << callNode(node.in0, "q", "FAR") << ";\n";
        } break;

        case SdfNodeOp::RotateY: {
          ss << "  float ang = radians(a.x);\n";
          ss << "  vec3 pivot = vec3(a.y, a.z, a.w);\n";
          ss << "  vec3 q = p - pivot;\n";
          ss << "  q = rotY(q, -ang);\n";
          ss << "  q += pivot;\n";
          ss << "  return " << callNode(node.in0, "q", "FAR") << ";\n";
        } break;

        case SdfNodeOp::RotateZ: {
          ss << "  float ang = radians(a.x);\n";
          ss << "  vec3 pivot = vec3(a.y, a.z, a.w);\n";
          ss << "  vec3 q = p - pivot;\n";
          ss << "  q = rotZ(q, -ang);\n";
          ss << "  q += pivot;\n";
          ss << "  return " << callNode(node.in0, "q", "FAR") << ";\n";
        } break;

        case SdfNodeOp::Scale: {
          ss << "  float sc = max(1e-6, abs(a.x));\n";
          ss << "  vec3 pivot = vec3(a.y, a.z, a.w);\n";
          ss << "  vec3 q = (p - pivot) / sc + pivot;\n";
          ss << "  float d = " << callNode(node.in0, "q", "FAR") << ";\n";
          ss << "  return d * sc;\n";
        } break;

        case SdfNodeOp::Repeat: {
          ss << "  vec3 period = max(vec3(a.x, a.y, a.z), vec3(1e-6));\n";
          ss << "  vec3 off = vec3(a.w, b.x, b.y);\n";
          ss << "  vec3 q = repeatP(p, period, off);\n";
          ss << "  return " << callNode(node.in0, "q", "FAR") << ";\n";
        } break;

        case SdfNodeOp::Mirror: {
          ss << "  vec3 mask = vec3(a.x, a.y, a.z);\n";
          ss << "  vec3 pivot = vec3(a.w, b.x, b.y);\n";
          ss << "  vec3 q = mirrorP(p, mask, pivot);\n";
          ss << "  return " << callNode(node.in0, "q", "FAR") << ";\n";
        } break;

        case SdfNodeOp::TwistY: {
          ss << "  float k = radians(a.x);\n";
          ss << "  vec3 pivot = vec3(a.y, a.z, a.w);\n";
          ss << "  vec3 q = p - pivot;\n";
          ss << "  float ang = k * q.y;\n";
          ss << "  q = rotY(q, -ang);\n";
          ss << "  q += pivot;\n";
          ss << "  return " << callNode(node.in0, "q", "FAR") << ";\n";
        } break;

        default: {
          ss << "  return FAR;\n";
        } break;
      }

      ss << "}\n\n";
    }
  }

  // --- Graph evaluation (dynamic output selection) ---
  ss << "\nfloat evalGraph(vec3 p) {\n";
  ss << "  int n = min(uNodeCount, MAX_NODES);\n";
  ss << "  if (n <= 0) return FAR;\n";
  ss << "  int oi = clamp(uOutput, 0, n - 1);\n";

  if (nodeCount <= 0) {
    ss << "  return FAR;\n";
  } else {
    for (int i = 0; i < nodeCount; ++i) {
      ss << "  if (oi == " << i << ") return n" << i << "(p);\n";
    }
    ss << "  return n" << (nodeCount - 1) << "(p);\n";
  }

  ss << "}\n\n";

  // Field function includes iso.
  ss << R"GLSL(
float field(vec3 p) {
  return evalGraph(p) - uIso;
}

bool inBounds(vec3 p) {
  float b = uBounds;
  return all(greaterThanEqual(p, vec3(-b))) && all(lessThanEqual(p, vec3(b)));
}

bool rayAabb(vec3 ro, vec3 rd, vec3 bmin, vec3 bmax, out float tmin, out float tmax) {
  vec3 inv = 1.0 / rd;
  vec3 t0 = (bmin - ro) * inv;
  vec3 t1 = (bmax - ro) * inv;
  vec3 tsm = min(t0, t1);
  vec3 tsM = max(t0, t1);
  tmin = max(max(tsm.x, tsm.y), tsm.z);
  tmax = min(min(tsM.x, tsM.y), tsM.z);
  return tmax >= max(tmin, 0.0);
}

vec3 calcNormal(vec3 p) {
  float e = uEps;
  vec2 h = vec2(e, 0);
  float dx = field(p + vec3(h.x, h.y, h.y)) - field(p - vec3(h.x, h.y, h.y));
  float dy = field(p + vec3(h.y, h.x, h.y)) - field(p - vec3(h.y, h.x, h.y));
  float dz = field(p + vec3(h.y, h.y, h.x)) - field(p - vec3(h.y, h.y, h.x));
  return normalize(vec3(dx, dy, dz));
}

float softShadow(vec3 ro, vec3 rd, float tMin, float tMax, float k) {
  float res = 1.0;
  float t = tMin;
  for (int i = 0; i < 64; ++i) {
    if (i >= uShadowSteps) break;
    float h = field(ro + rd * t);
    if (h < uEps) return 0.0;
    res = min(res, k * h / t);
    t += clamp(h, 0.01, 0.25);
    if (t > tMax) break;
  }
  return clamp(res, 0.0, 1.0);
}

float ambientOcclusion(vec3 p, vec3 n) {
  float occ = 0.0;
  float sca = 1.0;
  for (int i = 0; i < 32; ++i) {
    if (i >= uAoSteps) break;
    float h = uAoStepSize * float(i + 1);
    float d = field(p + n * h);
    occ += (h - d) * sca;
    sca *= 0.95;
  }
  return clamp(1.0 - uAoStrength * occ, 0.0, 1.0);
}

vec2 texUvFromPoint(vec3 p, vec3 n) {
  // Simple triplanar-like: choose projection based on normal.
  vec3 an = abs(n);
  vec2 uv;
  if (an.x > an.y && an.x > an.z) {
    uv = p.yz;
  } else if (an.y > an.z) {
    uv = p.xz;
  } else {
    uv = p.xy;
  }

  // Scale + rotate + offset.
  uv *= uTexScale;
  float a = radians(uTexRotateDeg);
  float c = cos(a), s = sin(a);
  uv = vec2(c * uv.x - s * uv.y, s * uv.x + c * uv.y);
  uv += uTexOffset;
  return uv;
}

vec3 sampleAlbedo(vec3 p, vec3 n) {
  if (uUseAlbedoTex == 0) return uBaseColor;
  vec2 uv = texUvFromPoint(p, n);
  return texture(uAlbedoTex, uv).rgb;
}

void main() {
  // Normalized pixel coordinates (-1..1) with aspect.
  vec2 uv = (gl_FragCoord.xy / uResolution.xy) * 2.0 - 1.0;
  uv.x *= uResolution.x / uResolution.y;

  // Camera ray.
  float tanHalf = tan(0.5 * uFovY);
  vec3 rd = normalize(uCamFwd + uv.x * uCamRight * tanHalf + uv.y * uCamUp * tanHalf);
  vec3 ro = uCamPos;

  float t = 0.0;
  float tEnd = uMaxDist;

  if (uUseBoundsAabb != 0) {
    float tmin, tmax;
    float b = uBounds;
    if (!rayAabb(ro, rd, vec3(-b), vec3(b), tmin, tmax)) {
      oColor = vec4(uBgColor, 1.0);
      return;
    }
    t = max(tmin, 0.0);
    tEnd = min(tmax, uMaxDist);
  }

  // Raymarch.
  int steps = 0;
  float d = 0.0;
  for (int i = 0; i < 2048; ++i) {
    if (i >= uMaxSteps) break;
    vec3 p = ro + rd * t;
    d = field(p);
    steps = i;
    if (abs(d) < uEps) break;
    t += d;
    if (t > tEnd) break;
  }

  if (t > tEnd || abs(d) >= uEps) {
    oColor = vec4(uBgColor, 1.0);
    return;
  }

  vec3 p = ro + rd * t;
  vec3 n = calcNormal(p);

  // Debug outputs.
  if (uDebug == 1) {
    oColor = vec4(n * 0.5 + 0.5, 1.0);
    return;
  } else if (uDebug == 2) {
    float s = float(steps) / float(max(uMaxSteps, 1));
    oColor = vec4(vec3(s), 1.0);
    return;
  } else if (uDebug == 3) {
    float dd = clamp((d + 1.0) * 0.5, 0.0, 1.0);
    oColor = vec4(vec3(dd), 1.0);
    return;
  }

  vec3 ld = normalize(uLightDir);
  float ndl = max(dot(n, ld), 0.0);

  vec3 albedo = sampleAlbedo(p, n);

  vec3 col = vec3(0.0);
  col += albedo * uAmbient;
  col += albedo * ndl * uDiffuse;

  // Specular.
  vec3 h = normalize(ld - rd);
  float ndh = max(dot(n, h), 0.0);
  col += pow(ndh, uShininess) * uSpecular;

  // Shadows.
  if (uSoftShadow != 0 && uShadowSteps > 0) {
    float sh = softShadow(p + n * (uEps * 2.0), ld, 0.02, uShadowMaxDist, uShadowK);
    col *= mix(1.0, sh, 0.9);
  }

  // AO.
  if (uAO != 0 && uAoSteps > 0) {
    float ao = ambientOcclusion(p, n);
    col *= ao;
  }

  // Fog based on distance.
  float fog = clamp(t / uMaxDist, 0.0, 1.0);
  col = mix(col, uBgColor, fog);

  oColor = vec4(col, 1.0);
}
)GLSL";

  return ss.str();
}


bool SdfRaymarcher::ensureShader(const SdfGraph& g, std::string* outError) {
  const core::u64 key = structureKey(g);
  stats_.shaderRebuilt = false;
  stats_.shaderBuildMs = 0.0;

  if (shader_.handle() != 0 && key == shaderKey_) return true;

  lastFragSrc_ = buildFragmentShader(g);

  auto t0 = std::chrono::high_resolution_clock::now();

  ShaderProgram newProg;
  std::string err;
  if (!newProg.build(kFullscreenTriVS, lastFragSrc_, &err)) {
    if (outError) *outError = err;
    return false;
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  stats_.shaderBuildMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
  stats_.shaderRebuilt = true;

  shader_ = std::move(newProg);
  shaderKey_ = key;

  if (vao_ == 0) {
    gl::GenVertexArrays(1, &vao_);
  }

  return true;
}

bool SdfRaymarcher::render(const SdfGraph& g,
                           int width,
                           int height,
                           float timeSec,
                           float bounds,
                           float iso,
                           const SdfRaymarchCamera& cam,
                           const SdfRaymarchSettings& s,
                           std::string* outError,
                           const SdfRaymarchMaterial* mat) {
  width = std::max(1, width);
  height = std::max(1, height);

  if (!target_.isInited() || target_.width() != width || target_.height() != height) {
    if (!target_.init(width, height)) {
      if (outError) *outError = "SdfRaymarcher: failed to init RenderTarget2D";
      return false;
    }
  }

  if (!ensureShader(g, outError)) return false;

  const int nodeCount = clampNodeCount((int)g.nodes.size());

  packedParams_.fill(0.0f);
  packedSeeds_.fill(0);

  for (int i = 0; i < nodeCount; ++i) {
    const auto& n = g.nodes[(std::size_t)i];
    const int base = i * 8;
    packedParams_[base + 0] = n.p0;
    packedParams_[base + 1] = n.p1;
    packedParams_[base + 2] = n.p2;
    packedParams_[base + 3] = n.p3;
    packedParams_[base + 4] = n.p4;
    packedParams_[base + 5] = n.p5;
    packedParams_[base + 6] = n.p6;
    packedParams_[base + 7] = n.p7;

    // Deterministic per-node seed: mix graph seed + node seed + index.
    const core::u64 seed64 = core::hashCombine(g.seed, core::hashCombine((core::u64)i, n.seed));
    packedSeeds_[i] = (int)(seed64 & 0x7fffffff);
  }

  // --- Save GL state we touch (important: we render during ImGui UI rendering). ---
  GLint prevFbo = 0;
  GLint prevProg = 0;
  GLint prevVao = 0;
  GLint prevViewport[4] = {0, 0, 0, 0};

  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
  glGetIntegerv(GL_CURRENT_PROGRAM, &prevProg);
  glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &prevVao);
  glGetIntegerv(GL_VIEWPORT, prevViewport);

  // Texture state (we may bind an albedo texture for the preview).
  GLint prevActiveTex = 0;
  GLint prevTex0 = 0;
  glGetIntegerv(GL_ACTIVE_TEXTURE, &prevActiveTex);
  gl::ActiveTexture(GL_TEXTURE0);
  glGetIntegerv(GL_TEXTURE_BINDING_2D, &prevTex0);

  const GLboolean depthEnabled = glIsEnabled(GL_DEPTH_TEST);
  const GLboolean blendEnabled = glIsEnabled(GL_BLEND);
  const GLboolean cullEnabled = glIsEnabled(GL_CULL_FACE);

  auto t0 = std::chrono::high_resolution_clock::now();

  target_.begin();
  glViewport(0, 0, width, height);

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  glDisable(GL_CULL_FACE);

  glClearColor(0, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT);

  shader_.bind();
  shader_.setUniform2f("uResolution", (float)width, (float)height);
  shader_.setUniform1f("uTime", timeSec);
  shader_.setUniform1i("uNodeCount", nodeCount);
  shader_.setUniform1i("uOutput", g.output);
  shader_.setUniform1f("uBounds", bounds);
  shader_.setUniform1f("uIso", iso);
  shader_.setUniform1i("uMaxSteps", std::max(1, s.maxSteps));
  shader_.setUniform1f("uEps", std::max(1e-6f, s.epsilon));
  shader_.setUniform1f("uMaxDist", std::max(0.1f, s.maxDistance));
  shader_.setUniform1i("uUseBoundsAabb", s.useBoundsAabb ? 1 : 0);

  shader_.setUniform3f("uCamPos", cam.pos[0], cam.pos[1], cam.pos[2]);
  shader_.setUniform3f("uCamFwd", cam.forward[0], cam.forward[1], cam.forward[2]);
  shader_.setUniform3f("uCamRight", cam.right[0], cam.right[1], cam.right[2]);
  shader_.setUniform3f("uCamUp", cam.up[0], cam.up[1], cam.up[2]);
  shader_.setUniform1f("uFovY", std::max(0.05f, cam.fovYRadians));

  // Normalize light dir on CPU side a bit (shader will normalize again).
  shader_.setUniform3f("uLightDir", s.lightDir[0], s.lightDir[1], s.lightDir[2]);
  shader_.setUniform3f("uBaseColor", s.baseColor[0], s.baseColor[1], s.baseColor[2]);
  const bool useAlbedo = (mat && mat->albedoTex && mat->albedoTex->handle() != 0);
  shader_.setUniform1i("uUseAlbedoTex", useAlbedo ? 1 : 0);
  shader_.setUniform1i("uAlbedoTex", 0);
  shader_.setUniform1f("uTexScale", useAlbedo ? mat->uvScale : 1.0f);
  shader_.setUniform1f("uTexRotateDeg", useAlbedo ? mat->uvRotateDeg : 0.0f);
  shader_.setUniform2f("uTexOffset",
                       useAlbedo ? mat->uvOffset[0] : 0.0f,
                       useAlbedo ? mat->uvOffset[1] : 0.0f);
  if (useAlbedo) {
    mat->albedoTex->bind(0);
  } else {
    gl::ActiveTexture(GL_TEXTURE0);
    gl::BindTexture(GL_TEXTURE_2D, 0);
  }
  shader_.setUniform3f("uBgColor", s.backgroundColor[0], s.backgroundColor[1], s.backgroundColor[2]);
  shader_.setUniform1f("uAmbient", s.ambient);
  shader_.setUniform1f("uDiffuse", s.diffuse);
  shader_.setUniform1f("uSpecular", s.specular);
  shader_.setUniform1f("uShininess", s.shininess);

  shader_.setUniform1i("uSoftShadow", s.softShadows ? 1 : 0);
  shader_.setUniform1i("uShadowSteps", std::max(0, s.shadowSteps));
  shader_.setUniform1f("uShadowMaxDist", std::max(0.0f, s.shadowMaxDistance));
  shader_.setUniform1f("uShadowK", std::max(0.0f, s.shadowK));

  shader_.setUniform1i("uAO", s.ambientOcclusion ? 1 : 0);
  shader_.setUniform1i("uAoSteps", std::max(0, s.aoSteps));
  shader_.setUniform1f("uAoStepSize", std::max(0.0f, s.aoStepSize));
  shader_.setUniform1f("uAoStrength", std::max(0.0f, s.aoStrength));

  shader_.setUniform1i("uDebug", (int)s.debug);

  // Upload arrays.
  shader_.setUniform4fv("uP[0]", kSdfGraphMaxNodes * 2, packedParams_.data());
  shader_.setUniform1iv("uS[0]", kSdfGraphMaxNodes, packedSeeds_.data());

  gl::BindVertexArray(vao_);
  glDrawArrays(GL_TRIANGLES, 0, 3);

  target_.end();

  auto t1 = std::chrono::high_resolution_clock::now();
  stats_.drawMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

  // --- Restore GL state ---
  // Restore texture bindings.
  gl::ActiveTexture(GL_TEXTURE0);
  gl::BindTexture(GL_TEXTURE_2D, (GLuint)prevTex0);
  gl::ActiveTexture((GLenum)prevActiveTex);

  gl::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
  glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);

  gl::UseProgram((GLuint)prevProg);
  gl::BindVertexArray((GLuint)prevVao);

  if (depthEnabled) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
  if (blendEnabled) glEnable(GL_BLEND); else glDisable(GL_BLEND);
  if (cullEnabled) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);

  return true;
}

} // namespace stellar::render
