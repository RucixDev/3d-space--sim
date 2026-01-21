#include "stellar/render/ProceduralGraph.h"

#include "stellar/core/Hash.h"
#include "stellar/render/Gl.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace stellar::render {

namespace {

int clampPalCount(int n) {
  return std::clamp(n, 2, kProcGraphMaxPaletteStops);
}

int clampNodeCount(int n) {
  return std::clamp(n, 0, kProcGraphMaxNodes);
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
  vec2 p = kPos[gl_VertexID];
  gl_Position = vec4(p, 0.0, 1.0);
  vUv = p * 0.5 + 0.5;
}
)GLSL";

} // namespace

const char* procNodeOpName(ProcNodeOp op) {
  switch (op) {
    case ProcNodeOp::Constant: return "Constant";
    case ProcNodeOp::UvX: return "UvX";
    case ProcNodeOp::UvY: return "UvY";
    case ProcNodeOp::Add: return "Add";
    case ProcNodeOp::Sub: return "Sub";
    case ProcNodeOp::Mul: return "Mul";
    case ProcNodeOp::Div: return "Div";
    case ProcNodeOp::Min: return "Min";
    case ProcNodeOp::Max: return "Max";
    case ProcNodeOp::Abs: return "Abs";
    case ProcNodeOp::Invert: return "Invert";
    case ProcNodeOp::Fract: return "Fract";
    case ProcNodeOp::Clamp01: return "Clamp01";
    case ProcNodeOp::Smoothstep: return "Smoothstep";
    case ProcNodeOp::Pow: return "Pow";
    case ProcNodeOp::Sine: return "Sine";
    case ProcNodeOp::Noise2D: return "Noise2D";
    case ProcNodeOp::Fbm2D: return "FBm2D";
    case ProcNodeOp::Perlin2D: return "Perlin2D";
    case ProcNodeOp::FbmPerlin2D: return "FBmPerlin2D";
    case ProcNodeOp::RidgedFbmPerlin2D: return "RidgedFBmPerlin2D";
    case ProcNodeOp::Voronoi2D: return "Voronoi2D";
    case ProcNodeOp::Warp: return "Warp";
    case ProcNodeOp::Pan: return "Pan";
    default: return "Unknown";
  }
}

static std::string lowerAscii(std::string s) {
  for (char& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

std::optional<ProcNodeOp> procNodeOpFromString(const std::string& s) {
  const std::string k = lowerAscii(s);
  if (k == "constant") return ProcNodeOp::Constant;
  if (k == "uvx" || k == "u") return ProcNodeOp::UvX;
  if (k == "uvy" || k == "v") return ProcNodeOp::UvY;

  if (k == "add" || k == "+") return ProcNodeOp::Add;
  if (k == "sub" || k == "subtract" || k == "-") return ProcNodeOp::Sub;
  if (k == "mul" || k == "multiply" || k == "*") return ProcNodeOp::Mul;
  if (k == "div" || k == "divide" || k == "/") return ProcNodeOp::Div;
  if (k == "min") return ProcNodeOp::Min;
  if (k == "max") return ProcNodeOp::Max;
  if (k == "abs") return ProcNodeOp::Abs;
  if (k == "invert" || k == "inv") return ProcNodeOp::Invert;
  if (k == "fract" || k == "frac") return ProcNodeOp::Fract;
  if (k == "clamp01" || k == "clamp" || k == "saturate") return ProcNodeOp::Clamp01;
  if (k == "smoothstep" || k == "smooth") return ProcNodeOp::Smoothstep;
  if (k == "pow" || k == "power") return ProcNodeOp::Pow;
  if (k == "sine" || k == "sin") return ProcNodeOp::Sine;

  if (k == "noise2d" || k == "noise") return ProcNodeOp::Noise2D;
  if (k == "fbm2d" || k == "fbm" || k == "fbm2") return ProcNodeOp::Fbm2D;
  if (k == "perlin2d" || k == "perlin" || k == "gradnoise2d" || k == "gradient2d") return ProcNodeOp::Perlin2D;
  if (k == "fbmperlin2d" || k == "fbmperlin" || k == "perlinfbm2d") return ProcNodeOp::FbmPerlin2D;
  if (k == "ridgedfbmperlin2d" || k == "ridgedfbmperlin" || k == "ridgedperlin2d") return ProcNodeOp::RidgedFbmPerlin2D;
  if (k == "voronoi2d" || k == "voronoi") return ProcNodeOp::Voronoi2D;
  if (k == "warp") return ProcNodeOp::Warp;
  if (k == "pan") return ProcNodeOp::Pan;
  return std::nullopt;
}

const char* procGraphPresetName(ProcGraphPreset preset) {
  switch (preset) {
    case ProcGraphPreset::Nebula: return "Nebula";
    case ProcGraphPreset::Marble: return "Marble";
    case ProcGraphPreset::Lava: return "Lava";
    case ProcGraphPreset::AlienCircuit: return "Alien Circuit";
    case ProcGraphPreset::Rocky: return "Rocky";
    default: return "Unknown";
  }
}

ProcGraph ProcGraph::makeDefault() {
  return makeProceduralGraphPreset(ProcGraphPreset::Nebula, 0xC0FFEEULL);
}

ProcGraph makeProceduralGraphPreset(ProcGraphPreset preset, core::u64 seed) {
  ProcGraph g{};
  g.seed = seed;
  g.nodes.clear();

  auto pal = [&](std::initializer_list<ProcPaletteStop> stops) {
    g.paletteCount = std::min<int>((int)stops.size(), kProcGraphMaxPaletteStops);
    int i = 0;
    for (auto& s : stops) {
      if (i >= g.paletteCount) break;
      g.palette[i++] = s;
    }
    // Fill remainder for safety.
    for (; i < kProcGraphMaxPaletteStops; ++i) g.palette[i] = g.palette[g.paletteCount - 1];
  };

  auto add = [&](ProcNode n) -> int {
    g.nodes.push_back(n);
    return (int)g.nodes.size() - 1;
  };

  switch (preset) {
    case ProcGraphPreset::Nebula: {
      // FBm -> Warp(FBm) -> Pow -> Clamp01
      // Palette: deep blue -> purple -> pink -> pale
      // Use Perlin FBm to reduce grid artifacts in large smooth gradients.
      int n0 = add(ProcNode{ProcNodeOp::FbmPerlin2D, -1, -1, -1, 3.2f, 5.0f, 2.0f, 0.52f, 0});
      int n1 = add(ProcNode{ProcNodeOp::Warp, n0, -1, -1, 0.32f, 2.4f, 4.0f, 0.55f, 0});
      int n2 = add(ProcNode{ProcNodeOp::Pow, n1, -1, -1, 1.6f, 0, 0, 0, 0});
      int n3 = add(ProcNode{ProcNodeOp::Clamp01, n2, -1, -1, 0, 0, 0, 0, 0});
      g.output = n3;

      g.usePalette = true;
      pal({
          {0.00f, 0.02f, 0.03f, 0.06f},
          {0.35f, 0.10f, 0.05f, 0.22f},
          {0.65f, 0.60f, 0.10f, 0.55f},
          {1.00f, 0.90f, 0.82f, 1.00f},
      });
    } break;

    case ProcGraphPreset::Marble: {
      // sin( uv.x * freq + fbm*amp )
      int uvx = add(ProcNode{ProcNodeOp::UvX, -1, -1, -1, 0, 0, 0, 0, 0});
      int fbm = add(ProcNode{ProcNodeOp::Fbm2D, -1, -1, -1, 6.0f, 5.0f, 2.0f, 0.50f, 0});
      int amp = add(ProcNode{ProcNodeOp::Constant, -1, -1, -1, 0.35f, 0, 0, 0, 0});
      int fbmScaled = add(ProcNode{ProcNodeOp::Mul, fbm, amp, -1, 0, 0, 0, 0, 0});
      int base = add(ProcNode{ProcNodeOp::Add, uvx, fbmScaled, -1, 0, 0, 0, 0, 0});
      int s = add(ProcNode{ProcNodeOp::Sine, base, -1, -1, 45.0f, 0.0f, 0, 0, 0});
      g.output = s;

      g.usePalette = true;
      pal({
          {0.00f, 0.05f, 0.05f, 0.06f},
          {0.45f, 0.55f, 0.55f, 0.58f},
          {0.75f, 0.85f, 0.83f, 0.80f},
          {1.00f, 1.00f, 0.98f, 0.95f},
      });
    } break;

    case ProcGraphPreset::Lava: {
      // Base fbm + crack mask from voronoi edges
      int base = add(ProcNode{ProcNodeOp::Fbm2D, -1, -1, -1, 4.0f, 6.0f, 2.1f, 0.55f, 0});
      int contrast = add(ProcNode{ProcNodeOp::Pow, base, -1, -1, 1.8f, 0, 0, 0, 0});
      int vor = add(ProcNode{ProcNodeOp::Voronoi2D, -1, -1, -1, 10.0f, 1.00f, 0, 0, 0});
      int cracks = add(ProcNode{ProcNodeOp::Invert, vor, -1, -1, 0, 0, 0, 0, 0});
      int crackLines = add(ProcNode{ProcNodeOp::Smoothstep, cracks, -1, -1, 0.55f, 0.75f, 0, 0, 0});
      // Blend in cracks via max()
      int out = add(ProcNode{ProcNodeOp::Max, contrast, crackLines, -1, 0, 0, 0, 0, 0});
      g.output = out;

      g.usePalette = true;
      pal({
          {0.00f, 0.02f, 0.01f, 0.01f},
          {0.40f, 0.30f, 0.02f, 0.01f},
          {0.70f, 0.85f, 0.15f, 0.02f},
          {1.00f, 1.00f, 0.95f, 0.35f},
      });
    } break;

    case ProcGraphPreset::AlienCircuit: {
      // Voronoi distance with smoothstep to highlight "traces"
      int v = add(ProcNode{ProcNodeOp::Voronoi2D, -1, -1, -1, 18.0f, 1.00f, 0, 0, 0});
      int inv = add(ProcNode{ProcNodeOp::Invert, v, -1, -1, 0, 0, 0, 0, 0});
      int edge = add(ProcNode{ProcNodeOp::Smoothstep, inv, -1, -1, 0.55f, 0.85f, 0, 0, 0});
      int out = add(ProcNode{ProcNodeOp::Pow, edge, -1, -1, 1.6f, 0, 0, 0, 0});
      g.output = out;

      g.usePalette = true;
      pal({
          {0.00f, 0.01f, 0.03f, 0.02f},
          {0.55f, 0.02f, 0.18f, 0.08f},
          {0.80f, 0.15f, 0.90f, 0.55f},
          {1.00f, 0.92f, 1.00f, 0.95f},
      });
    } break;

    case ProcGraphPreset::Rocky: {
      // FBm base warped slightly by another noise (via Warp node).
      // Ridged Perlin FBm tends to make better "craggy" features.
      int fbm = add(ProcNode{ProcNodeOp::RidgedFbmPerlin2D, -1, -1, -1, 7.0f, 5.0f, 2.0f, 0.55f, 0});
      int warp = add(ProcNode{ProcNodeOp::Warp, fbm, -1, -1, 0.10f, 4.0f, 3.0f, 0.55f, 0});
      int pits = add(ProcNode{ProcNodeOp::Voronoi2D, -1, -1, -1, 9.0f, 1.00f, 0, 0, 0});
      int pitsInv = add(ProcNode{ProcNodeOp::Invert, pits, -1, -1, 0, 0, 0, 0, 0});
      int mix = add(ProcNode{ProcNodeOp::Add, warp, pitsInv, -1, 0, 0, 0, 0, 0});
      int out = add(ProcNode{ProcNodeOp::Clamp01, mix, -1, -1, 0, 0, 0, 0, 0});
      g.output = out;

      g.usePalette = true;
      pal({
          {0.00f, 0.06f, 0.05f, 0.04f},
          {0.55f, 0.25f, 0.20f, 0.15f},
          {0.85f, 0.55f, 0.52f, 0.45f},
          {1.00f, 0.80f, 0.78f, 0.70f},
      });
    } break;
  }

  // Ensure palette positions are monotonic.
  g.paletteCount = clampPalCount(g.paletteCount);
  for (int i = 0; i < g.paletteCount; ++i) {
    g.palette[i].pos = std::clamp(g.palette[i].pos, 0.0f, 1.0f);
  }
  std::sort(g.palette.begin(), g.palette.begin() + g.paletteCount,
            [](const ProcPaletteStop& a, const ProcPaletteStop& b) { return a.pos < b.pos; });

  return g;
}

ProceduralGraphBaker::~ProceduralGraphBaker() {
  if (vao_ != 0) {
    gl::DeleteVertexArrays(1, &vao_);
    vao_ = 0;
  }
}

core::u64 ProceduralGraphBaker::structureKey(const ProcGraph& g) {
  // Shader structure key: depends only on topology (op + connections + output index).
  // Parameters and seeds are uploaded via uniforms and do NOT require a recompile.
  core::u64 h = core::fnv1a64("ProcGraphShaderV2");
  const int nodeCount = clampNodeCount((int)g.nodes.size());
  h = core::hashCombine(h, (core::u64)nodeCount);
  for (int i = 0; i < nodeCount; ++i) {
    const auto& n = g.nodes[(std::size_t)i];
    const core::u64 packed =
        (core::u64)n.op |
        ((core::u64)(core::u32)(n.in0 + 2) << 8) |
        ((core::u64)(core::u32)(n.in1 + 2) << 20) |
        ((core::u64)(core::u32)(n.in2 + 2) << 32);
    h = core::hashCombine(h, packed);
  }
  h = core::hashCombine(h, (core::u64)(g.output + 1));
  return h;
}

std::string ProceduralGraphBaker::buildFragmentShader(const ProcGraph& g) {
  const int nodeCount = clampNodeCount((int)g.nodes.size());

  auto call = [&](int idx, const char* uvExpr = "uv") -> std::string {
    if (idx >= 0 && idx < nodeCount) {
      std::ostringstream oss;
      oss << "n" << idx << "(" << uvExpr << ")";
      return oss.str();
    }
    return "0.0";
  };

  std::ostringstream ss;
  ss << "#version 330 core\n";
  ss << "in vec2 vUv;\n";
  ss << "out vec4 oColor;\n\n";

  ss << "#define MAX_NODES " << kProcGraphMaxNodes << "\n";
  ss << "#define MAX_PAL " << kProcGraphMaxPaletteStops << "\n\n";

  ss << "uniform vec2 uResolution;\n";
  ss << "uniform float uTime;\n";
  ss << "uniform int uGraphSeed;\n";
  ss << "uniform float uDitherStrength;\n";
  ss << "uniform int uPackHeightInAlpha;\n";
  ss << "uniform vec4 uP[MAX_NODES];\n";
  ss << "uniform int  uS[MAX_NODES];\n";
  ss << "uniform int  uUsePalette;\n";
  ss << "uniform int  uPalCount;\n";
  ss << "uniform vec4 uPal[MAX_PAL]; // rgb + pos(a)\n\n";

  // Small helper library (value noise + FBm + Voronoi distance).
  ss << R"GLSL(
float hash21(vec2 p, int seed) {
  // Deterministic hash for grid points (no textures, no assets).
  // Note: sin-based hashing is OK here; this shader is used for baking.
  float h = dot(p, vec2(127.1, 311.7)) + float(seed) * 0.0131;
  return fract(sin(h) * 43758.5453123);
}

vec2 hash22(vec2 p, int seed) {
  float x = hash21(p + vec2(0.0, 0.0), seed);
  float y = hash21(p + vec2(19.19, 73.73), seed + 17);
  return vec2(x, y);
}

float fade(float t) {
  // Quintic fade
  return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

float valueNoise(vec2 x, int seed) {
  vec2 i = floor(x);
  vec2 f = fract(x);

  float a = hash21(i, seed);
  float b = hash21(i + vec2(1.0, 0.0), seed);
  float c = hash21(i + vec2(0.0, 1.0), seed);
  float d = hash21(i + vec2(1.0, 1.0), seed);

  vec2 u = vec2(fade(f.x), fade(f.y));
  return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

float fbm(vec2 x, int seed, int octaves, float lacunarity, float gain) {
  float sum = 0.0;
  float amp = 0.5;
  float freq = 1.0;

  for (int i = 0; i < 8; ++i) {
    if (i >= octaves) break;
    sum += amp * valueNoise(x * freq, seed + i * 1013);
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

vec2 grad2(vec2 p, int seed) {
  // Random unit vector per lattice point.
  float a = hash21(p, seed) * 6.28318530718;
  return vec2(cos(a), sin(a));
}

float perlinNoise(vec2 x, int seed) {
  vec2 i = floor(x);
  vec2 f = fract(x);

  vec2 g00 = grad2(i, seed);
  vec2 g10 = grad2(i + vec2(1.0, 0.0), seed);
  vec2 g01 = grad2(i + vec2(0.0, 1.0), seed);
  vec2 g11 = grad2(i + vec2(1.0, 1.0), seed);

  float d00 = dot(g00, f - vec2(0.0, 0.0));
  float d10 = dot(g10, f - vec2(1.0, 0.0));
  float d01 = dot(g01, f - vec2(0.0, 1.0));
  float d11 = dot(g11, f - vec2(1.0, 1.0));

  vec2 u = vec2(fade(f.x), fade(f.y));
  float a0 = mix(d00, d10, u.x);
  float b0 = mix(d01, d11, u.x);
  float n = mix(a0, b0, u.y);

  // Typical normalization for 2D Perlin; keeps output mostly in [-1,1].
  n *= 0.70710678;
  return clamp(0.5 + 0.5 * n, 0.0, 1.0);
}

float fbmPerlin(vec2 x, int seed, int octaves, float lacunarity, float gain) {
  float sum = 0.0;
  float amp = 0.5;
  float freq = 1.0;

  for (int i = 0; i < 8; ++i) {
    if (i >= octaves) break;
    sum += amp * perlinNoise(x * freq, seed + i * 1013);
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

float ridgedFbmPerlin(vec2 x, int seed, int octaves, float lacunarity, float gain) {
  float sum = 0.0;
  float amp = 0.5;
  float freq = 1.0;

  for (int i = 0; i < 8; ++i) {
    if (i >= octaves) break;
    float n = perlinNoise(x * freq, seed + i * 1013);
    float ridge = 1.0 - abs(2.0 * n - 1.0);
    ridge *= ridge;
    sum += amp * ridge;
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

float voronoiDist(vec2 x, int seed, float jitter) {
  vec2 i = floor(x);
  vec2 f = fract(x);

  float minD = 1e9;
  for (int y = -1; y <= 1; ++y) {
    for (int x0 = -1; x0 <= 1; ++x0) {
      vec2 g = vec2(float(x0), float(y));
      vec2 o = hash22(i + g, seed) - 0.5;
      o *= jitter;
      vec2 r = g + o - f;
      float d = dot(r, r);
      minD = min(minD, d);
    }
  }
  return sqrt(minD);
}

vec3 palette(float t) {
  t = clamp(t, 0.0, 1.0);
  int n = clamp(uPalCount, 2, MAX_PAL);

  // Ensure endpoints are usable even if user gives degenerate positions.
  vec4 s0 = uPal[0];
  vec4 s1 = uPal[n - 1];
  if (t <= s0.a) return s0.rgb;
  if (t >= s1.a) return s1.rgb;

  for (int i = 1; i < MAX_PAL; ++i) {
    if (i >= n) break;
    vec4 a = uPal[i - 1];
    vec4 b = uPal[i];
    if (t <= b.a) {
      float denom = max(1e-6, (b.a - a.a));
      float k = clamp((t - a.a) / denom, 0.0, 1.0);
      return mix(a.rgb, b.rgb, k);
    }
  }
  return uPal[n - 1].rgb;
}

// --- Dithering helpers ---
// The procedural baker writes to RGBA8. A tiny amount of deterministic noise helps reduce
// visible banding in smooth gradients.
uint hashU32(uint x) {
  x ^= x >> 16;
  x *= 2246822519u;
  x ^= x >> 13;
  x *= 3266489917u;
  x ^= x >> 16;
  return x;
}

float dither1(ivec2 px, int seed) {
  uint h = uint(px.x) * 1973u + uint(px.y) * 9277u + uint(seed) * 26699u;
  h = hashU32(h);
  // [0,1)
  return float(h & 0x00ffffffu) / float(0x01000000u);
}

vec3 ditherRgb(ivec2 px, int seed) {
  float r = dither1(px, seed);
  float g = dither1(px + ivec2(17, 11), seed + 19);
  float b = dither1(px + ivec2(37, 29), seed + 47);
  return vec3(r, g, b) - 0.5;
}
)GLSL";

  ss << "\n";

  // Node functions (scalar outputs). Structure is compiled; params/seeds are uniforms.
  for (int i = 0; i < nodeCount; ++i) {
    const auto& n = g.nodes[(std::size_t)i];
    ss << "float n" << i << "(vec2 uv) {\n";
    ss << "  vec4 p = uP[" << i << "];\n";
    ss << "  int seed = uS[" << i << "];\n";

    auto a = [&](const char* uvExpr = "uv") { return call(n.in0, uvExpr); };
    auto b = [&](const char* uvExpr = "uv") { return call(n.in1, uvExpr); };
    auto c = [&](const char* uvExpr = "uv") { return call(n.in2, uvExpr); };

    switch (n.op) {
      case ProcNodeOp::Constant:
        ss << "  return p.x;\n";
        break;

      case ProcNodeOp::UvX:
        ss << "  return uv.x;\n";
        break;

      case ProcNodeOp::UvY:
        ss << "  return uv.y;\n";
        break;

      case ProcNodeOp::Add:
        ss << "  return " << a() << " + " << b() << ";\n";
        break;

      case ProcNodeOp::Sub:
        ss << "  return " << a() << " - " << b() << ";\n";
        break;

      case ProcNodeOp::Mul:
        ss << "  return " << a() << " * " << b() << ";\n";
        break;

      case ProcNodeOp::Div:
        ss << "  return " << a() << " / max(1e-6, " << b() << ");\n";
        break;

      case ProcNodeOp::Min:
        ss << "  return min(" << a() << ", " << b() << ");\n";
        break;

      case ProcNodeOp::Max:
        ss << "  return max(" << a() << ", " << b() << ");\n";
        break;

      case ProcNodeOp::Abs:
        ss << "  return abs(" << a() << ");\n";
        break;

      case ProcNodeOp::Invert:
        ss << "  return 1.0 - (" << a() << ");\n";
        break;

      case ProcNodeOp::Fract:
        ss << "  return fract(" << a() << ");\n";
        break;

      case ProcNodeOp::Clamp01:
        ss << "  return clamp(" << a() << ", 0.0, 1.0);\n";
        break;

      case ProcNodeOp::Smoothstep:
        ss << "  return smoothstep(p.x, p.y, " << a() << ");\n";
        break;

      case ProcNodeOp::Pow:
        ss << "  return pow(max(" << a() << ", 0.0), max(p.x, 1e-6));\n";
        break;

      case ProcNodeOp::Sine:
        ss << "  return 0.5 + 0.5 * sin((" << a() << ") * p.x + p.y);\n";
        break;

      case ProcNodeOp::Noise2D: {
        ss << "  float freq = max(p.x, 1e-4);\n";
        ss << "  float speed = p.y;\n";
        ss << "  float amp = (p.z == 0.0) ? 1.0 : p.z;\n";
        ss << "  float bias = p.w;\n";
        ss << "  vec2 q = uv * freq + vec2(uTime * speed, uTime * speed * 0.37);\n";
        ss << "  return valueNoise(q, seed) * amp + bias;\n";
      } break;

      case ProcNodeOp::Fbm2D: {
        ss << "  float freq = max(p.x, 1e-4);\n";
        ss << "  int oct = clamp(int(p.y + 0.5), 1, 8);\n";
        ss << "  float lac = max(p.z, 1.0);\n";
        ss << "  float gain = clamp(p.w, 0.0, 1.0);\n";
        ss << "  return fbm(uv * freq, seed, oct, lac, gain);\n";
      } break;

      case ProcNodeOp::Perlin2D: {
        ss << "  float freq = max(p.x, 1e-4);\n";
        ss << "  float speed = p.y;\n";
        ss << "  float amp = (p.z == 0.0) ? 1.0 : p.z;\n";
        ss << "  float bias = p.w;\n";
        ss << "  vec2 q = uv * freq + vec2(uTime * speed, uTime * speed * 0.37);\n";
        ss << "  return perlinNoise(q, seed) * amp + bias;\n";
      } break;

      case ProcNodeOp::FbmPerlin2D: {
        ss << "  float freq = max(p.x, 1e-4);\n";
        ss << "  int oct = clamp(int(p.y + 0.5), 1, 8);\n";
        ss << "  float lac = max(p.z, 1.0);\n";
        ss << "  float gain = clamp(p.w, 0.0, 1.0);\n";
        ss << "  return fbmPerlin(uv * freq, seed, oct, lac, gain);\n";
      } break;

      case ProcNodeOp::RidgedFbmPerlin2D: {
        ss << "  float freq = max(p.x, 1e-4);\n";
        ss << "  int oct = clamp(int(p.y + 0.5), 1, 8);\n";
        ss << "  float lac = max(p.z, 1.0);\n";
        ss << "  float gain = clamp(p.w, 0.0, 1.0);\n";
        ss << "  return ridgedFbmPerlin(uv * freq, seed, oct, lac, gain);\n";
      } break;

      case ProcNodeOp::Voronoi2D: {
        ss << "  float freq = max(p.x, 1e-4);\n";
        ss << "  float jitter = clamp(p.y, 0.0, 1.0);\n";
        ss << "  float d = voronoiDist(uv * freq, seed, jitter);\n";
        ss << "  // Normalize a bit so distance fields are easier to palette-map.\n";
        ss << "  return clamp(d * 1.41421356, 0.0, 1.0);\n";
      } break;

      case ProcNodeOp::Warp: {
        // p0 strength, p1 freq, p2 octaves, p3 gain
        ss << "  float strength = p.x;\n";
        ss << "  float freq = max(p.y, 1e-4);\n";
        ss << "  int oct = clamp(int(p.z + 0.5), 1, 8);\n";
        ss << "  float gain = clamp(p.w, 0.0, 1.0);\n";
        ss << "  vec2 w;\n";
        ss << "  w.x = fbm(uv * freq + vec2(12.3, 4.56), seed, oct, 2.0, gain);\n";
        ss << "  w.y = fbm(uv * freq + vec2(78.9, 1.23), seed + 17, oct, 2.0, gain);\n";
        ss << "  w = (w - 0.5) * strength;\n";
        ss << "  return " << a("uv + w") << ";\n";
      } break;

      case ProcNodeOp::Pan: {
        // p0 vx, p1 vy
        ss << "  vec2 v = vec2(p.x, p.y);\n";
        ss << "  return " << a("uv + v * uTime") << ";\n";
      } break;

      default:
        ss << "  return 0.0;\n";
        break;
    }

    ss << "}\n\n";
  }

  // Main.
  ss << "void main() {\n";
  ss << "  vec2 uv = vUv;\n";

  if (nodeCount <= 0) {
    ss << "  float t = 0.0;\n";
  } else {
    int outIdx = g.output;
    if (outIdx < 0 || outIdx >= nodeCount) outIdx = nodeCount - 1;
    ss << "  float t = n" << outIdx << "(uv);\n";
  }

  ss << "  t = clamp(t, 0.0, 1.0);\n";
  ss << "  vec3 col = (uUsePalette != 0) ? palette(t) : vec3(t);\n";
  ss << "  if (uDitherStrength > 0.0) {\n";
  ss << "    vec3 dn = ditherRgb(ivec2(gl_FragCoord.xy), uGraphSeed) * (uDitherStrength / 255.0);\n";
  ss << "    col = clamp(col + dn, 0.0, 1.0);\n";
  ss << "  }\n";
  ss << "  float a = (uPackHeightInAlpha != 0) ? t : 1.0;\n";
  ss << "  oColor = vec4(col, a);\n";
  ss << "}\n";

  return ss.str();
}

bool ProceduralGraphBaker::ensureShader(const ProcGraph& g, std::string* outError) {
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

bool ProceduralGraphBaker::bake(const ProcGraph& g, int width, int height, float timeSec, std::string* outError) {
  width = std::max(1, width);
  height = std::max(1, height);

  // Ensure FBO / shader are available.
  // Note: the bake is a fullscreen pass; we don't need depth.
  RenderTarget2DDesc rtDesc{};
  rtDesc.hasDepth = false;
  rtDesc.generateMips = generateMips_;

  std::string rtErr;
  if (!target_.init(width, height, rtDesc, &rtErr)) {
    if (outError) *outError = std::string("ProceduralGraphBaker: ") + rtErr;
    return false;
  }

  if (!ensureShader(g, outError)) return false;

  const int nodeCount = clampNodeCount((int)g.nodes.size());

  // Pack parameters + per-node seeds.
  packedParams_.fill(0.0f);
  packedSeeds_.fill(0);

  for (int i = 0; i < nodeCount; ++i) {
    const auto& n = g.nodes[(std::size_t)i];
    packedParams_[i * 4 + 0] = n.p0;
    packedParams_[i * 4 + 1] = n.p1;
    packedParams_[i * 4 + 2] = n.p2;
    packedParams_[i * 4 + 3] = n.p3;

    // Deterministic per-node seed: mix graph seed + node seed + index.
    const core::u64 s = core::hashCombine(g.seed, core::hashCombine((core::u64)i, n.seed));
    packedSeeds_[i] = (int)(s & 0x7fffffff);
  }

  // Pack palette.
  packedPalette_.fill(0.0f);
  int palCount = clampPalCount(g.paletteCount);

  std::array<ProcPaletteStop, kProcGraphMaxPaletteStops> pal = g.palette;
  for (int i = 0; i < palCount; ++i) pal[i].pos = std::clamp(pal[i].pos, 0.0f, 1.0f);
  std::sort(pal.begin(), pal.begin() + palCount,
            [](const ProcPaletteStop& a, const ProcPaletteStop& b) { return a.pos < b.pos; });

  for (int i = 0; i < kProcGraphMaxPaletteStops; ++i) {
    const ProcPaletteStop& s = pal[std::min(i, palCount - 1)];
    packedPalette_[i * 4 + 0] = std::clamp(s.r, 0.0f, 1.0f);
    packedPalette_[i * 4 + 1] = std::clamp(s.g, 0.0f, 1.0f);
    packedPalette_[i * 4 + 2] = std::clamp(s.b, 0.0f, 1.0f);
    packedPalette_[i * 4 + 3] = s.pos;
  }

  // --- Save GL state we touch (important: we bake during ImGui UI rendering). ---
  GLint prevFbo = 0;
  GLint prevProg = 0;
  GLint prevVao = 0;
  GLint prevViewport[4] = {0, 0, 0, 0};

  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
  glGetIntegerv(GL_CURRENT_PROGRAM, &prevProg);
  glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &prevVao);
  glGetIntegerv(GL_VIEWPORT, prevViewport);

  const GLboolean depthEnabled = glIsEnabled(GL_DEPTH_TEST);
  const GLboolean blendEnabled = glIsEnabled(GL_BLEND);
  const GLboolean cullEnabled = glIsEnabled(GL_CULL_FACE);

  // --- Render ---
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
  shader_.setUniform1i("uUsePalette", g.usePalette ? 1 : 0);
  shader_.setUniform1i("uPalCount", palCount);
  shader_.setUniform1i("uGraphSeed", (int)(g.seed & 0x7fffffffULL));
  shader_.setUniform1f("uDitherStrength", ditherStrength_);
  shader_.setUniform1i("uPackHeightInAlpha", packHeightInAlpha_ ? 1 : 0);

  // Array uniforms: upload in one call each.
  shader_.setUniform4fv("uP[0]", kProcGraphMaxNodes, packedParams_.data());
  shader_.setUniform1iv("uS[0]", kProcGraphMaxNodes, packedSeeds_.data());
  shader_.setUniform4fv("uPal[0]", kProcGraphMaxPaletteStops, packedPalette_.data());

  gl::BindVertexArray(vao_);
  glDrawArrays(GL_TRIANGLES, 0, 3);

  target_.end();

  auto t1 = std::chrono::high_resolution_clock::now();
  stats_.drawMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

  stats_.mipsGenerated = false;
  stats_.mipGenMs = 0.0;

  if (generateMips_) {
    auto tm0 = std::chrono::high_resolution_clock::now();
    target_.generateMips();
    auto tm1 = std::chrono::high_resolution_clock::now();
    stats_.mipsGenerated = true;
    stats_.mipGenMs = std::chrono::duration<double, std::milli>(tm1 - tm0).count();
  }

  // --- Restore GL state ---
  gl::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
  glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);

  gl::UseProgram((GLuint)prevProg);
  gl::BindVertexArray((GLuint)prevVao);

  if (depthEnabled) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
  if (blendEnabled) glEnable(GL_BLEND); else glDisable(GL_BLEND);
  if (cullEnabled) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);

  return true;
}

bool saveProcGraphToFile(const ProcGraph& g, const std::string& path, std::string* outError) {
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

  const int nodeCount = clampNodeCount((int)g.nodes.size());
  const int palCount = clampPalCount(g.paletteCount);

  f.setf(std::ios::fixed);
  f.precision(6);

  f << "StellarForgeProcGraph 1\n";
  f << "seed " << (unsigned long long)g.seed << "\n";
  f << "output " << g.output << "\n";
  f << "usePalette " << (g.usePalette ? 1 : 0) << "\n";
  f << "paletteCount " << g.paletteCount << "\n";
  f << "paletteStops " << palCount << "\n";
  for (int i = 0; i < palCount; ++i) {
    const ProcPaletteStop& s = g.palette[(std::size_t)i];
    f << "stop " << s.pos << " " << s.r << " " << s.g << " " << s.b << "\n";
  }

  f << "nodes " << nodeCount << "\n";
  for (int i = 0; i < nodeCount; ++i) {
    const ProcNode& n = g.nodes[(std::size_t)i];
    f << "node " << procNodeOpName(n.op)
      << " " << n.in0 << " " << n.in1 << " " << n.in2
      << " " << n.p0 << " " << n.p1 << " " << n.p2 << " " << n.p3
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

bool loadProcGraphFromFile(const std::string& path, ProcGraph& out, std::string* outError) {
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
  if (header != "StellarForgeProcGraph") {
    if (outError) *outError = "Not a StellarForgeProcGraph file.";
    return false;
  }

  ProcGraph g{};
  g.seed = 0xC0FFEEULL;
  g.output = -1;
  g.usePalette = true;
  g.paletteCount = 4;
  // Default palette: grayscale ramp.
  for (int i = 0; i < kProcGraphMaxPaletteStops; ++i) {
    const float t = (float)i / (float)std::max(1, kProcGraphMaxPaletteStops - 1);
    g.palette[i] = ProcPaletteStop{t, t, t, t};
  }
  g.nodes.clear();

  (void)version; // Currently only v1.

  std::vector<ProcPaletteStop> stops;

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
    } else if (key == "usepalette") {
      int v = 1;
      ss >> v;
      g.usePalette = (v != 0);
    } else if (key == "palettecount") {
      ss >> g.paletteCount;
    } else if (key == "stop" || key == "palstop" || key == "pal") {
      ProcPaletteStop s{};
      if (ss >> s.pos >> s.r >> s.g >> s.b) {
        stops.push_back(s);
      }
    } else if (key == "node") {
      ProcNode n{};
      std::string opTok;
      if (!(ss >> opTok)) continue;

      if (auto op = procNodeOpFromString(opTok)) {
        n.op = *op;
      } else {
        int opInt = 0;
        if (parseI32(opTok, opInt)) {
          n.op = (ProcNodeOp)std::clamp(opInt, 0, 255);
        } else {
          // Unknown op; skip the node.
          continue;
        }
      }

      ss >> n.in0 >> n.in1 >> n.in2;
      ss >> n.p0 >> n.p1 >> n.p2 >> n.p3;
      std::string seedTok;
      if (ss >> seedTok) {
        core::u64 tmp = 0;
        if (parseU64(seedTok, tmp)) n.seed = tmp;
      }

      if ((int)g.nodes.size() < kProcGraphMaxNodes) {
        g.nodes.push_back(n);
      }
    }
  }

  // Apply palette from stops (sorted by pos).
  g.paletteCount = clampPalCount(g.paletteCount);
  if (!stops.empty()) {
    std::sort(stops.begin(), stops.end(), [](const ProcPaletteStop& a, const ProcPaletteStop& b) {
      return a.pos < b.pos;
    });

    const int n = std::min((int)stops.size(), kProcGraphMaxPaletteStops);
    for (int i = 0; i < n; ++i) {
      g.palette[i] = stops[(std::size_t)i];
    }
    // Fill remainder for safety.
    for (int i = n; i < kProcGraphMaxPaletteStops; ++i) {
      g.palette[i] = g.palette[std::max(0, n - 1)];
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

} // namespace stellar::render
