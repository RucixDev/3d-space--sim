#include "ProceduralShaderLabWindow.h"

#include "imgui.h"

#include "stellar/core/Hash.h"
#include "stellar/math/Math.h"
#include "stellar/math/Vec3.h"

#include <SDL_opengl.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace stellar::game {

namespace fs = std::filesystem;

namespace {

// -----------------------------
// Small utilities
// -----------------------------

template <size_t N>
static std::string_view cbufView(const std::array<char, N>& buf) {
  size_t len = 0;
  while (len < N && buf[len] != '\0') ++len;
  return std::string_view(buf.data(), len);
}

template <size_t N>
static void copyToBuf(std::array<char, N>& dst, std::string_view src) {
  const size_t n = std::min(N - 1, src.size());
  std::memcpy(dst.data(), src.data(), n);
  dst[n] = '\0';
  if (n + 1 < N) {
    std::memset(dst.data() + n + 1, 0, N - (n + 1));
  }
}

static int seed32(core::u64 s) {
  // Stable-ish mapping to 32-bit for GLSL.
  const core::u64 h = core::hashCombine(s, core::fnv1a64("shaderLab"));
  return (int)(h & 0x7FFFFFFFu);
}

static std::uint64_t fileStamp(const fs::path& p) {
  std::error_code ec;
  if (!fs::exists(p, ec)) return 0;
  const auto t = fs::last_write_time(p, ec);
  if (ec) return 0;
  // Coarse stamp is good enough.
  return (std::uint64_t)t.time_since_epoch().count();
}

static bool readTextFile(const fs::path& p, std::string& outText, std::string& outErr) {
  std::ifstream f(p, std::ios::binary);
  if (!f) {
    outErr = "Failed to open file.";
    return false;
  }

  std::ostringstream ss;
  ss << f.rdbuf();
  outText = ss.str();
  return true;
}

static bool writeTextFile(const fs::path& p, std::string_view text, std::string& outErr) {
  std::ofstream f(p, std::ios::binary);
  if (!f) {
    outErr = "Failed to open file for writing.";
    return false;
  }
  f.write(text.data(), (std::streamsize)text.size());
  if (!f) {
    outErr = "Write failed.";
    return false;
  }
  return true;
}

static bool icontains(std::string_view hay, std::string_view needle) {
  if (needle.empty()) return true;
  if (hay.empty()) return false;
  if (needle.size() > hay.size()) return false;

  auto lo = [](unsigned char c) { return (char)std::tolower(c); };

  for (size_t i = 0; i + needle.size() <= hay.size(); ++i) {
    bool ok = true;
    for (size_t j = 0; j < needle.size(); ++j) {
      if (lo((unsigned char)hay[i + j]) != lo((unsigned char)needle[j])) {
        ok = false;
        break;
      }
    }
    if (ok) return true;
  }
  return false;
}

static core::u64 splitmix64(core::u64& state) {
  // Simple deterministic RNG suitable for procedural textures.
  // https://prng.di.unimi.it/splitmix64.c (public domain reference implementation)
  core::u64 z = (state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

static float rand01(core::u64& state) {
  // 24-bit mantissa in [0,1).
  const core::u64 r = splitmix64(state);
  const unsigned int x = (unsigned int)(r & 0xFFFFFFu);
  return (float)x / 16777216.0f;
}

// -----------------------------
// Pass / channel mapping helpers
// -----------------------------

static const char* passLabel(int passIdx) {
  switch (passIdx) {
    case 0: return "Buffer A";
    case 1: return "Buffer B";
    case 2: return "Buffer C";
    case 3: return "Buffer D";
    default: return "Image";
  }
}

static stellar::render::ShaderToyPass passEnum(int passIdx) {
  passIdx = std::clamp(passIdx, 0, ProceduralShaderLabWindowState::kPassCount - 1);
  return static_cast<stellar::render::ShaderToyPass>(passIdx);
}

static const char* channelLabel(stellar::render::ShaderToyChannelSource s) {
  using S = stellar::render::ShaderToyChannelSource;
  switch (s) {
    case S::None: return "None";
    case S::BufferA: return "Buffer A";
    case S::BufferB: return "Buffer B";
    case S::BufferC: return "Buffer C";
    case S::BufferD: return "Buffer D";
    case S::External0: return "Checker";
    case S::External1: return "Noise";
    case S::External2: return "Fluid";
    case S::External3: return "Flow";
  }
  return "None";
}

static stellar::render::ShaderToyChannelSource channelFromToken(std::string_view t) {
  // Accept a few spellings.
  auto up = [](std::string_view s) {
    std::string r;
    r.reserve(s.size());
    for (char c : s) r.push_back((char)std::toupper((unsigned char)c));
    return r;
  };

  const std::string u = up(t);
  using S = stellar::render::ShaderToyChannelSource;
  if (u == "NONE" || u == "-" || u == "0") return S::None;
  if (u == "A" || u == "BUFFERA" || u == "BUFFER_A") return S::BufferA;
  if (u == "B" || u == "BUFFERB" || u == "BUFFER_B") return S::BufferB;
  if (u == "C" || u == "BUFFERC" || u == "BUFFER_C") return S::BufferC;
  if (u == "D" || u == "BUFFERD" || u == "BUFFER_D") return S::BufferD;

  // External textures (Shader Lab binds these to useful built-ins).
  if (u == "EXT0" || u == "EXTERNAL0" || u == "CHECKER" || u == "CHK") return S::External0;
  if (u == "EXT1" || u == "EXTERNAL1" || u == "NOISE" || u == "RNG") return S::External1;
  if (u == "EXT2" || u == "EXTERNAL2" || u == "FLUID") return S::External2;
  if (u == "EXT3" || u == "EXTERNAL3" || u == "FLOW" || u == "VEL" || u == "VELOCITY") return S::External3;
  return S::None;
}

static const char* channelToken(stellar::render::ShaderToyChannelSource s) {
  using S = stellar::render::ShaderToyChannelSource;
  switch (s) {
    case S::None: return "None";
    case S::BufferA: return "A";
    case S::BufferB: return "B";
    case S::BufferC: return "C";
    case S::BufferD: return "D";
    case S::External0: return "Checker";
    case S::External1: return "Noise";
    case S::External2: return "Fluid";
    case S::External3: return "Flow";
  }
  return "None";
}

// -----------------------------
// ShaderToyGraph file format (.stoy)
// -----------------------------

static bool looksLikeGraphFile(const std::string& text) {
  return text.find("PASS ") != std::string::npos && text.find("CODE_BEGIN") != std::string::npos;
}

static int passIndexFromName(std::string_view name) {
  auto up = [](std::string_view s) {
    std::string r;
    r.reserve(s.size());
    for (char c : s) r.push_back((char)std::toupper((unsigned char)c));
    return r;
  };
  const std::string u = up(name);
  if (u == "A" || u == "BUFFERA" || u == "BUFFER A") return 0;
  if (u == "B" || u == "BUFFERB" || u == "BUFFER B") return 1;
  if (u == "C" || u == "BUFFERC" || u == "BUFFER C") return 2;
  if (u == "D" || u == "BUFFERD" || u == "BUFFER D") return 3;
  if (u == "IMAGE" || u == "IMG") return 4;
  return -1;
}

static bool loadGraphText(ProceduralShaderLabWindowState& st, const std::string& text, std::string& outErr) {
  if (!looksLikeGraphFile(text)) {
    // Fallback: treat as a single Image snippet.
    copyToBuf(st.code[4], text);
    st.passes[4].dirty = true;
    st.presetIndex = 0;
    st.requestCompileAll = true;
    st.requestResetBuffers = true;
    return true;
  }

  int curPass = -1;
  bool inCode = false;
  std::vector<std::string> codeLines;

  bool inParams = false;
  std::unordered_map<std::string, std::array<float, 4>> loadedParams;

  auto flushCode = [&]() {
    if (curPass < 0 || curPass >= ProceduralShaderLabWindowState::kPassCount) return;
    std::string joined;
    for (size_t i = 0; i < codeLines.size(); ++i) {
      joined += codeLines[i];
      if (i + 1 < codeLines.size()) joined += "\n";
    }
    copyToBuf(st.code[curPass], joined);
    st.passes[curPass].dirty = true;
    codeLines.clear();
  };

  std::istringstream iss(text);
  std::string line;
  while (std::getline(iss, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();

    if (inCode) {
      if (line == "CODE_END") {
        inCode = false;
        flushCode();
        continue;
      }
      codeLines.push_back(line);
      continue;
    }

    // Params section
    if (inParams) {
      if (line == "PARAMS_END") {
        inParams = false;
        continue;
      }
      if (line.rfind("PARAM ", 0) == 0) {
        // PARAM <name> <x> <y> <z> <w>
        std::istringstream ls(line);
        std::string cmd;
        std::string name;
        float x = 0, y = 0, z = 0, w = 0;
        ls >> cmd >> name >> x >> y >> z >> w;
        if (!name.empty()) {
          loadedParams[name] = {x, y, z, w};
        }
      }
      continue;
    }

    // Skip comments / empty lines.
    const char* c = line.c_str();
    while (*c && std::isspace((unsigned char)*c)) ++c;
    if (*c == '\0' || *c == '#') continue;

    if (line == "PARAMS_BEGIN") {
      inParams = true;
      continue;
    }

    if (line.rfind("PASS ", 0) == 0) {
      flushCode();
      const std::string name = line.substr(5);
      curPass = passIndexFromName(name);
      if (curPass < 0) {
        outErr = "Unknown pass name: " + name;
        return false;
      }
      continue;
    }

    if (line.rfind("ENABLED ", 0) == 0) {
      if (curPass < 0) continue;
      const int v = std::atoi(line.c_str() + 8);
      st.passes[curPass].enabled = (v != 0);
      continue;
    }

    if (line.rfind("SCALE ", 0) == 0) {
      if (curPass < 0) continue;
      const int s = std::atoi(line.c_str() + 6);
      st.passes[curPass].resolutionScale = (s <= 1) ? 1 : ((s <= 2) ? 2 : 4);
      continue;
    }

    if (line.rfind("CHANNELS ", 0) == 0) {
      if (curPass < 0) continue;
      std::istringstream ls(line.substr(9));
      for (int i = 0; i < 4; ++i) {
        std::string tok;
        if (!(ls >> tok)) tok = "None";
        st.passes[curPass].channels[i] = channelFromToken(tok);
      }
      continue;
    }

    if (line == "CODE_BEGIN") {
      inCode = true;
      codeLines.clear();
      continue;
    }
  }

  // If file ended inside a code block.
  if (inCode) {
    outErr = "Unexpected EOF inside CODE_BEGIN block.";
    return false;
  }

  if (inParams) {
    outErr = "Unexpected EOF inside PARAMS_BEGIN block.";
    return false;
  }

  flushCode();

  // Parse params from the loaded code and apply any serialized overrides.
  {
    std::array<std::string_view, ProceduralShaderLabWindowState::kPassCount> src{};
    for (int p = 0; p < ProceduralShaderLabWindowState::kPassCount; ++p) {
      src[p] = cbufView(st.code[p]);
    }
    st.params = stellar::render::parseShaderToyParamsFromSources(src, /*preserveValuesFrom=*/nullptr);
    for (const auto& [name, v] : loadedParams) {
      (void)st.params.setValue(name, v);
    }
  }

  st.requestCompileAll = true;
  st.requestResetBuffers = true;
  return true;
}

static std::string saveGraphText(const ProceduralShaderLabWindowState& st) {
  std::ostringstream ss;
  ss << "# stellar ShaderToyGraph v2\n";
  ss << "# PASS <Buffer A|B|C|D|Image>\n";
  ss << "# ENABLED <0|1> (buffers only; image always enabled)\n";
  ss << "# SCALE <1|2|4>\n";
  ss << "# CHANNELS <None|A|B|C|D|Checker|Noise|Fluid|Flow> x4\n";
  ss << "# PARAMS_BEGIN/END are optional overrides for `// @param ...` values\n";
  ss << "# PARAM <name> <x> <y> <z> <w>\n";
  ss << "\n";

  for (int p = 0; p < ProceduralShaderLabWindowState::kPassCount; ++p) {
    ss << "PASS " << passLabel(p) << "\n";
    const bool enabled = (p == 4) ? true : st.passes[p].enabled;
    ss << "ENABLED " << (enabled ? 1 : 0) << "\n";
    ss << "SCALE " << st.passes[p].resolutionScale << "\n";
    ss << "CHANNELS";
    for (int c = 0; c < 4; ++c) {
      ss << " " << channelToken(st.passes[p].channels[c]);
    }
    ss << "\n";
    ss << "CODE_BEGIN\n";
    ss << cbufView(st.code[p]) << "\n";
    ss << "CODE_END\n\n";
  }

  if (!st.params.empty()) {
    ss << "PARAMS_BEGIN\n";
    const int n = (int)st.params.defs.size();
    for (int i = 0; i < n && i < (int)st.params.values.size(); ++i) {
      const auto& d = st.params.defs[i];
      const auto& v = st.params.values[i];
      ss << "PARAM " << d.name << " " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << "\n";
    }
    ss << "PARAMS_END\n";
  }

  return ss.str();
}

// -----------------------------
// Presets
// -----------------------------

struct PresetPass {
  const char* code;
  bool enabled;
  int scale;
  std::array<stellar::render::ShaderToyChannelSource, 4> channels;
};

struct GraphPreset {
  const char* name;
  PresetPass pass[ProceduralShaderLabWindowState::kPassCount];
};

static const char* kBlackSnippet = R"GLSL(
vec4 shaderMain(vec2 uv) {
  return vec4(0.0);
}
)GLSL";

static const char* kInkflowBufferA = R"GLSL(
// Inkflow (Buffer A)
//
// Demonstrates external iChannel textures:
//   iChannel1: Flow (External3)  -- velocity encoded in RG (0.5 = 0)
//   iChannel2: Fluid dye (External2)
//   iChannel3: Noise (External1)
//
// @group Inkflow
// @param float ink_advect 1.20 0.0 5.0 0.01
// @param float ink_flowAmp 1.00 0.0 3.0 0.01
// @param float ink_decay 0.985 0.0 1.0 0.0005
// @param float ink_inject 0.22 0.0 1.0 0.01
// @param float ink_noise 0.08 0.0 0.5 0.001
// @endgroup

vec2 flowAt(vec2 uv) {
  // External3 encodes velocity in RG, 0.5=center.
  vec2 f = texture(iChannel1, uv).rg;
  return (f * 2.0 - 1.0);
}

vec4 shaderMain(vec2 uv) {
  // Bootstrap with the dye field.
  if (iFrame < 1) {
    return vec4(texture(iChannel2, uv).rgb, 1.0);
  }

  float dt = clamp(iTimeDelta, 0.0, 0.1);

  vec2 flow = ink_flowAmp * flowAt(uv);
  // This constant is tuned against the CPU fluid's default velocity encoding.
  vec2 uvPrev = uv - dt * ink_advect * flow * 0.25;

  vec3 prev = texture(iChannel0, uvPrev).rgb;
  vec3 dye  = texture(iChannel2, uv).rgb;

  float n = texture(iChannel3, uv * 2.0 + 0.05 * iTime).r - 0.5;

  vec3 col = prev * ink_decay;
  col = mix(col, dye, ink_inject);
  col += ink_noise * n;
  col = max(col, vec3(0.0));
  return vec4(col, 1.0);
}
)GLSL";

static const char* kInkflowImage = R"GLSL(
// Inkflow (Image)
// Visualizes the feedback buffer with a simple palette + vignette.

// @group Look
// @param float look_contrast 1.25 0.2 2.5 0.01
// @param float look_vignette 0.35 0.0 1.2 0.01
// @endgroup

vec4 shaderMain(vec2 uv) {
  vec3 ink = texture(iChannel0, uv).rgb;
  float l = clamp(dot(ink, vec3(0.3333)), 0.0, 1.0);

  // A palette that reads well for both faint and saturated dye.
  vec3 col = palette(l,
                     vec3(0.20, 0.22, 0.25),
                     vec3(0.55, 0.55, 0.55),
                     vec3(1.00, 1.00, 1.00),
                     vec3(0.00, 0.10, 0.55));

  // Bring back some of the original color.
  col = mix(col, ink, 0.25);
  col *= look_contrast;

  // Soft vignette.
  vec2 p = uv * 2.0 - 1.0;
  p.x *= iAspect;
  float r2 = dot(p, p);
  float vig = smoothstep(1.2, look_vignette, r2);
  col *= vig;

  return vec4(tonemapSimple(col), 1.0);
}
)GLSL";

static const GraphPreset kPresets[] = {
  {
    "Template (single-pass)",
    {
      // Buffer A
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Buffer B
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Buffer C
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Buffer D
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Image
      {
        R"GLSL(
// Wrapper-injected uniforms (subset of ShaderToy):
//   uniform vec2  iResolution;   // (w,h) pixels
//   uniform float iTime;         // seconds
//   uniform float iTimeDelta;    // seconds
//   uniform int   iFrame;        // frame counter
//   uniform vec4  iMouse;        // x,y current, z,w click origin (pixels)
//   uniform int   iSeed;
//   uniform int   iPass;         // 0..3 buffers, 4 image
//   uniform sampler2D iChannel0..3;
//   uniform vec3  iChannelResolution[4];
//
// And helper functions:
//   noise2(), fbm2(), worley2(), warp2(), palette(), sdSphere(), sdBox(), rayDirFromUv(), tonemapSimple()

vec4 shaderMain(vec2 uv) {
  // uv in [0,1]
  vec2 p = uv * 2.0 - 1.0;
  p.x *= iAspect;

  float v = 0.5 + 0.5 * sin(iTime + 5.0 * length(p));
  vec3 col = palette(v,
                     vec3(0.35, 0.35, 0.35),
                     vec3(0.55, 0.55, 0.55),
                     vec3(1.00, 1.00, 1.00),
                     vec3(0.00, 0.33, 0.67));
  return vec4(col, 1.0);
}
)GLSL",
        true,
        1,
        {stellar::render::ShaderToyChannelSource::None,
         stellar::render::ShaderToyChannelSource::None,
         stellar::render::ShaderToyChannelSource::None,
         stellar::render::ShaderToyChannelSource::None},
      },
    },
  },

  {
    "Reaction Diffusion (Buffer A feedback)",
    {
      // Buffer A (simulation)
      {
        R"GLSL(
// Gray-Scott reaction diffusion.
// State packed as: (u, v) in RG.

// @group Simulation
// @param float rd_F 0.0367 0.0 0.08 0.0001
// @param float rd_K 0.0649 0.0 0.08 0.0001
// @param float rd_du 1.0 0.0 2.0 0.01
// @param float rd_dv 0.5 0.0 2.0 0.01
// @param float rd_paintRadius 0.08 0.0 0.25 0.001
// @param float rd_paintStrength 1.0 0.0 1.0 0.01
// @param float rd_modAmp 0.015 0.0 0.05 0.0005
// @endgroup

vec2 laplace(vec2 uv) {
  vec2 px = 1.0 / iResolution;

  vec2 c  = texture(iChannel0, uv).xy;
  vec2 n  = texture(iChannel0, uv + vec2(0.0,  px.y)).xy;
  vec2 s  = texture(iChannel0, uv + vec2(0.0, -px.y)).xy;
  vec2 e  = texture(iChannel0, uv + vec2( px.x, 0.0)).xy;
  vec2 w  = texture(iChannel0, uv + vec2(-px.x, 0.0)).xy;
  vec2 ne = texture(iChannel0, uv + vec2( px.x,  px.y)).xy;
  vec2 nw = texture(iChannel0, uv + vec2(-px.x,  px.y)).xy;
  vec2 se = texture(iChannel0, uv + vec2( px.x, -px.y)).xy;
  vec2 sw = texture(iChannel0, uv + vec2(-px.x, -px.y)).xy;

  // 9-tap kernel
  vec2 lap = (n + s + e + w) * 0.20 + (ne + nw + se + sw) * 0.05 - c;
  return lap;
}

vec4 shaderMain(vec2 uv) {
  // Init on frame 0.
  if (iFrame < 1) {
    vec2 st = vec2(1.0, 0.0); // u=1, v=0

    // Seed with sparse spots.
    float r = fract(sin(dot((uv * iResolution) + float(iSeed), vec2(12.9898, 78.233))) * 43758.5453);
    float s = step(0.9965, r);
    st = mix(st, vec2(0.0, 1.0), s);

    // A central blob.
    float d = length(uv - vec2(0.5));
    st = mix(st, vec2(0.0, 1.0), smoothstep(0.06, 0.03, d));

    return vec4(st, 0.0, 1.0);
  }

  vec2 st = texture(iChannel0, uv).xy;
  float u = st.x;
  float v = st.y;

  // Mouse paints more V.
  if (iMouse.z > 0.0) {
    vec2 m = iMouse.xy / iResolution;
    float d = length(uv - m);
    float k = smoothstep(rd_paintRadius, 0.0, d) * rd_paintStrength;
    v = mix(v, 1.0, clamp(k, 0.0, 1.0));
    u = mix(u, 0.0, clamp(k, 0.0, 1.0));
  }

  vec2 lap = laplace(uv);

  // Parameters.
  float F = rd_F;
  float K = rd_K;

  // A tiny spatial modulation so it doesn't settle too quickly.
  float modF = rd_modAmp * fbm2(uv * 3.0 + 0.02 * iTime);
  F += modF;

  float du = rd_du;
  float dv = rd_dv;

  float uvv = u * v * v;

  float dt = clamp(iTimeDelta, 0.0, 0.1);

  u += (du * lap.x - uvv + F * (1.0 - u)) * dt;
  v += (dv * lap.y + uvv - (F + K) * v) * dt;

  u = clamp(u, 0.0, 1.0);
  v = clamp(v, 0.0, 1.0);

  return vec4(u, v, 0.0, 1.0);
}
)GLSL",
        true,
        1,
        {stellar::render::ShaderToyChannelSource::BufferA,
         stellar::render::ShaderToyChannelSource::None,
         stellar::render::ShaderToyChannelSource::None,
         stellar::render::ShaderToyChannelSource::None},
      },
      // Buffer B
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Buffer C
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Buffer D
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Image (visualize)
      {
        R"GLSL(
// @group Visualization
// @param float rd_edgeGain 3.0 0.0 10.0 0.1
// @param color rd_edgeColor 1.0 1.0 1.0
// @endgroup

vec4 shaderMain(vec2 uv) {
  vec2 st = texture(iChannel0, uv).xy;
  float u = st.x;
  float v = st.y;

  float edge = clamp((v - u) * rd_edgeGain, 0.0, 1.0);
  vec3 base = palette(v,
                      vec3(0.2, 0.2, 0.25),
                      vec3(0.6, 0.6, 0.65),
                      vec3(1.0, 1.0, 1.0),
                      vec3(0.0, 0.2, 0.6));

  vec3 col = mix(base, rd_edgeColor, edge);
  col = tonemapSimple(col);
  return vec4(col, 1.0);
}
)GLSL",
        true,
        1,
        {stellar::render::ShaderToyChannelSource::BufferA,
         stellar::render::ShaderToyChannelSource::None,
         stellar::render::ShaderToyChannelSource::None,
         stellar::render::ShaderToyChannelSource::None},
      },
    },
  },

  {
    "Inkflow (External fluid flow)",
    {
      // Buffer A (feedback ink)
      {kInkflowBufferA, true, 1, {stellar::render::ShaderToyChannelSource::BufferA,
                                 stellar::render::ShaderToyChannelSource::External3, // Flow
                                 stellar::render::ShaderToyChannelSource::External2, // Dye
                                 stellar::render::ShaderToyChannelSource::External1}}, // Noise
      // Buffer B
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Buffer C
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Buffer D
      {kBlackSnippet, false, 1, {stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None,
                                stellar::render::ShaderToyChannelSource::None}},
      // Image
      {kInkflowImage, true, 1, {stellar::render::ShaderToyChannelSource::BufferA,
                               stellar::render::ShaderToyChannelSource::None,
                               stellar::render::ShaderToyChannelSource::None,
                               stellar::render::ShaderToyChannelSource::None}},
    },
  },
};

static constexpr int kPresetCount = (int)(sizeof(kPresets) / sizeof(kPresets[0]));

// -----------------------------
// Graph / compilation helpers
// -----------------------------

static void syncGraphSettings(ProceduralShaderLabWindowState& st) {
  for (int p = 0; p < ProceduralShaderLabWindowState::kPassCount; ++p) {
    stellar::render::ShaderToyGraphPassSettings s;
    s.enabled = (p == 4) ? true : st.passes[p].enabled;
    s.resolutionScale = st.passes[p].resolutionScale;
    s.channels = st.passes[p].channels;
    st.graph.setPassSettings(passEnum(p), s);
  }
}

static bool updateParamsFromCode(ProceduralShaderLabWindowState& st, bool preserveExistingValues = true) {
  // Gather all pass sources.
  std::array<std::string_view, ProceduralShaderLabWindowState::kPassCount> src{};
  for (int p = 0; p < ProceduralShaderLabWindowState::kPassCount; ++p) {
    src[p] = cbufView(st.code[p]);
  }

  const stellar::render::ShaderToyParamSet* preserve = preserveExistingValues ? &st.params : nullptr;
  stellar::render::ShaderToyParamSet parsed = stellar::render::parseShaderToyParamsFromSources(src, preserve);

  const bool schemaChanged = !st.params.schemaEquals(parsed);
  st.params = std::move(parsed);

  if (st.graphInited) {
    st.graph.setUserHeader(st.params.buildUniformDecls());
  }

  return schemaChanged;
}

static void applyPreset(ProceduralShaderLabWindowState& st, int presetIdx, const ToastFn& toast) {
  presetIdx = std::clamp(presetIdx, 0, kPresetCount - 1);
  st.presetIndex = presetIdx;

  const GraphPreset& pr = kPresets[presetIdx];
  for (int p = 0; p < ProceduralShaderLabWindowState::kPassCount; ++p) {
    st.passes[p].enabled = pr.pass[p].enabled;
    st.passes[p].resolutionScale = pr.pass[p].scale;
    st.passes[p].channels = pr.pass[p].channels;
    st.passes[p].dirty = true;
    copyToBuf(st.code[p], pr.pass[p].code ? std::string_view(pr.pass[p].code) : std::string_view(kBlackSnippet));
  }

  st.editPassIndex = 4;
  st.previewPassIndex = 4;

  // New graph = reset buffers.
  st.requestResetBuffers = true;

  // Refresh the parameter list immediately so the UI updates even before a compile.
  (void)updateParamsFromCode(st, /*preserveExistingValues=*/false);

  if (st.autoCompileOnPreset) {
    st.requestCompileAll = true;
  }

  if (toast) toast(std::string("Applied preset: ") + pr.name, 1.6);
}

static bool ensureInit(ProceduralShaderLabWindowState& st) {
  if (!st.previewInited) {
    std::string err;
    if (!st.previewTarget.init(st.previewResolution, st.previewResolution, &err)) {
      st.previewInitError = err;
      return false;
    }
    st.previewInited = true;
  }

  if (!st.graphInited) {
    std::string err;
    if (!st.graph.init(&err)) {
      st.previewInitError = "ShaderToyGraph init failed: " + err;
      return false;
    }
    st.graphInited = true;
  }

  return true;
}

static void compilePassNow(ProceduralShaderLabWindowState& st, int passIdx, const ToastFn& toast) {
  if (!st.graphInited) return;

  // Keep the global parameter schema in sync. If it changed, we need to
  // recompile all passes so every pass gets the same injected uniform block.
  const bool schemaChanged = updateParamsFromCode(st, /*preserveExistingValues=*/true);
  if (schemaChanged) {
    // Delegate to compileAll to avoid subtle "some passes know about params" states.
    if (toast) toast("Parameters changed; recompiling all passes...", 1.6);
    // compileAllNow will call updateParamsFromCode again, but that's fine.
    st.requestCompileAll = true;
    return;
  }

  const auto p = passEnum(passIdx);
  std::string err;
  const bool ok = st.graph.buildPass(p, cbufView(st.code[passIdx]), &err);
  if (!ok) {
    if (toast) toast("Compile failed (see error).", 3.0);
  } else {
    st.passes[passIdx].dirty = false;
    if (toast) toast(std::string("Compiled ") + passLabel(passIdx) + ".", 1.3);
  }
}

static void compileAllNow(ProceduralShaderLabWindowState& st, const ToastFn& toast) {
  if (!st.graphInited) return;

  // Refresh parameter schema & header once for the whole graph.
  (void)updateParamsFromCode(st, /*preserveExistingValues=*/true);

  bool anyFail = false;
  for (int p = 0; p < ProceduralShaderLabWindowState::kPassCount; ++p) {
    std::string err;
    const bool ok = st.graph.buildPass(passEnum(p), cbufView(st.code[p]), &err);
    if (ok) {
      st.passes[p].dirty = false;
    } else {
      anyFail = true;
    }
  }

  if (toast) {
    if (anyFail) toast("Compile all: some passes failed (see errors).", 3.0);
    else toast("Compiled all passes.", 1.5);
  }
}

static void updateLiveReload(ProceduralShaderLabWindowState& st, const ToastFn& toast) {
  if (!st.liveReload) return;

  const fs::path p(st.filePath);
  const std::uint64_t stamp = fileStamp(p);

  if (stamp == 0) {
    if (!st.liveReloadMissingOk) {
      st.fileError = "Live reload: file missing.";
    }
    return;
  }

  if (st.liveReloadStamp == 0) {
    st.liveReloadStamp = stamp;
    return;
  }

  if (stamp == st.liveReloadStamp) return;
  st.liveReloadStamp = stamp;

  std::string text;
  std::string err;
  if (!readTextFile(p, text, err)) {
    st.fileError = "Live reload read failed: " + err;
    return;
  }

  if (!loadGraphText(st, text, err)) {
    st.fileError = "Live reload parse failed: " + err;
    return;
  }

  if (toast) toast("Live reloaded graph file.", 1.5);
}

// -----------------------------
// External iChannel textures (External0..3)
// -----------------------------

static void rebuildChecker(ProceduralShaderLabWindowState::ExternalTextures& e) {
  e.checkerRes = std::clamp(e.checkerRes, 64, 1024);
  e.checkerSize = std::clamp(e.checkerSize, 2, std::max(2, e.checkerRes / 4));
  e.checkerTex.createChecker(e.checkerRes, e.checkerRes, e.checkerSize);
}

static void rebuildNoise(ProceduralShaderLabWindowState::ExternalTextures& e) {
  e.noiseRes = std::clamp(e.noiseRes, 64, 1024);
  const int w = e.noiseRes;
  const int h = e.noiseRes;
  e.noiseRGBA.resize((size_t)w * (size_t)h * 4u);

  core::u64 rng = core::hashCombine(e.noiseSeed, core::fnv1a64("shaderLabNoise"));
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      const float r = rand01(rng);
      const float g = rand01(rng);
      const float b = rand01(rng);
      const size_t i = ((size_t)y * (size_t)w + (size_t)x) * 4u;
      e.noiseRGBA[i + 0] = (unsigned char)std::clamp((int)std::lround(r * 255.0f), 0, 255);
      e.noiseRGBA[i + 1] = (unsigned char)std::clamp((int)std::lround(g * 255.0f), 0, 255);
      e.noiseRGBA[i + 2] = (unsigned char)std::clamp((int)std::lround(b * 255.0f), 0, 255);
      e.noiseRGBA[i + 3] = 255;
    }
  }

  if (!e.noiseTex.handle() || e.noiseTex.width() != w || e.noiseTex.height() != h) {
    // Repeatable tiling noise.
    e.noiseTex.createRGBA(w, h, e.noiseRGBA.data(), /*generateMips=*/true, /*nearest=*/false, /*clampToEdge=*/false);
  } else {
    e.noiseTex.updateRGBA(0, 0, w, h, e.noiseRGBA.data());
  }
  e.requestRegenNoise = false;
}

static void seedFluid(stellar::proc::FluidSim2D& sim, core::u64 seed) {
  sim.clear();

  core::u64 rng = core::hashCombine(seed, core::fnv1a64("shaderLabFluidSeed"));
  for (int i = 0; i < 10; ++i) {
    const float x = 0.15f + 0.70f * rand01(rng);
    const float y = 0.15f + 0.70f * rand01(rng);
    const float radius = 0.03f + 0.05f * rand01(rng);

    const float ang = 6.2831853f * rand01(rng);
    const float spd = 30.0f + 80.0f * rand01(rng);
    const float vx = std::cos(ang) * spd;
    const float vy = std::sin(ang) * spd;

    const float r = 0.5f + 1.5f * rand01(rng);
    const float g = 0.5f + 1.5f * rand01(rng);
    const float b = 0.5f + 1.5f * rand01(rng);

    sim.splat(x, y, radius, vx, vy, r * 12.0f, g * 12.0f, b * 12.0f);
  }
}

static void ensureFluid(ProceduralShaderLabWindowState::ExternalTextures& e) {
  if (e.fluidGrid < 32) e.fluidGrid = 32;
  if (e.fluidGrid > 512) e.fluidGrid = 512;

  if (e.fluid.gridSize() != e.fluidGrid) {
    e.fluid.resize(e.fluidGrid);

    // Pleasant defaults for a "living texture".
    auto& p = e.fluid.params();
    p.viscosity = 0.00015f;
    p.diffusion = 0.00005f;
    p.dyeDissipation = 0.08f;
    p.velocityDamping = 0.0f;
    p.vorticityConfinement = 22.0f;
    p.curlNoiseStrength = 7.0f;
    p.curlNoiseFrequency = 2.2f;
    p.curlNoiseTimeScale = 0.20f;
    p.maxSpeed = 220.0f;
    p.maxDye = 80.0f;

    seedFluid(e.fluid, e.fluidSeed);
  }

  e.flowRes = e.fluid.gridSize();

  // Allocate GPU textures on first use (or if resized).
  const int n = e.fluid.gridSize();
  if (!e.fluidTex.handle() || e.fluidTex.width() != n || e.fluidTex.height() != n) {
    e.fluidTex.allocateRGBA(n, n, /*generateMips=*/true, /*nearest=*/false, /*clampToEdge=*/true);
  }
  if (!e.flowTex.handle() || e.flowTex.width() != n || e.flowTex.height() != n) {
    e.flowTex.allocateRGBA(n, n, /*generateMips=*/true, /*nearest=*/false, /*clampToEdge=*/true);
  }
}

static inline float saturate(float x) { return std::clamp(x, 0.0f, 1.0f); }

static void updateFluidTextures(ProceduralShaderLabWindowState::ExternalTextures& e) {
  const int n = e.fluid.gridSize();
  if (n <= 0) return;

  e.fluidRGBA.resize((size_t)n * (size_t)n * 4u);
  e.flowRGBA.resize((size_t)n * (size_t)n * 4u);

  const auto& u = e.fluid.u();
  const auto& v = e.fluid.v();
  const auto& r = e.fluid.dyeR();
  const auto& g = e.fluid.dyeG();
  const auto& b = e.fluid.dyeB();
  const int stride = e.fluid.paddedStride();

  auto tonemap = [](float x, float exposure) {
    // Simple filmic-ish curve for pleasant display.
    return 1.0f - std::exp(-exposure * std::max(0.0f, x));
  };

  for (int y = 0; y < n; ++y) {
    for (int x = 0; x < n; ++x) {
      const int ix = x + 1;
      const int iy = y + 1;
      const int id = ix + stride * iy;

      const size_t o = ((size_t)y * (size_t)n + (size_t)x) * 4u;

      // Dye (RGB) -> External2
      const float rr = tonemap(r[id], e.fluidExposure);
      const float gg = tonemap(g[id], e.fluidExposure);
      const float bb = tonemap(b[id], e.fluidExposure);
      e.fluidRGBA[o + 0] = (unsigned char)std::clamp((int)std::lround(saturate(rr) * 255.0f), 0, 255);
      e.fluidRGBA[o + 1] = (unsigned char)std::clamp((int)std::lround(saturate(gg) * 255.0f), 0, 255);
      e.fluidRGBA[o + 2] = (unsigned char)std::clamp((int)std::lround(saturate(bb) * 255.0f), 0, 255);
      e.fluidRGBA[o + 3] = 255;

      // Flow (velocity) -> External3
      const float vx = u[id] * e.flowVizScale;
      const float vy = v[id] * e.flowVizScale;
      const float sp = std::sqrt(vx * vx + vy * vy);

      const float fr = saturate(0.5f + 0.5f * vx);
      const float fg = saturate(0.5f + 0.5f * vy);
      const float fb = saturate(0.25f + 0.75f * sp);

      e.flowRGBA[o + 0] = (unsigned char)std::clamp((int)std::lround(fr * 255.0f), 0, 255);
      e.flowRGBA[o + 1] = (unsigned char)std::clamp((int)std::lround(fg * 255.0f), 0, 255);
      e.flowRGBA[o + 2] = (unsigned char)std::clamp((int)std::lround(fb * 255.0f), 0, 255);
      e.flowRGBA[o + 3] = 255;
    }
  }

  e.fluidTex.updateRGBA(0, 0, n, n, e.fluidRGBA.data());
  e.flowTex.updateRGBA(0, 0, n, n, e.flowRGBA.data());
}

static void ensureExternalTextures(ProceduralShaderLabWindowState& st) {
  auto& e = st.ext;
  if (!e.inited) {
    rebuildChecker(e);
    rebuildNoise(e);
    ensureFluid(e);
    updateFluidTextures(e);
    e.inited = true;
  }

  // Lazy rebuilds.
  if (e.requestRegenNoise) {
    rebuildNoise(e);
  }

  // Keep graph bindings current.
  st.graph.setExternalTexture(0, &e.checkerTex);
  st.graph.setExternalTexture(1, &e.noiseTex);
  st.graph.setExternalTexture(2, e.fluidEnabled ? &e.fluidTex : nullptr);
  st.graph.setExternalTexture(3, e.fluidEnabled ? &e.flowTex : nullptr);
}

// -----------------------------
// Rendering
// -----------------------------

static void renderPreview(ProceduralShaderLabWindowState& st, float timeRealSec) {
  if (!st.previewInited || !st.graphInited) return;

  // Save GL state we touch.
  GLint prevFbo = 0;
  GLint prevViewport[4] = {0, 0, 0, 0};
  GLint prevProg = 0;
  GLint prevVao = 0;
  glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
  glGetIntegerv(GL_VIEWPORT, prevViewport);
  glGetIntegerv(GL_CURRENT_PROGRAM, &prevProg);
  glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &prevVao);

  const bool prevDepth = glIsEnabled(GL_DEPTH_TEST);
  const bool prevCull  = glIsEnabled(GL_CULL_FACE);
  const bool prevBlend = glIsEnabled(GL_BLEND);
  const bool prevSciss = glIsEnabled(GL_SCISSOR_TEST);

  // Make sure preview RT matches requested size.
  st.previewTarget.ensureSize(st.previewResolution, st.previewResolution);

  // Push pass settings into the graph.
  syncGraphSettings(st);

  // Ensure internal buffer RTs exist.
  {
    std::string err;
    if (!st.graph.ensureSize(st.previewTarget.width(), st.previewTarget.height(), &err)) {
      st.previewInitError = "ShaderToyGraph ensureSize failed: " + err;
      return;
    }
  }

  // Reset buffers when requested.
  if (st.requestResetBuffers) {
    st.graph.resetBuffers();
    st.frame = 0;
    st.requestResetBuffers = false;
  }

  // Time integration.
  if (!st.timeInit) {
    st.timeInit = true;
    st.lastRealTimeSec = timeRealSec;
    st.timeAccumSec = 0.0f;
  }

  const bool advance = (!st.paused) || st.stepFrame;

  float dt = 0.0f;
  float t = 0.0f;

  if (st.useRealTime) {
    float realDt = timeRealSec - st.lastRealTimeSec;
    if (realDt < 0.0f) realDt = 0.0f;
    st.lastRealTimeSec = timeRealSec;

    if (advance) {
      dt = st.stepFrame ? st.fixedStepSec : (realDt * st.timeScale);
      st.timeAccumSec += dt;
    }

    t = st.timeAccumSec;
  } else {
    if (advance) {
      dt = st.fixedStepSec;
      st.timeOverride += dt;
    }
    t = st.timeOverride;
  }

  // Update external iChannel textures before running the graph so passes can
  // reference them this frame.
  ensureExternalTextures(st);
  if (st.ext.fluidEnabled) {
    // Keep fluid buffers sized if the user changed settings.
    ensureFluid(st.ext);

    // Optional auto-injection to keep the texture alive.
    if (advance && st.ext.fluidAutoInject && st.ext.fluidInjectEveryNFrames > 0) {
      if ((st.frame % st.ext.fluidInjectEveryNFrames) == 0) {
        core::u64 rng = core::hashCombine(st.ext.fluidSeed,
          core::hashCombine((core::u64)st.frame, core::fnv1a64("shaderLabAutoInject")));

        const float x = 0.10f + 0.80f * rand01(rng);
        const float y = 0.10f + 0.80f * rand01(rng);

        const float ang = 6.2831853f * rand01(rng);
        const float spd = st.ext.fluidInjectForce * (0.4f + 0.6f * rand01(rng));
        const float vx = std::cos(ang) * spd;
        const float vy = std::sin(ang) * spd;

        const float r = 0.4f + 1.2f * rand01(rng);
        const float g = 0.4f + 1.2f * rand01(rng);
        const float b = 0.4f + 1.2f * rand01(rng);

        st.ext.fluid.splat(x, y,
                           st.ext.fluidInjectRadius,
                           vx, vy,
                           r * st.ext.fluidInjectDye,
                           g * st.ext.fluidInjectDye,
                           b * st.ext.fluidInjectDye);
      }
    }

    // Mouse injection (only when the preview is hovered/clicked).
    if (st.ext.fluidMouseInject && st.mouseDown) {
      const float u = st.mouseU;
      const float v = st.mouseV;
      float dx = 0.0f;
      float dy = 0.0f;
      if (st.ext.fluidWasDown) {
        dx = (u - st.ext.fluidLastU);
        dy = (v - st.ext.fluidLastV);
      }
      st.ext.fluidLastU = u;
      st.ext.fluidLastV = v;
      st.ext.fluidWasDown = true;

      // Map drag delta to a strong local impulse.
      const float vx = dx * st.ext.fluidInjectForce * 200.0f;
      const float vy = dy * st.ext.fluidInjectForce * 200.0f;

      // Dye color: warm ramp that reads well under tonemapping.
      st.ext.fluid.splat(u, v,
                         st.ext.fluidInjectRadius,
                         vx, vy,
                         1.2f * st.ext.fluidInjectDye,
                         0.6f * st.ext.fluidInjectDye,
                         0.25f * st.ext.fluidInjectDye);
    } else {
      st.ext.fluidWasDown = false;
    }

    // Advance sim.
    if (advance && !st.ext.fluidPaused) {
      const float dts = std::clamp(dt * st.ext.fluidTimeScale, 0.0f, st.ext.fluidMaxDt);
      if (dts > 0.0f) {
        st.ext.fluid.step(dts, st.ext.fluidIterations, st.ext.fluidSeed, t);
      }
    }

    // Always refresh textures if something could have written to the sim.
    updateFluidTextures(st.ext);
  }

  // Orbit camera around origin.
  float yaw = (float)stellar::math::degToRad((double)st.yawDeg);
  if (st.autoOrbit) yaw += (float)stellar::math::degToRad((double)(st.orbitDegPerSec * t));
  const float pitch = (float)stellar::math::degToRad((double)st.pitchDeg);

  const float cy = std::cos(yaw);
  const float sy = std::sin(yaw);
  const float cp = std::cos(pitch);
  const float sp = std::sin(pitch);

  const stellar::math::Vec3d eye{(double)(st.distance * cp * cy),
                                (double)(st.distance * sp),
                                (double)(st.distance * cp * sy)};

  stellar::math::Vec3d fwd = (-eye).normalized();
  const stellar::math::Vec3d worldUp{0, 1, 0};
  stellar::math::Vec3d right = stellar::math::cross(fwd, worldUp).normalized();
  if (right.lengthSq() < 1e-10) right = stellar::math::Vec3d{1, 0, 0};
  stellar::math::Vec3d up = stellar::math::cross(right, fwd).normalized();

  const float fovY = (float)stellar::math::degToRad((double)st.fovYDeg);
  const float tanHalfFovY = std::tan(0.5f * fovY);

  stellar::render::ShaderToyInputs in{};
  in.width = st.previewTarget.width();
  in.height = st.previewTarget.height();
  in.timeSec = t;
  in.timeDeltaSec = advance ? dt : 0.0f;
  in.frame = st.frame;
  in.seed = seed32(st.seed);
  in.userParams = &st.params;

  // iMouse: pixels, origin lower-left.
  const float mx = st.mouseU * (float)in.width;
  const float my = st.mouseV * (float)in.height;
  in.mouse[0] = mx;
  in.mouse[1] = my;
  in.mouse[2] = st.mouseDown ? (st.mouseDownU * (float)in.width) : 0.0f;
  in.mouse[3] = st.mouseDown ? (st.mouseDownV * (float)in.height) : 0.0f;

  in.camPos[0] = (float)eye.x;
  in.camPos[1] = (float)eye.y;
  in.camPos[2] = (float)eye.z;

  in.camRight[0] = (float)right.x;
  in.camRight[1] = (float)right.y;
  in.camRight[2] = (float)right.z;

  in.camUp[0] = (float)up.x;
  in.camUp[1] = (float)up.y;
  in.camUp[2] = (float)up.z;

  in.camForward[0] = (float)fwd.x;
  in.camForward[1] = (float)fwd.y;
  in.camForward[2] = (float)fwd.z;

  in.tanHalfFovY = tanHalfFovY;

  // Render graph.
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);
  glDisable(GL_BLEND);
  glDisable(GL_SCISSOR_TEST);

  st.graph.render(in, st.previewTarget);

  // Advance frame counter after rendering.
  if (advance) {
    st.frame += 1;
  }
  if (st.stepFrame) {
    st.stepFrame = false;
  }

  // Restore state.
  if (prevDepth) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
  if (prevCull) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);
  if (prevBlend) glEnable(GL_BLEND); else glDisable(GL_BLEND);
  if (prevSciss) glEnable(GL_SCISSOR_TEST); else glDisable(GL_SCISSOR_TEST);

  glUseProgram((GLuint)prevProg);
  glBindVertexArray((GLuint)prevVao);
  glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
  glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);
}

static void updatePreviewMouse(ProceduralShaderLabWindowState& st) {
  // Use the last drawn ImGui::Image rectangle.
  if (!ImGui::IsItemHovered()) {
    if (!ImGui::IsMouseDown(0)) st.mouseDown = false;
    return;
  }

  const ImVec2 min = ImGui::GetItemRectMin();
  const ImVec2 max = ImGui::GetItemRectMax();
  const ImVec2 mp = ImGui::GetMousePos();

  const float w = std::max(1.0f, max.x - min.x);
  const float h = std::max(1.0f, max.y - min.y);

  float u = (mp.x - min.x) / w;
  float v = 1.0f - (mp.y - min.y) / h;
  u = std::clamp(u, 0.0f, 1.0f);
  v = std::clamp(v, 0.0f, 1.0f);

  st.mouseU = u;
  st.mouseV = v;

  if (ImGui::IsMouseClicked(0)) {
    st.mouseDown = true;
    st.mouseDownU = u;
    st.mouseDownV = v;
  }

  if (!ImGui::IsMouseDown(0)) {
    st.mouseDown = false;
  }
}

} // namespace

void drawProceduralShaderLabWindow(ProceduralShaderLabWindowState& st, float timeSec, const ToastFn& toast) {
  if (!st.open) return;

  // Initialize default content.
  if (st.code[4][0] == '\0') {
    applyPreset(st, 0, nullptr);
  }

  // Lazy init.
  (void)ensureInit(st);

  updateLiveReload(st, toast);

  ImGui::SetNextWindowSize(ImVec2(1200, 820), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Procedural Shader Lab", &st.open)) {
    ImGui::End();
    return;
  }

  if (!ensureInit(st)) {
    ImGui::TextWrapped("Preview init failed: %s", st.previewInitError.c_str());
    ImGui::End();
    return;
  }

  // Handle queued actions from last frame.
  if (st.requestCompileAll) {
    compileAllNow(st, toast);
    st.requestCompileAll = false;
  }

  // Render preview before UI draws it.
  renderPreview(st, timeSec);

  if (ImGui::BeginTable("##shaderlab", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
    ImGui::TableSetupColumn("Controls", ImGuiTableColumnFlags_WidthStretch, 0.0f);
    ImGui::TableSetupColumn("Preview", ImGuiTableColumnFlags_WidthFixed, 500.0f);

    // ---------------- Controls ----------------
    ImGui::TableNextColumn();

    // Preset picker
    {
      const char* preview = kPresets[std::clamp(st.presetIndex, 0, kPresetCount - 1)].name;
      if (ImGui::BeginCombo("Preset", preview)) {
        for (int i = 0; i < kPresetCount; ++i) {
          const bool selected = (i == st.presetIndex);
          if (ImGui::Selectable(kPresets[i].name, selected)) {
            applyPreset(st, i, toast);
          }
          if (selected) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
      }
      ImGui::SameLine();
      if (ImGui::Button("Re-apply")) {
        applyPreset(st, st.presetIndex, toast);
      }
      ImGui::SameLine();
      ImGui::Checkbox("Auto-compile preset", &st.autoCompileOnPreset);
    }

    // External iChannel textures (External0..3)
    {
      ImGui::Separator();
      if (ImGui::CollapsingHeader("External Textures (iChannel)", ImGuiTreeNodeFlags_DefaultOpen)) {
        // Ensure textures exist for previews even if the graph is paused.
        ensureExternalTextures(st);

        ImGui::TextWrapped(
          "Route these via iChannel0..3 using the channel selectors: "
          "Checker (Ext0), Noise (Ext1), Fluid (Ext2), Flow (Ext3).  "
          "Flow encodes velocity in RG with 0.5=center (no motion)."
        );

        const float thumb = 86.0f;

        // Thumbnails
        if (ImGui::BeginTable("##ext_tex_thumbs", 2, ImGuiTableFlags_SizingFixedFit)) {
          auto thumbCell = [&](const char* label, const stellar::render::Texture2D& tex) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted(label);
            if (tex.handle()) {
              ImGui::Image((ImTextureID)(intptr_t)tex.handle(), ImVec2(thumb, thumb));
            } else {
              ImGui::Dummy(ImVec2(thumb, thumb));
            }
          };

          ImGui::TableNextRow();
          thumbCell("Checker (Ext0)", st.ext.checkerTex);
          thumbCell("Noise   (Ext1)", st.ext.noiseTex);

          ImGui::TableNextRow();
          if (st.ext.fluidEnabled) {
            thumbCell("Fluid   (Ext2)", st.ext.fluidTex);
            thumbCell("Flow    (Ext3)", st.ext.flowTex);
          } else {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Fluid   (Ext2)");
            ImGui::Dummy(ImVec2(thumb, thumb));
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("Flow    (Ext3)");
            ImGui::Dummy(ImVec2(thumb, thumb));
          }

          ImGui::EndTable();
        }

        ImGui::Spacing();

        // Checker controls
        {
          bool changed = false;
          changed |= ImGui::SliderInt("Checker size", &st.ext.checkerSize, 2, 64);
          changed |= ImGui::SliderInt("Checker res", &st.ext.checkerRes, 64, 1024);
          if (changed) {
            rebuildChecker(st.ext);
          }
        }

        // Noise controls
        {
          ImGui::SeparatorText("Noise (Ext1)");
          if (ImGui::SliderInt("Noise res", &st.ext.noiseRes, 64, 1024)) {
            st.ext.requestRegenNoise = true;
          }

          // Use a 64-bit seed input via two 32-bit halves for portability.
          std::uint32_t seedLo = (std::uint32_t)(st.ext.noiseSeed & 0xFFFFFFFFULL);
          std::uint32_t seedHi = (std::uint32_t)((st.ext.noiseSeed >> 32) & 0xFFFFFFFFULL);
          bool seedChanged = false;
          seedChanged |= ImGui::InputScalar("Noise seed (lo)", ImGuiDataType_U32, &seedLo);
          seedChanged |= ImGui::InputScalar("Noise seed (hi)", ImGuiDataType_U32, &seedHi);
          if (seedChanged) {
            st.ext.noiseSeed = ((core::u64)seedHi << 32) | (core::u64)seedLo;
            st.ext.requestRegenNoise = true;
          }

          if (ImGui::Button("Regenerate noise")) {
            st.ext.requestRegenNoise = true;
          }
          ImGui::SameLine();
          if (ImGui::Button("Randomize seed")) {
            core::u64 rng = core::hashCombine(st.ext.noiseSeed, core::hashCombine((core::u64)st.frame, core::fnv1a64("noiseSeedUI")));
            st.ext.noiseSeed = splitmix64(rng);
            st.ext.requestRegenNoise = true;
          }

          if (st.ext.requestRegenNoise) {
            rebuildNoise(st.ext);
          }
        }

        // Fluid controls
        {
          ImGui::SeparatorText("Fluid (Ext2) + Flow (Ext3)");
          ImGui::Checkbox("Enable fluid externals", &st.ext.fluidEnabled);
          ImGui::SameLine();
          ImGui::Checkbox("Pause fluid", &st.ext.fluidPaused);
          ImGui::SameLine();
          ImGui::Checkbox("Mouse inject", &st.ext.fluidMouseInject);

          if (ImGui::SliderInt("Grid", &st.ext.fluidGrid, 32, 512)) {
            ensureFluid(st.ext);
            updateFluidTextures(st.ext);
          }
          ImGui::SliderInt("Iterations", &st.ext.fluidIterations, 4, 60);
          ImGui::SliderFloat("Time scale", &st.ext.fluidTimeScale, 0.0f, 3.0f, "%.2f");
          ImGui::SliderFloat("Max dt", &st.ext.fluidMaxDt, 1.0f / 240.0f, 1.0f / 10.0f, "%.4f");
          ImGui::SliderFloat("Exposure", &st.ext.fluidExposure, 0.005f, 0.25f, "%.3f");
          ImGui::SliderFloat("Flow viz scale", &st.ext.flowVizScale, 0.001f, 0.08f, "%.4f");

          ImGui::Checkbox("Auto inject", &st.ext.fluidAutoInject);
          ImGui::SameLine();
          ImGui::SliderInt("Every N frames", &st.ext.fluidInjectEveryNFrames, 1, 240);
          ImGui::SliderFloat("Inject radius", &st.ext.fluidInjectRadius, 0.005f, 0.20f, "%.3f");
          ImGui::SliderFloat("Inject dye", &st.ext.fluidInjectDye, 0.0f, 40.0f, "%.1f");
          ImGui::SliderFloat("Inject force", &st.ext.fluidInjectForce, 0.0f, 200.0f, "%.1f");

          // Seed controls
          std::uint32_t fSeedLo = (std::uint32_t)(st.ext.fluidSeed & 0xFFFFFFFFULL);
          std::uint32_t fSeedHi = (std::uint32_t)((st.ext.fluidSeed >> 32) & 0xFFFFFFFFULL);
          bool fSeedChanged = false;
          fSeedChanged |= ImGui::InputScalar("Fluid seed (lo)", ImGuiDataType_U32, &fSeedLo);
          fSeedChanged |= ImGui::InputScalar("Fluid seed (hi)", ImGuiDataType_U32, &fSeedHi);
          if (fSeedChanged) {
            st.ext.fluidSeed = ((core::u64)fSeedHi << 32) | (core::u64)fSeedLo;
          }

          if (ImGui::Button("Reseed fluid")) {
            ensureFluid(st.ext);
            seedFluid(st.ext.fluid, st.ext.fluidSeed);
            updateFluidTextures(st.ext);
          }
          ImGui::SameLine();
          if (ImGui::Button("Clear")) {
            ensureFluid(st.ext);
            st.ext.fluid.clear();
            updateFluidTextures(st.ext);
          }
        }
      }
    }

    // Pass tabs
    {
      if (ImGui::BeginTabBar("##pass_tabs")) {
        for (int p = 0; p < ProceduralShaderLabWindowState::kPassCount; ++p) {
          const bool selected = (st.editPassIndex == p);
          ImGuiTabItemFlags flags = 0;
          if (st.passes[p].dirty) flags |= ImGuiTabItemFlags_UnsavedDocument;

          if (ImGui::BeginTabItem(passLabel(p), nullptr, flags)) {
            st.editPassIndex = p;
            ImGui::EndTabItem();
          } else if (selected) {
            // Keep selection stable.
          }
        }
        ImGui::EndTabBar();
      }
    }

    // Top action buttons
    {
      if (ImGui::Button("Compile Pass")) {
        compilePassNow(st, st.editPassIndex, toast);
      }
      ImGui::SameLine();
      if (ImGui::Button("Compile All")) {
        compileAllNow(st, toast);
      }
      ImGui::SameLine();
      if (ImGui::Button("Reset Buffers")) {
        st.requestResetBuffers = true;
      }
      ImGui::SameLine();
      ImGui::Checkbox("Show generated GLSL", &st.showGeneratedSource);
    }

    // File I/O
    {
      ImGui::Separator();
      ImGui::TextUnformatted("File");
      ImGui::InputText("Path", st.filePath, (int)sizeof(st.filePath));
      ImGui::SameLine();
      if (ImGui::Button("Load")) {
        std::string text;
        std::string err;
        if (!readTextFile(fs::path(st.filePath), text, err)) {
          st.fileError = "Load failed: " + err;
        } else {
          if (!loadGraphText(st, text, err)) {
            st.fileError = "Parse failed: " + err;
          } else {
            st.fileError.clear();
            if (toast) toast("Loaded graph.", 1.4);
          }
        }
      }
      ImGui::SameLine();
      if (ImGui::Button("Save")) {
        // Ensure parameter overrides are up to date in the saved graph.
        (void)updateParamsFromCode(st, /*preserveExistingValues=*/true);
        std::string err;
        const std::string text = saveGraphText(st);
        if (!writeTextFile(fs::path(st.filePath), text, err)) {
          st.fileError = "Save failed: " + err;
        } else {
          st.fileError.clear();
          if (toast) toast("Saved graph.", 1.4);
        }
      }
      ImGui::SameLine();
      ImGui::Checkbox("Live reload", &st.liveReload);
      ImGui::SameLine();
      ImGui::Checkbox("Missing OK", &st.liveReloadMissingOk);

      if (!st.fileError.empty()) {
        ImGui::TextColored(ImVec4(1.0f, 0.65f, 0.35f, 1.0f), "%s", st.fileError.c_str());
      }
    }

    // Time controls
    {
      ImGui::Separator();
      ImGui::TextUnformatted("Time");
      ImGui::Checkbox("Use real time", &st.useRealTime);
      ImGui::SameLine();
      ImGui::Checkbox("Paused", &st.paused);
      ImGui::SameLine();
      if (ImGui::Button("Step")) {
        st.stepFrame = true;
      }

      if (st.useRealTime) {
        ImGui::SliderFloat("Time scale", &st.timeScale, -4.0f, 4.0f, "%.2fx");
      } else {
        ImGui::SliderFloat("Time", &st.timeOverride, 0.0f, 120.0f, "%.2f s");
      }
      ImGui::SliderFloat("Fixed step", &st.fixedStepSec, 1.0f / 240.0f, 1.0f / 10.0f, "%.4f s");
      ImGui::Text("iFrame: %d", st.frame);
    }

    // User parameters (parsed from shader code)
    {
      ImGui::Separator();
      ImGui::SetNextItemOpen(st.paramsOpen, ImGuiCond_FirstUseEver);
      const bool open = ImGui::CollapsingHeader("Parameters");
      st.paramsOpen = open;
      if (open) {
        ImGui::InputText("Filter", st.paramsFilter, (int)sizeof(st.paramsFilter));
        ImGui::SameLine();
        ImGui::Checkbox("Advanced", &st.paramsShowAdvanced);
        ImGui::SameLine();
        if (ImGui::Button("Reset to defaults##params")) {
          st.params.resetToDefaults();
        }

        if (st.params.empty()) {
          ImGui::TextWrapped("(No parameters parsed yet. Add `// @param ...` directives and compile.)");
        } else {
          ImGui::TextDisabled("%d parameters", (int)st.params.defs.size());

          const std::string_view needle = st.paramsFilter;

          // Group parameters by the first-seen order of `@group`.
          std::vector<std::string> groupOrder;
          std::unordered_map<std::string, int> groupToIdx;
          std::vector<std::vector<int>> buckets;
          groupOrder.reserve(8);
          buckets.reserve(8);

          for (int i = 0; i < (int)st.params.defs.size(); ++i) {
            const std::string& g = st.params.defs[i].group;
            auto it = groupToIdx.find(g);
            if (it == groupToIdx.end()) {
              const int idx = (int)buckets.size();
              groupToIdx.emplace(g, idx);
              groupOrder.push_back(g);
              buckets.push_back({});
              buckets.back().reserve(16);
              it = groupToIdx.find(g);
            }
            buckets[it->second].push_back(i);
          }

          auto drawParam = [&](int i) {
            const auto& d = st.params.defs[i];
            auto& v = st.params.values[i];

            // Filter by name/label.
            if (!needle.empty() && !icontains(d.name, needle) && !icontains(d.label, needle)) return;

            // Unique ImGui ID per param.
            const std::string uiLabel = d.label + "##" + d.name;

            bool changed = false;
            switch (d.type) {
              case stellar::render::ShaderToyParamType::Float: {
                float x = v[0];
                changed = ImGui::DragFloat(uiLabel.c_str(), &x, d.step, d.minValue[0], d.maxValue[0], "%.4f");
                if (changed) v[0] = std::clamp(x, d.minValue[0], d.maxValue[0]);
                break;
              }
              case stellar::render::ShaderToyParamType::Int: {
                int x = (int)std::lround(v[0]);
                changed = ImGui::DragInt(uiLabel.c_str(), &x, (float)std::max(1.0f, d.step), (int)d.minValue[0], (int)d.maxValue[0]);
                if (changed) v[0] = (float)std::clamp(x, (int)d.minValue[0], (int)d.maxValue[0]);
                break;
              }
              case stellar::render::ShaderToyParamType::Bool: {
                bool b = (v[0] >= 0.5f);
                changed = ImGui::Checkbox(uiLabel.c_str(), &b);
                if (changed) v[0] = b ? 1.0f : 0.0f;
                break;
              }
              case stellar::render::ShaderToyParamType::Vec2: {
                float a[2] = {v[0], v[1]};
                changed = ImGui::DragFloat2(uiLabel.c_str(), a, d.step, d.minValue[0], d.maxValue[0], "%.4f");
                if (changed) {
                  v[0] = std::clamp(a[0], d.minValue[0], d.maxValue[0]);
                  v[1] = std::clamp(a[1], d.minValue[1], d.maxValue[1]);
                }
                break;
              }
              case stellar::render::ShaderToyParamType::Vec3: {
                float a[3] = {v[0], v[1], v[2]};
                changed = ImGui::DragFloat3(uiLabel.c_str(), a, d.step, d.minValue[0], d.maxValue[0], "%.4f");
                if (changed) {
                  v[0] = std::clamp(a[0], d.minValue[0], d.maxValue[0]);
                  v[1] = std::clamp(a[1], d.minValue[1], d.maxValue[1]);
                  v[2] = std::clamp(a[2], d.minValue[2], d.maxValue[2]);
                }
                break;
              }
              case stellar::render::ShaderToyParamType::Color3: {
                float a[3] = {v[0], v[1], v[2]};
                changed = ImGui::ColorEdit3(uiLabel.c_str(), a);
                if (changed) {
                  v[0] = std::clamp(a[0], 0.0f, 1.0f);
                  v[1] = std::clamp(a[1], 0.0f, 1.0f);
                  v[2] = std::clamp(a[2], 0.0f, 1.0f);
                }
                break;
              }
            }

            if (st.paramsShowAdvanced) {
              ImGui::SameLine();
              ImGui::TextDisabled("[%s %s]", stellar::render::shaderToyParamGlslType(d.type), d.name.c_str());
            }
          };

          for (int gi = 0; gi < (int)groupOrder.size(); ++gi) {
            const std::string& g = groupOrder[gi];
            const char* header = g.empty() ? "Ungrouped" : g.c_str();
            ImGui::SeparatorText(header);
            for (int idx : buckets[gi]) {
              drawParam(idx);
            }
          }

          if (st.paramsShowAdvanced) {
            ImGui::Spacing();
            if (ImGui::Button("Copy uniform block##params")) {
              const std::string u = st.params.buildUniformDecls();
              ImGui::SetClipboardText(u.c_str());
              if (toast) toast("Copied uniform block to clipboard.", 1.2);
            }
          }
        }
      }
    }

    // Camera controls
    {
      ImGui::Separator();
      ImGui::TextUnformatted("Camera (for raymarch sketches)");
      ImGui::SliderFloat("Yaw", &st.yawDeg, -180.0f, 180.0f, "%.1f deg");
      ImGui::SliderFloat("Pitch", &st.pitchDeg, -89.0f, 89.0f, "%.1f deg");
      ImGui::SliderFloat("Distance", &st.distance, 0.25f, 20.0f, "%.2f");
      ImGui::SliderFloat("FovY", &st.fovYDeg, 15.0f, 120.0f, "%.0f deg");
      ImGui::Checkbox("Auto orbit", &st.autoOrbit);
      ImGui::SliderFloat("Orbit speed", &st.orbitDegPerSec, -60.0f, 60.0f, "%.1f deg/s");
    }

    // Pass-specific settings
    {
      ImGui::Separator();
      const int p = st.editPassIndex;
      ImGui::Text("Editing: %s", passLabel(p));

      if (p < 4) {
        ImGui::Checkbox("Enabled", &st.passes[p].enabled);
        ImGui::SameLine();
      }

      // Resolution scale.
      {
        int scaleIdx = 0;
        const int s = st.passes[p].resolutionScale;
        if (s == 1) scaleIdx = 0;
        else if (s == 2) scaleIdx = 1;
        else scaleIdx = 2;

        const char* items = "Full (1x)\0Half (1/2)\0Quarter (1/4)\0";
        if (ImGui::Combo("Resolution", &scaleIdx, items)) {
          st.passes[p].resolutionScale = (scaleIdx == 0) ? 1 : ((scaleIdx == 1) ? 2 : 4);
          st.requestResetBuffers = true; // changing resolution implies re-init.
        }
      }

      // Channels.
      {
        const char* items =
          "None\0"
          "Buffer A\0"
          "Buffer B\0"
          "Buffer C\0"
          "Buffer D\0"
          "Checker (Ext0)\0"
          "Noise   (Ext1)\0"
          "Fluid   (Ext2)\0"
          "Flow    (Ext3)\0";
        for (int c = 0; c < 4; ++c) {
          int idx = 0;
          switch (st.passes[p].channels[c]) {
            case stellar::render::ShaderToyChannelSource::None: idx = 0; break;
            case stellar::render::ShaderToyChannelSource::BufferA: idx = 1; break;
            case stellar::render::ShaderToyChannelSource::BufferB: idx = 2; break;
            case stellar::render::ShaderToyChannelSource::BufferC: idx = 3; break;
            case stellar::render::ShaderToyChannelSource::BufferD: idx = 4; break;
            case stellar::render::ShaderToyChannelSource::External0: idx = 5; break;
            case stellar::render::ShaderToyChannelSource::External1: idx = 6; break;
            case stellar::render::ShaderToyChannelSource::External2: idx = 7; break;
            case stellar::render::ShaderToyChannelSource::External3: idx = 8; break;
          }

          char label[32];
          std::snprintf(label, sizeof(label), "iChannel%d", c);
          if (ImGui::Combo(label, &idx, items)) {
            using S = stellar::render::ShaderToyChannelSource;
            st.passes[p].channels[c] =
              (idx == 0) ? S::None :
              (idx == 1) ? S::BufferA :
              (idx == 2) ? S::BufferB :
              (idx == 3) ? S::BufferC :
              (idx == 4) ? S::BufferD :
              (idx == 5) ? S::External0 :
              (idx == 6) ? S::External1 :
              (idx == 7) ? S::External2 :
              S::External3;
          }
        }
      }

      // Error display for this pass.
      const auto passE = passEnum(p);
      if (!st.graph.passValid(passE) && !st.graph.passError(passE).empty()) {
        ImGui::Spacing();
        ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "Compile error:");
        ImGui::BeginChild("##pass_err", ImVec2(0, 120), true);
        ImGui::TextUnformatted(st.graph.passError(passE).c_str());
        ImGui::EndChild();
      }

      // Code editor.
      ImGui::Spacing();
      ImGui::TextUnformatted("Code");
      ImGuiInputTextFlags flags = ImGuiInputTextFlags_AllowTabInput;
      if (ImGui::InputTextMultiline("##code", st.code[p].data(), (int)ProceduralShaderLabWindowState::kCodeCapacity, ImVec2(0, 260), flags)) {
        st.passes[p].dirty = true;
      }

      if (st.showGeneratedSource) {
        ImGui::Spacing();
        ImGui::TextUnformatted("Generated GLSL (wrapped)\n");
        ImGui::BeginChild("##glsl", ImVec2(0, 190), true);
        ImGui::TextUnformatted(st.graph.shader(passE).lastFragmentSource().c_str());
        ImGui::EndChild();
      }
    }

    // ---------------- Preview ----------------
    ImGui::TableNextColumn();

    {
      ImGui::TextUnformatted("Preview");
      ImGui::SliderInt("Resolution", &st.previewResolution, 128, 1024);

      const char* passItems = "Buffer A\0Buffer B\0Buffer C\0Buffer D\0Image\0";
      ImGui::Combo("Show", &st.previewPassIndex, passItems);

      // Show the selected texture.
      ImTextureID texId = nullptr;
      if (st.previewPassIndex == 4) {
        texId = (ImTextureID)(intptr_t)st.previewTarget.color().handle();
      } else {
        const auto tex = st.graph.bufferTexture(passEnum(st.previewPassIndex));
        if (tex) texId = (ImTextureID)(intptr_t)tex->handle();
      }

      const float avail = ImGui::GetContentRegionAvail().x;
      const float size = std::min(avail, 480.0f);
      if (texId) {
        ImGui::Image(texId, ImVec2(size, size), ImVec2(0, 1), ImVec2(1, 0));
        updatePreviewMouse(st);
      } else {
        ImGui::TextWrapped("(No texture available for this pass yet)");
      }

      ImGui::Spacing();
      ImGui::Text("Mouse: (%.3f, %.3f)%s", st.mouseU, st.mouseV, st.mouseDown ? " [down]" : "");
      ImGui::Text("Seed: 0x%llX", (unsigned long long)st.seed);
    }

    ImGui::EndTable();
  }

  ImGui::End();
}

} // namespace stellar::game
