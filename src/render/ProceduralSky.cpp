#include "stellar/render/ProceduralSky.h"

#include "stellar/render/Gl.h"

#include <algorithm>

namespace stellar::render {

ProceduralSky::~ProceduralSky() {
  if (vao_) {
    gl::DeleteVertexArrays(1, &vao_);
    vao_ = 0;
  }
}

static const char* kVs = R"GLSL(
#version 330 core
out vec2 vUv;

// Full-screen triangle using gl_VertexID (no VBO needed).
// Covers the viewport with vertices: (-1,-1), (3,-1), (-1,3)
void main() {
  vec2 p;
  if (gl_VertexID == 0) p = vec2(-1.0, -1.0);
  else if (gl_VertexID == 1) p = vec2( 3.0, -1.0);
  else p = vec2(-1.0,  3.0);
  vUv = p * 0.5 + 0.5;
  gl_Position = vec4(p, 0.0, 1.0);
}
)GLSL";

static const char* kFs = R"GLSL(
#version 330 core
in vec2 vUv;
out vec4 oColor;

uniform vec2  uResolution;
uniform float uAspect;
uniform float uTanHalfFovY;
uniform float uTime;

uniform vec3 uCamRight;
uniform vec3 uCamUp;
uniform vec3 uCamForward;
uniform vec3 uCamPos;

uniform int   uSeed;
uniform float uIntensity;

// Stars
uniform float uStarIntensity;
uniform float uStarDensity;
uniform float uStarProb;
uniform float uStarSize;
uniform float uStarTwinkle;

// Star shaping
uniform float uStarClusterStrength;
uniform float uStarClusterFrequency;
uniform float uStarSpikeStrength;

// Milky Way
uniform int   uMilkyWayEnabled;
uniform float uMilkyWayIntensity;
uniform float uMilkyWayFrequency;
uniform float uMilkyWayBandPower;
uniform float uMilkyWayDust;

// Nebula
uniform int   uNebulaEnabled;
uniform float uNebulaIntensity;
uniform float uNebulaFrequency;
uniform float uNebulaThreshold;
uniform float uNebulaSoftness;
uniform float uNebulaBandPower;
uniform float uNebulaParallax;
uniform int   uNebulaSteps;
uniform float uNebulaWorldScale;
uniform float uNebulaDrift;

// -----------------------------------------------------------------------------
// Hash / noise helpers
// -----------------------------------------------------------------------------

float hash12(vec2 p) {
  // Cheap but decent hash for procedural visuals.
  // (sin-based; deterministic)
  return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453123);
}

vec2 hash22(vec2 p) {
  float n = hash12(p);
  return vec2(n, hash12(p + n + 19.19));
}

float hash13(vec3 p) {
  return fract(sin(dot(p, vec3(127.1, 311.7, 74.7))) * 43758.5453123);
}

float fade(float t) {
  // Quintic fade (6t^5 - 15t^4 + 10t^3).
  // Smooth 1st/2nd derivatives at cell boundaries.
  return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
}

// NOTE: GLSL defines deprecated built-in functions named noise1/noise2/noise3/noise4.
// Some drivers still expose these in core profiles. Avoid colliding with them by
// using a project-specific name.
float stellarNoise3(vec3 x) {
  vec3 i = floor(x);
  vec3 f = fract(x);
  vec3 u = vec3(fade(f.x), fade(f.y), fade(f.z));

  // 8 corners
  float n000 = hash13(i + vec3(0,0,0));
  float n100 = hash13(i + vec3(1,0,0));
  float n010 = hash13(i + vec3(0,1,0));
  float n110 = hash13(i + vec3(1,1,0));
  float n001 = hash13(i + vec3(0,0,1));
  float n101 = hash13(i + vec3(1,0,1));
  float n011 = hash13(i + vec3(0,1,1));
  float n111 = hash13(i + vec3(1,1,1));

  float nx00 = mix(n000, n100, u.x);
  float nx10 = mix(n010, n110, u.x);
  float nx01 = mix(n001, n101, u.x);
  float nx11 = mix(n011, n111, u.x);

  float nxy0 = mix(nx00, nx10, u.y);
  float nxy1 = mix(nx01, nx11, u.y);

  return mix(nxy0, nxy1, u.z);
}

float fbm(vec3 p) {
  float sum = 0.0;
  float amp = 0.55;
  float freq = 1.0;
  for (int i = 0; i < 5; ++i) {
    sum += amp * stellarNoise3(p * freq);
    freq *= 2.0;
    amp *= 0.5;
  }
  return sum;
}

vec3 hsv2rgb(vec3 c) {
  vec4 K = vec4(1.0, 2.0/3.0, 1.0/3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 blackbody(float kelvin) {
  // Approximate blackbody temperature (Kelvin) to linear-ish RGB.
  // This is a compact GPU adaptation of common "color temperature" approximations.
  kelvin = clamp(kelvin, 1000.0, 40000.0);
  float t = kelvin / 100.0;
  float r;
  float g;
  float b;

  if (t <= 66.0) {
    r = 1.0;
    g = clamp((99.4708025861 * log(t) - 161.1195681661) / 255.0, 0.0, 1.0);
    if (t <= 19.0) b = 0.0;
    else b = clamp((138.5177312231 * log(t - 10.0) - 305.0447927307) / 255.0, 0.0, 1.0);
  } else {
    r = clamp((329.698727446 * pow(t - 60.0, -0.1332047592)) / 255.0, 0.0, 1.0);
    g = clamp((288.1221695283 * pow(t - 60.0, -0.0755148492)) / 255.0, 0.0, 1.0);
    b = 1.0;
  }
  return vec3(r, g, b);
}

vec3 unitFromSeed(float seedf, float salt) {
  // Deterministic unit vector from seed (spherical sampling).
  vec2 h = hash22(vec2(seedf + salt, 17.17 + salt));
  float z = h.x * 2.0 - 1.0;
  float a = 6.28318530718 * h.y;
  float r = sqrt(max(0.0, 1.0 - z * z));
  return vec3(r * cos(a), z, r * sin(a));
}

mat3 galBasis(float seedf) {
  // A deterministic "galactic" frame that varies with seed.
  // Columns: right, forward, pole.
  vec3 pole = unitFromSeed(seedf, 91.7);
  vec3 ref = (abs(pole.y) < 0.95) ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
  vec3 right = normalize(cross(ref, pole));
  vec3 fwd = cross(pole, right);
  return mat3(right, fwd, pole);
}

vec3 toGal(vec3 dirW, mat3 B) {
  // Express direction in the galactic basis (z is "galactic north").
  return vec3(dot(dirW, B[0]), dot(dirW, B[1]), dot(dirW, B[2]));
}

// -----------------------------------------------------------------------------

vec3 computeRayWorld(vec2 uv) {
  // Map 0..1 UV to NDC -1..1
  vec2 ndc = uv * 2.0 - 1.0;
  vec3 rayCam = normalize(vec3(ndc.x * uAspect * uTanHalfFovY,
                               ndc.y * uTanHalfFovY,
                               1.0));

  // Transform from camera space to world space via basis.
  return normalize(uCamRight * rayCam.x + uCamUp * rayCam.y + uCamForward * rayCam.z);
}

vec3 stars(vec3 dirW, vec3 dirG) {
  // Equirectangular projection (cheap; OK for stars).
  const float PI = 3.14159265359;
  float u = atan(dirW.z, dirW.x) / (2.0 * PI) + 0.5;
  float v = asin(clamp(dirW.y, -1.0, 1.0)) / PI + 0.5;

  vec2 suv = vec2(u, v);
  vec2 p = suv * max(1.0, uStarDensity);
  vec2 cell = floor(p);
  vec2 f = fract(p);

  // Wrap horizontally to avoid a seam at u=0/1.
  // (v wrap is not needed; poles are naturally small.)
  cell.x = mod(cell.x, max(1.0, uStarDensity));

  float seedf = float(uSeed);

  // Coherent star clustering: a sparse fBm field on the unit sphere,
  // biased toward the galactic plane (dirG.z ~= latitude).
  float band = pow(clamp(1.0 - abs(dirG.z), 0.0, 1.0), max(1.0, uMilkyWayBandPower));
  float cStr = max(0.0, uStarClusterStrength);
  float cFreq = max(1e-4, uStarClusterFrequency);
  vec3 cp = dirG * cFreq + vec3(seedf * 0.001, seedf * 0.002, seedf * 0.0017);
  float c0 = fbm(cp * 1.7);
  float cluster = pow(clamp(c0, 0.0, 1.0), 3.0); // sparse peaks
  float clusterBoost = 1.0 + cStr * 5.0 * cluster * (0.25 + 0.75 * band);

  float baseProb = clamp(uStarProb, 0.0, 1.0);
  float prob = clamp(baseProb * clusterBoost, 0.0, 1.0);

  float r0 = hash12(cell + seedf);
  float present = step(1.0 - prob, r0);
  if (present < 0.5) return vec3(0.0);

  vec2 sp = hash22(cell + seedf * 0.123) * 0.92 + 0.04;
  vec2 df = f - sp;
  float d = length(df);

  float sizeRnd = hash12(cell + seedf * 1.731);
  float rad = mix(0.010, 0.045, pow(sizeRnd, 4.0)) / max(0.35, uStarSize);
  // Slightly tighten stars in very dense clusters so clusters read as "sparkle".
  rad *= 1.0 / (1.0 + 0.35 * cStr * cluster);
  float core = smoothstep(rad, 0.0, d);

  // Rare bright stars.
  float big = step(0.997, hash12(cell + seedf * 9.91));
  float halo = smoothstep(rad * (5.0 + 8.0 * big), 0.0, d) * (0.25 + 0.75 * big);

  // Diffraction spikes for the rare bright ones (cheap anisotropic halo).
  float spike = 0.0;
  float spikeStr = clamp(uStarSpikeStrength, 0.0, 1.0);
  if (big > 0.5 && spikeStr > 0.001) {
    float ang = 6.28318530718 * hash12(cell + seedf * 12.17);
    float ca = cos(ang);
    float sa = sin(ang);
    vec2 r = vec2(ca * df.x - sa * df.y, sa * df.x + ca * df.y);
    float w = max(1e-5, rad * 0.25);
    float sx = 1.0 / (1.0 + abs(r.x) / w);
    float sy = 1.0 / (1.0 + abs(r.y) / w);
    spike = (sx + sy) * smoothstep(rad * 10.0, 0.0, d) * spikeStr;
  }

  float br = mix(0.55, 2.20, pow(hash12(cell + seedf * 2.41), 7.0));
  br *= mix(1.0, 4.0, big);
  br *= (0.9 + 2.0 * cStr * cluster * (0.2 + 0.8 * band));

  // Blackbody-ish color temperature variation.
  float tRand = hash12(cell + seedf * 3.17);
  float tempK = mix(1800.0, 12000.0, pow(tRand, 2.8)); // bias toward cooler stars
  vec3 col = blackbody(tempK);
  float tempN = clamp((tempK - 1800.0) / 10200.0, 0.0, 1.0);
  br *= mix(0.85, 1.25, tempN);

  // Twinkle: subtle intensity wobble.
  float tw = 1.0;
  float twA = clamp(uStarTwinkle, 0.0, 1.0);
  if (twA > 0.0) {
    float spd = mix(0.7, 3.5, hash12(cell + seedf * 5.13));
    float ph = 6.28318 * hash12(cell + seedf * 7.77);
    tw += twA * 0.35 * sin(uTime * spd + ph);
  }

  float s = (core + halo + spike) * br * tw;
  return col * s * uStarIntensity;
}

vec3 milkyWay(vec3 dirG) {
  if (uMilkyWayEnabled == 0) return vec3(0.0);

  const float PI = 3.14159265359;
  float seedf = float(uSeed);

  // In galactic space, z acts like latitude: 0 near the galactic plane, +-1 near the poles.
  float lat = clamp(dirG.z, -1.0, 1.0);
  float band = pow(clamp(1.0 - abs(lat), 0.0, 1.0), max(1.0, uMilkyWayBandPower));

  // Seamless longitude coordinates via sin/cos embedding.
  float lon = atan(dirG.y, dirG.x) / (2.0 * PI) + 0.5;
  float ang = 2.0 * PI * lon;
  vec3 q = vec3(cos(ang), sin(ang), lat);

  float freq = max(1e-4, uMilkyWayFrequency);
  vec3 p = q * freq + vec3(seedf * 0.001, seedf * 0.002, seedf * 0.0017);

  // Multi-scale clumps + filaments.
  float n0 = fbm(p * 1.15);
  float n1 = fbm(p * 2.35 + 7.1);
  float cl = mix(n0, n1, 0.45);
  cl = smoothstep(0.35, 0.95, cl);

  // Dust lanes: ridged noise that subtracts from glow.
  float dn = fbm(p * 4.6 + 19.7);
  float ridge = 1.0 - abs(2.0 * dn - 1.0);
  ridge = pow(ridge, 3.0);
  float dust = smoothstep(0.25, 0.85, ridge);
  cl *= (1.0 - clamp(uMilkyWayDust, 0.0, 1.0) * dust);

  // Seeded palette (gentle, desaturated galactic hues).
  float h0 = hash12(vec2(seedf, 0.41));
  float h1 = fract(h0 + 0.08 + 0.45 * hash12(vec2(seedf, 0.99)));

  vec3 cA = hsv2rgb(vec3(h0, 0.35, 0.85));
  vec3 cB = hsv2rgb(vec3(h1, 0.45, 0.95));

  float hueMix = smoothstep(0.0, 1.0, n1);
  vec3 col = mix(cA, cB, hueMix);

  float glow = band * (0.20 + 1.15 * cl);
  return col * glow * uMilkyWayIntensity;
}

vec3 nebula(vec3 dirW, vec3 dirG) {
  if (uNebulaEnabled == 0) return vec3(0.0);

  float seedf = float(uSeed);

  // Palette based on seed.
  float h0 = hash12(vec2(seedf, 0.17));
  float h1 = fract(h0 + 0.12 + 0.45 * hash12(vec2(seedf, 0.73)));

  vec3 cA = hsv2rgb(vec3(h0, 0.55, 0.95));
  vec3 cB = hsv2rgb(vec3(h1, 0.60, 0.80));

  // Emphasize a "galactic plane" band (dirG.z ~= latitude).
  float band = pow(clamp(1.0 - abs(dirG.z), 0.0, 1.0), max(1.0, uNebulaBandPower));

  // Camera movement parallax (0..1). Convert camera position into a tiny noise-space offset.
  float par = clamp(uNebulaParallax, 0.0, 1.0);
  vec3 base = (uCamPos * (1.0 - par)) * uNebulaWorldScale;

  // Time drift as a small translation.
  float drift = uNebulaDrift;
  base += vec3(0.07, 0.03, 0.05) * (uTime * drift);

  int steps = clamp(uNebulaSteps, 2, 24);
  float invSteps = 1.0 / float(steps);

  // Per-pixel jitter to break raymarch banding at low step counts.
  vec2 px = floor(vUv * uResolution);
  float jitter = hash12(px + seedf * 0.17) - 0.5;

  vec3 acc = vec3(0.0);
  float trans = 1.0;

  // A short ray-march through noise space. This is not a physically-correct volume,
  // but produces convincing depth in the background.
  for (int i = 0; i < 24; ++i) {
    if (i >= steps) break;

    float fi = (float(i) + 0.5 + jitter) * invSteps;
    fi = clamp(fi, 0.0, 1.0);

    // Spread samples along the ray.
    float t = mix(0.35, 2.30, fi);

    vec3 p = base + dirW * t;
    p += vec3(seedf * 0.001, seedf * 0.002, seedf * 0.0017);
    p *= max(0.001, uNebulaFrequency);

    float n = fbm(p);
    // Add a little extra structure.
    float n2 = fbm(p * 2.03 + 13.7);
    float d = mix(n, n2, 0.35);

    float th = uNebulaThreshold;
    float soft = max(1e-4, uNebulaSoftness);
    float dens = smoothstep(th, th + soft, d);
    dens *= band;
    dens *= mix(1.10, 0.85, fi); // slight front-bias so the volume reads with depth

    // Color variation along the ray.
    float hueMix = smoothstep(0.15, 0.95, d);
    vec3 cc = mix(cA, cB, hueMix);

    // Cheap front-to-back compositing (adds depth and prevents over-bright fog).
    float a = clamp(dens * (1.15 * invSteps), 0.0, 1.0);
    acc += cc * a * trans;
    trans *= (1.0 - a);

    if (trans < 0.02) break;
  }

  return acc * uNebulaIntensity;
}

void main() {
  vec2 uv = clamp(vUv, 0.0, 1.0);
  vec3 dirW = computeRayWorld(uv);

  float seedf = float(uSeed);
  mat3 G = galBasis(seedf);
  vec3 dirG = toGal(dirW, G);

  vec3 col = vec3(0.0035, 0.0042, 0.0060); // deep space base
  col += milkyWay(dirG);
  col += nebula(dirW, dirG);
  col += stars(dirW, dirG);

  col *= max(0.0, uIntensity);

  oColor = vec4(col, 1.0);
}
)GLSL";

bool ProceduralSky::init(std::string* outError) {
  if (!shader_.build(kVs, kFs, outError)) return false;

  gl::GenVertexArrays(1, &vao_);
  return true;
}

void ProceduralSky::draw(int w,
                         int h,
                         const ProceduralSkySettings& s,
                         const float camRight[3],
                         const float camUp[3],
                         const float camForward[3],
                         const float camPosU[3],
                         float tanHalfFovY,
                         float aspect,
                         float timeSeconds) const {
  if (!vao_) return;
  if (w <= 0 || h <= 0) return;

  shader_.bind();

  shader_.setUniform2f("uResolution", (float)w, (float)h);
  shader_.setUniform1f("uAspect", aspect);
  shader_.setUniform1f("uTanHalfFovY", tanHalfFovY);
  shader_.setUniform1f("uTime", timeSeconds);

  shader_.setUniform3f("uCamRight", camRight[0], camRight[1], camRight[2]);
  shader_.setUniform3f("uCamUp", camUp[0], camUp[1], camUp[2]);
  shader_.setUniform3f("uCamForward", camForward[0], camForward[1], camForward[2]);
  shader_.setUniform3f("uCamPos", camPosU[0], camPosU[1], camPosU[2]);

  shader_.setUniform1i("uSeed", s.seed);
  shader_.setUniform1f("uIntensity", s.intensity);

  shader_.setUniform1f("uStarIntensity", s.starIntensity);
  shader_.setUniform1f("uStarDensity", s.starDensity);
  shader_.setUniform1f("uStarProb", s.starProbability);
  shader_.setUniform1f("uStarSize", s.starSize);
  shader_.setUniform1f("uStarTwinkle", s.starTwinkle);

  shader_.setUniform1f("uStarClusterStrength", s.starClusterStrength);
  shader_.setUniform1f("uStarClusterFrequency", s.starClusterFrequency);
  shader_.setUniform1f("uStarSpikeStrength", s.starSpikeStrength);

  shader_.setUniform1i("uMilkyWayEnabled", s.milkyWayEnabled ? 1 : 0);
  shader_.setUniform1f("uMilkyWayIntensity", s.milkyWayIntensity);
  shader_.setUniform1f("uMilkyWayFrequency", s.milkyWayFrequency);
  shader_.setUniform1f("uMilkyWayBandPower", s.milkyWayBandPower);
  shader_.setUniform1f("uMilkyWayDust", s.milkyWayDust);

  shader_.setUniform1i("uNebulaEnabled", s.nebulaEnabled ? 1 : 0);
  shader_.setUniform1f("uNebulaIntensity", s.nebulaIntensity);
  shader_.setUniform1f("uNebulaFrequency", s.nebulaFrequency);
  shader_.setUniform1f("uNebulaThreshold", s.nebulaThreshold);
  shader_.setUniform1f("uNebulaSoftness", s.nebulaSoftness);
  shader_.setUniform1f("uNebulaBandPower", s.nebulaBandPower);
  shader_.setUniform1f("uNebulaParallax", s.nebulaParallax);
  shader_.setUniform1i("uNebulaSteps", s.nebulaSteps);
  shader_.setUniform1f("uNebulaWorldScale", s.nebulaWorldScale);
  shader_.setUniform1f("uNebulaDrift", s.nebulaDrift);

  gl::BindVertexArray(vao_);
  glDrawArrays(GL_TRIANGLES, 0, 3);
  gl::BindVertexArray(0);
}

} // namespace stellar::render
