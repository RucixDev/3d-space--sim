#include "stellar/render/GpuSurfaceCache.h"

#include "stellar/core/Hash.h"

#include "stellar/render/Gl.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace stellar::render {

namespace {

static const char* kVS = R"GLSL(
#version 330 core

out vec2 vUv;

void main() {
  vec2 verts[3] = vec2[3](vec2(-1.0, -1.0), vec2(3.0, -1.0), vec2(-1.0, 3.0));
  vec2 p = verts[gl_VertexID];
  gl_Position = vec4(p, 0.0, 1.0);
  vUv = 0.5 * p + vec2(0.5);
}
)GLSL";

// A single fragment shader is used for both albedo and tangent-space normal maps.
//
// uMode:
//   0 -> output RGBA albedo
//   1 -> output RGBA tangent-space normal (XYZ encoded into RGB, A=1)
static const char* kFS = R"GLSL(
#version 330 core

in vec2 vUv;

uniform int uKind;    // render::SurfaceKind (0..)
uniform int uMode;    // 0=albedo, 1=normal
uniform int uSeedLo;  // bitwise seed (lower 32)
uniform int uSeedHi;  // bitwise seed (upper 32)
uniform vec2 uInvSize; // (1/w, 1/h)

out vec4 FragColor;

const float PI = 3.1415926535897932384626433832795;
const float TWO_PI = 6.2831853071795864769252867665590;

float clamp01(float x) { return clamp(x, 0.0, 1.0); }

// Quintic fade (Perlin-style): 6t^5 - 15t^4 + 10t^3
float fade(float t) {
  return t*t*t*(t*(t*6.0 - 15.0) + 10.0);
}

float smoothstep01(float e0, float e1, float x) {
  float t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3.0 - 2.0 * t);
}

uint pcg_hash(uint v) {
  // "PCG hash" permutation (good quality vs cost).
  uint state = v * 747796405u + 2891336453u;
  uint word  = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
  return (word >> 22u) ^ word;
}

uint seedMix() {
  // Fold two 32-bit halves into one 32-bit seed.
  uint lo = uint(uSeedLo);
  uint hi = uint(uSeedHi);
  return pcg_hash(lo ^ pcg_hash(hi + 0x9E3779B9u));
}

float hash01(uvec3 p, uint seed) {
  uint h = pcg_hash(p.x + pcg_hash(p.y + pcg_hash(p.z + seed)));
  return float(h) / 4294967295.0;
}

float smoothNoise3D(vec3 x, uint seed) {
  vec3 p0 = floor(x);
  vec3 f  = fract(x);

  uvec3 i0 = uvec3(p0);
  uvec3 i1 = i0 + uvec3(1,0,0);
  uvec3 j1 = i0 + uvec3(0,1,0);
  uvec3 k1 = i0 + uvec3(0,0,1);

  float n000 = hash01(i0, seed);
  float n100 = hash01(i1, seed);
  float n010 = hash01(j1, seed);
  float n110 = hash01(i0 + uvec3(1,1,0), seed);
  float n001 = hash01(k1, seed);
  float n101 = hash01(i0 + uvec3(1,0,1), seed);
  float n011 = hash01(i0 + uvec3(0,1,1), seed);
  float n111 = hash01(i0 + uvec3(1,1,1), seed);

  float u = fade(f.x);
  float v = fade(f.y);
  float w = fade(f.z);

  float x00 = mix(n000, n100, u);
  float x10 = mix(n010, n110, u);
  float x01 = mix(n001, n101, u);
  float x11 = mix(n011, n111, u);

  float y0v = mix(x00, x10, v);
  float y1v = mix(x01, x11, v);
  return mix(y0v, y1v, w);
}

float fbm3D(vec3 p, uint seed, int octaves, float lacunarity, float gain) {
  float amp = 0.5;
  float freq = 1.0;
  float sum = 0.0;
  for (int i = 0; i < 8; ++i) {
    if (i >= octaves) break;
    sum += amp * smoothNoise3D(p * freq, seed + uint(i) * 1013u);
    freq *= lacunarity;
    amp *= gain;
  }
  return sum;
}

float ampSum(int octaves, float gain) {
  float amp = 0.5;
  float s = 0.0;
  for (int i = 0; i < 8; ++i) {
    if (i >= octaves) break;
    s += amp;
    amp *= gain;
  }
  return max(1.0e-6, s);
}

float fbm01(vec3 p, uint seed, int octaves, float lacunarity, float gain) {
  return clamp01(fbm3D(p, seed, octaves, lacunarity, gain) / ampSum(octaves, gain));
}

float ridged01(float n01) {
  return 1.0 - abs(n01 * 2.0 - 1.0);
}

vec3 noiseVec3(vec3 p, uint seed, float freq) {
  vec3 q = p * freq;
  float x = fbm01(q, seed ^ 0xA11CE5u, 4, 2.0, 0.5) - 0.5;
  float y = fbm01(q, seed ^ 0xBADC0FFEu, 4, 2.0, 0.5) - 0.5;
  float z = fbm01(q, seed ^ 0xC0FFEEu, 4, 2.0, 0.5) - 0.5;
  return vec3(x, y, z);
}

vec3 warpDir(vec3 p, uint seed, float freq, float amp) {
  vec3 w = noiseVec3(p, seed, freq);
  return normalize(p + w * amp);
}

vec3 latLonToSphere(float lat, float lon) {
  float cosLat = cos(lat);
  float sinLat = sin(lat);
  float cosLon = cos(lon);
  float sinLon = sin(lon);
  return vec3(cosLat * cosLon, sinLat, cosLat * sinLon);
}

vec3 randDir(uint seed, uint idx, uint salt) {
  float u = hash01(uvec3(idx, 0u, salt), seed);
  float v = hash01(uvec3(idx, 1u, salt), seed);
  float z = u * 2.0 - 1.0;
  float phi = TWO_PI * v;
  float r = sqrt(max(0.0, 1.0 - z*z));
  return vec3(r * cos(phi), z, r * sin(phi));
}

// -----------------------------------------------------------------------------
// Tectonic plates (spherical Voronoi), craters, vortices

const uint SALT_PLATE  = 0x504C4154u; // 'PLAT'
const uint SALT_CRATER = 0x43524154u; // 'CRAT'
const uint SALT_VGAS   = 0x56474153u; // 'VGAS'
const uint SALT_VCLD   = 0x56434C44u; // 'VCLD'

int plateBaseCount(int kind) {
  if (kind == 0) return 26; // Rocky
  if (kind == 1) return 22; // Desert
  if (kind == 2) return 18; // Ocean
  if (kind == 3) return 20; // Ice
  if (kind == 4) return 12; // Gas
  if (kind == 6) return 12; // Clouds
  return 0;
}

int plateCount(int kind, uint seed) {
  int base = plateBaseCount(kind);
  if (base <= 0) return 0;
  int jitter = int(pcg_hash(seed ^ (SALT_PLATE + uint(kind) * 101u)) % 11u) - 4; // -4..+6
  return clamp(base + jitter, 8, 48);
}

float plateHeightBias(int kind, uint seed, int idx) {
  float t = hash01(uvec3(uint(idx), 2u, SALT_PLATE), seed);
  if (kind == 2) {
    // Bias towards smaller values so water dominates.
    return pow(t, 1.35);
  }
  if (kind == 4 || kind == 6) {
    return mix(0.35, 0.65, t);
  }
  return mix(0.25, 0.85, t);
}

float plateTempBias(int kind, uint seed, int idx) {
  float t = hash01(uvec3(uint(idx), 3u, SALT_PLATE), seed);
  if (kind == 2) return mix(-0.12, 0.12, t);
  if (kind == 4 || kind == 6) return mix(-0.08, 0.08, t);
  return mix(-0.10, 0.10, t);
}

float plateMoistBias(int kind, uint seed, int idx) {
  float t = hash01(uvec3(uint(idx), 4u, SALT_PLATE), seed);
  if (kind == 2) return mix(-0.20, 0.20, t);
  if (kind == 4 || kind == 6) return mix(-0.08, 0.08, t);
  return mix(-0.10, 0.10, t);
}

struct PlateSample {
  int idx;
  float boundary01;
  float bestDot;
  float secondDot;
};

PlateSample samplePlates(vec3 p, uint seed, int kind) {
  PlateSample s;
  s.idx = 0;
  s.boundary01 = 0.0;
  s.bestDot = -2.0;
  s.secondDot = -2.0;

  int count = plateCount(kind, seed);
  if (count <= 0) return s;

  float best = -2.0;
  float second = -2.0;
  int bestIdx = 0;

  uint plateSeed = seed ^ SALT_PLATE;
  for (int i = 0; i < 48; ++i) {
    if (i >= count) break;
    vec3 d = randDir(plateSeed, uint(i), SALT_PLATE);
    float dd = dot(p, d);
    if (dd > best) {
      second = best;
      best = dd;
      bestIdx = i;
    } else if (dd > second) {
      second = dd;
    }
  }

  float diff = best - second;
  float t = clamp01((0.11 - diff) / 0.11);
  s.idx = bestIdx;
  s.bestDot = best;
  s.secondDot = second;
  s.boundary01 = smoothstep01(0.0, 1.0, t);
  return s;
}

int craterCount(int kind, uint seed) {
  // Oceans/gas/stars/clouds do not use this crater field.
  if (kind == 2 || kind == 4 || kind == 5 || kind == 6) return 0;

  int base = 48;
  if (kind == 1) base = 38;
  if (kind == 3) base = 52;

  int jitter = int(pcg_hash(seed ^ SALT_CRATER) % 27u) - 10; // -10..+16
  return clamp(base + jitter, 0, 96);
}

float applyCraters(float h, vec3 p, uint seed, int kind, float strengthMul,
                  out float craterInterior01, out float craterRim01) {
  craterInterior01 = 0.0;
  craterRim01 = 0.0;

  int count = craterCount(kind, seed);
  if (count <= 0 || strengthMul <= 0.0) return h;

  float rMinDeg = 4.0;
  float rMaxDeg = 18.0;
  float depthBase = 0.085;
  float rimBase = 0.25;

  if (kind == 1) { // Desert
    rMinDeg = 5.0;
    rMaxDeg = 20.0;
    depthBase = 0.070;
    rimBase = 0.22;
  } else if (kind == 3) { // Ice
    rMinDeg = 6.0;
    rMaxDeg = 22.0;
    depthBase = 0.080;
    rimBase = 0.28;
  }

  uint craterSeed = seed ^ SALT_CRATER;

  for (int i = 0; i < 96; ++i) {
    if (i >= count) break;
    uint idx = uint(i);

    vec3 dir = randDir(craterSeed, idx, SALT_CRATER);

    float tR = hash01(uvec3(idx, 2u, SALT_CRATER), seed);
    float rad = radians(mix(rMinDeg, rMaxDeg, tR));
    float cosR = cos(rad);
    float invSpan = 1.0 / max(1.0e-6, (1.0 - cosR));

    float d = dot(p, dir);
    if (d <= cosR) continue;

    float t = clamp01((d - cosR) * invSpan); // 0=edge -> 1=center
    float w = smoothstep01(0.0, 1.0, t);

    float depth = depthBase * mix(0.65, 1.15, hash01(uvec3(idx, 3u, SALT_CRATER), seed));
    float rim = clamp(rimBase * mix(0.85, 1.15, hash01(uvec3(idx, 4u, SALT_CRATER), seed)), 0.0, 1.0);

    h -= depth * strengthMul * w;
    craterInterior01 = max(craterInterior01, w);

    if (rim > 0.0) {
      float ring = smoothstep01(0.06, 0.22, t) * (1.0 - smoothstep01(0.22, 0.58, t));
      h += depth * strengthMul * rim * ring;
      craterRim01 = max(craterRim01, ring);
    }
  }

  return clamp01(h);
}

int vortexCount(uint seed, uint salt) {
  // 3..8
  return int(3u + (pcg_hash(seed ^ salt) % 6u));
}

vec3 vortexTint(int i, uint seed) {
  int v = (i + int(seed & 7u)) % 5;
  if (v == 0) return vec3(0.95, 0.60, 0.45);
  if (v == 1) return vec3(0.85, 0.92, 1.00);
  if (v == 2) return vec3(1.00, 0.92, 0.70);
  if (v == 3) return vec3(0.72, 0.78, 0.86);
  return vec3(0.92, 0.82, 0.62);
}

// -----------------------------------------------------------------------------
// Biomes

vec3 biomeColor(float temp01, float moist01) {
  temp01 = clamp01(temp01);
  moist01 = clamp01(moist01);

  vec3 desert = vec3(0.86, 0.78, 0.45);
  vec3 savanna = vec3(0.62, 0.70, 0.32);
  vec3 grass = vec3(0.34, 0.60, 0.26);
  vec3 forest = vec3(0.16, 0.42, 0.18);
  vec3 jungle = vec3(0.12, 0.50, 0.22);
  vec3 steppe = vec3(0.55, 0.50, 0.26);
  vec3 tundra = vec3(0.62, 0.62, 0.54);

  if (temp01 < 0.22) {
    float s = smoothstep01(0.00, 0.20, (0.22 - temp01) / 0.22);
    return mix(tundra, vec3(0.90, 0.95, 1.00), s);
  }

  if (temp01 < 0.42) {
    if (moist01 > 0.55) return mix(steppe, forest, (moist01 - 0.55) / 0.45);
    return mix(tundra, steppe, moist01 / 0.55);
  }

  if (temp01 < 0.72) {
    if (moist01 > 0.70) return mix(forest, jungle, (moist01 - 0.70) / 0.30);
    if (moist01 > 0.45) return mix(grass, forest, (moist01 - 0.45) / 0.25);
    if (moist01 > 0.25) return mix(savanna, grass, (moist01 - 0.25) / 0.20);
    return mix(desert, savanna, moist01 / 0.25);
  }

  if (moist01 > 0.70) return jungle;
  if (moist01 > 0.48) return mix(savanna, jungle, (moist01 - 0.48) / 0.22);
  if (moist01 > 0.30) return mix(desert, savanna, (moist01 - 0.30) / 0.18);
  return desert;
}

// -----------------------------------------------------------------------------
// Surface functions

float surfaceCloudAlpha(uint seed, vec3 p, float lat, float lon) {
  vec3 q = warpDir(p, seed ^ 0xC10Du, 2.2, 0.22);

  float n1 = fbm01(q * 3.6, seed ^ 0xC10Du, 6, 2.0, 0.5);
  float n2 = fbm01(q * 10.5 + vec3(0.0, lat * 0.35, lon * 0.15), seed ^ 0xBEEFu, 4, 2.1, 0.55);
  float n = 0.65 * n1 + 0.35 * n2;
  float d = smoothstep01(0.55, 0.82, n);

  float bandWarp = fbm01(q * 2.8, seed ^ 0xC1A0u, 4, 2.0, 0.55);
  float band = 0.5 + 0.5 * sin(lat * (10.0 + 8.0 * bandWarp) + (n2 - 0.5) * 3.0);
  d *= 0.78 + 0.22 * band;

  // Storms.
  int count = vortexCount(seed, SALT_VCLD);
  uint vSeed = seed ^ SALT_VCLD;
  for (int i = 0; i < 10; ++i) {
    if (i >= count) break;

    vec3 dir0 = randDir(vSeed, uint(i), SALT_VCLD);
    dir0.y = clamp(dir0.y, -0.92, 0.92);
    vec3 dir = normalize(dir0);

    float radDeg = mix(8.0, 28.0, hash01(uvec3(uint(i), 2u, SALT_VCLD), seed));
    float cosR = cos(radians(radDeg));
    float invSpan = 1.0 / max(1.0e-6, (1.0 - cosR));

    float dd = dot(p, dir);
    if (dd <= cosR) continue;

    float u = clamp01((dd - cosR) * invSpan);
    float w = smoothstep01(0.0, 1.0, u);
    float str = mix(0.35, 0.85, hash01(uvec3(uint(i), 3u, SALT_VCLD), seed));

    float swirl = 0.5 + 0.5 * sin(lon * 4.5 + lat * 7.0 + float(seed & 7u));
    d += (swirl * 0.55 + 0.45) * w * str * 0.30;
  }

  return clamp01(d);
}

float surfaceHeight(int kind, uint seed, vec3 p, float lat, float lon) {
  // Height fields are intentionally approximate; only relative variation matters.

  // Shared plate warp for rocky/desert/ocean/ice.
  vec3 pWarp = p;
  if (kind == 0) pWarp = warpDir(p, seed ^ 0xA11CE5u, 2.1, 0.22);
  if (kind == 1) pWarp = warpDir(p, seed ^ 0xD35E7u, 1.8, 0.18);
  if (kind == 2) pWarp = warpDir(p, seed ^ 0x0CE4u, 1.6, 0.20);
  if (kind == 3) pWarp = warpDir(p, seed ^ 0x1CEu,  1.9, 0.16);

  if (kind == 0) { // Rocky
    PlateSample ps = samplePlates(pWarp, seed, kind);
    float plateH = plateHeightBias(kind, seed, ps.idx);

    float macro = fbm01(pWarp * 2.8, seed ^ 0xA11CE5u, 6, 2.0, 0.5);
    float rN = fbm01(pWarp * 7.5, seed ^ 0xBADC0FFEu, 4, 2.1, 0.55);
    float ridged = ridged01(rN);
    float dust = fbm01(pWarp * 15.0, seed ^ 0xC0FFEEu, 3, 2.4, 0.55);

    float mountain = pow(ps.boundary01, 1.65) * (0.18 + 0.22 * dust);

    float h = clamp01(0.52 * plateH + 0.30 * macro + 0.18 * ridged + mountain);
    h = clamp01(h * (0.92 + 0.10 * dust));

    float cIn, cRim;
    h = applyCraters(h, p, seed, kind, 1.0, cIn, cRim);
    return h;
  }

  if (kind == 1) { // Desert
    PlateSample ps = samplePlates(pWarp, seed, kind);
    float plateH = plateHeightBias(kind, seed, ps.idx);

    float macro = fbm01(pWarp * 2.2, seed ^ 0xD35E7u, 6, 2.0, 0.5);
    float warp = fbm01(pWarp * 3.4, seed ^ 0x51DEu, 4, 2.1, 0.55);

    float band = 0.5 + 0.5 * sin((lat * 16.0) + (lon * 2.2) + (warp - 0.5) * 3.2);

    float rock = fbm01(pWarp * 6.0, seed ^ 0x0BADC0DEu, 4, 2.1, 0.52);
    float rockMask = smoothstep01(0.70, 0.88, rock);

    float mountain = pow(ps.boundary01, 1.55) * (0.12 + 0.16 * rock);

    float h = clamp01(0.55 * plateH + 0.22 * macro + 0.20 * band + 0.08 * rockMask + mountain);
    float cIn, cRim;
    h = applyCraters(h, p, seed, kind, 0.75, cIn, cRim);
    return h;
  }

  if (kind == 2) { // Ocean
    PlateSample ps = samplePlates(pWarp, seed, kind);
    float plateH = plateHeightBias(kind, seed, ps.idx);

    float macro = fbm01(pWarp * 2.1, seed ^ 0x0CE4u, 6, 2.0, 0.5);
    float detail = fbm01(pWarp * 6.5, seed ^ 0xE1E7u, 5, 2.1, 0.55);

    float mountain = pow(ps.boundary01, 1.85) * (0.18 + 0.26 * detail);
    float elevRaw = clamp01(0.68 * plateH + 0.32 * macro + mountain);

    float tSea = float((uint(uSeedLo) >> 11u) & 255u) / 255.0;
    float sea = 0.50 + tSea * 0.09;

    float land = smoothstep01(sea - 0.02, sea + 0.02, elevRaw);
    float landElev = clamp01((elevRaw - sea) / max(1.0e-6, (1.0 - sea)));

    float wave = fbm01(pWarp * 18.0, seed ^ 0xA57EA11Cu, 3, 2.2, 0.55);
    float hOcean = 0.03 * (wave - 0.5) * (0.85 + 0.15 * abs(sin(lat)));
    float hLand = 0.18 + 0.82 * landElev;

    return clamp01(mix(hOcean, hLand, land));
  }

  if (kind == 3) { // Ice
    PlateSample ps = samplePlates(pWarp, seed, kind);
    float plateH = plateHeightBias(kind, seed, ps.idx);

    float n = fbm01(pWarp * 3.3, seed ^ 0x1CEu, 6, 2.0, 0.5);
    float crack = fbm01(pWarp * 20.0, seed ^ 0xC24Cu, 3, 2.4, 0.55);
    float crackMask = smoothstep01(0.78, 0.93, crack);

    float ridge = pow(ps.boundary01, 1.35) * 0.28;

    float h = clamp01(0.52 * plateH + 0.32 * n + 0.16 * crackMask + ridge);

    float eq = 1.0 - abs(sin(lat));
    float ripple = fbm01(pWarp * 9.0, seed ^ 0xF00Du, 4, 2.1, 0.55);
    h = clamp01(h + (ripple - 0.5) * 0.05 * eq);

    float cIn, cRim;
    h = applyCraters(h, p, seed, kind, 0.85, cIn, cRim);
    return h;
  }

  if (kind == 4) { // GasGiant
    vec3 q = warpDir(p, seed ^ 0x6A59u, 2.6, 0.25);
    float warp = fbm01(q * 2.6, seed ^ 0x6A59u, 5, 2.0, 0.55);

    float bandFreq = 20.0 + 10.0 * warp;
    float band = 0.5 + 0.5 * sin((lat * bandFreq) + (warp - 0.5) * 4.2);
    float turb = fbm01(q * 7.0, seed ^ 0x7B3Bu, 4, 2.1, 0.55);

    float h = 0.5 + 0.10 * (band - 0.5) + 0.08 * (turb - 0.5);

    int count = vortexCount(seed, SALT_VGAS);
    uint vSeed = seed ^ SALT_VGAS;
    for (int i = 0; i < 10; ++i) {
      if (i >= count) break;

      vec3 dir0 = randDir(vSeed, uint(i), SALT_VGAS);
      dir0.y = clamp(dir0.y, -0.92, 0.92);
      vec3 dir = normalize(dir0);

      float radDeg = mix(8.0, 28.0, hash01(uvec3(uint(i), 2u, SALT_VGAS), seed));
      float cosR = cos(radians(radDeg));
      float invSpan = 1.0 / max(1.0e-6, (1.0 - cosR));

      float dd = dot(p, dir);
      if (dd <= cosR) continue;

      float u = clamp01((dd - cosR) * invSpan);
      float w = smoothstep01(0.0, 1.0, u);
      float str = mix(0.35, 0.85, hash01(uvec3(uint(i), 3u, SALT_VGAS), seed));
      float swirl = 0.5 + 0.5 * sin(lon * 3.2 + lat * 5.1 + float(seed & 7u));
      h += (swirl - 0.5) * 0.08 * w * str;
    }

    return clamp01(h);
  }

  if (kind == 5) { // Star
    float coarse = fbm01(p * 5.0, seed ^ 0x57A4u, 5, 2.0, 0.55);
    float fine   = fbm01(p * 18.0, seed ^ 0xF1AEu, 3, 2.4, 0.55);
    float h = 0.55 + 0.55 * coarse + 0.20 * (fine - 0.5);
    return clamp01(h);
  }

  if (kind == 6) { // Clouds
    return surfaceCloudAlpha(seed, p, lat, lon);
  }

  return 0.5;
}

vec3 surfaceAlbedo(int kind, uint seed, vec3 p, float lat, float lon, out float outAlpha) {
  outAlpha = 1.0;

  int paletteIx = int((uint(uSeedLo) >> 21u) & 7u);

  vec3 pWarp = p;
  if (kind == 0) pWarp = warpDir(p, seed ^ 0xA11CE5u, 2.1, 0.22);
  if (kind == 1) pWarp = warpDir(p, seed ^ 0xD35E7u, 1.8, 0.18);
  if (kind == 2) pWarp = warpDir(p, seed ^ 0x0CE4u, 1.6, 0.20);
  if (kind == 3) pWarp = warpDir(p, seed ^ 0x1CEu,  1.9, 0.16);

  if (kind == 0) { // Rocky
    PlateSample ps = samplePlates(pWarp, seed, kind);
    float plateH = plateHeightBias(kind, seed, ps.idx);

    float macro = fbm01(pWarp * 2.8, seed ^ 0xA11CE5u, 6, 2.0, 0.5);
    float rN = fbm01(pWarp * 7.5, seed ^ 0xBADC0FFEu, 4, 2.1, 0.55);
    float ridged = ridged01(rN);
    float dust = fbm01(pWarp * 15.0, seed ^ 0xC0FFEEu, 3, 2.4, 0.55);

    float mountain = pow(ps.boundary01, 1.65) * (0.18 + 0.22 * dust);
    float h = clamp01(0.52 * plateH + 0.30 * macro + 0.18 * ridged + mountain);
    h = clamp01(h * (0.92 + 0.10 * dust));

    float craterIn, craterRim;
    h = applyCraters(h, p, seed, kind, 1.0, craterIn, craterRim);

    vec3 baseLo = vec3(0.18, 0.16, 0.15);
    vec3 baseHi = vec3(0.62, 0.56, 0.48);
    if ((paletteIx % 3) == 1) {
      baseLo = vec3(0.20, 0.12, 0.10);
      baseHi = vec3(0.70, 0.46, 0.32);
    } else if ((paletteIx % 3) == 2) {
      baseLo = vec3(0.10, 0.10, 0.12);
      baseHi = vec3(0.46, 0.46, 0.50);
    }

    vec3 col = mix(baseLo, baseHi, h);
    col = mix(col, col * 0.65, (1.0 - ridged) * 0.35);
    col *= (0.85 + 0.28 * dust);

    col = mix(col, vec3(0.74, 0.72, 0.70), ps.boundary01 * 0.18);

    col = mix(col, col * 0.55, craterIn * 0.55);
    col = mix(col, vec3(0.95, 0.93, 0.90), craterRim * 0.28);

    float speck = fbm01(pWarp * 26.0, seed ^ 0xFEEDBEEFu, 2, 2.5, 0.6);
    float speckMask = smoothstep01(0.84, 0.97, speck);
    col = mix(col, vec3(0.88, 0.84, 0.78), speckMask * 0.32);

    return col;
  }

  if (kind == 1) { // Desert
    PlateSample ps = samplePlates(pWarp, seed, kind);
    float plateH = plateHeightBias(kind, seed, ps.idx);

    float macro = fbm01(pWarp * 2.2, seed ^ 0xD35E7u, 6, 2.0, 0.5);
    float warp = fbm01(pWarp * 3.4, seed ^ 0x51DEu, 4, 2.1, 0.55);
    float band = 0.5 + 0.5 * sin((lat * 16.0) + (lon * 2.2) + (warp - 0.5) * 3.2);
    float rock = fbm01(pWarp * 6.0, seed ^ 0x0BADC0DEu, 4, 2.1, 0.52);
    float rockMask = smoothstep01(0.70, 0.88, rock);
    float mountain = pow(ps.boundary01, 1.55) * (0.12 + 0.16 * rock);

    float h = clamp01(0.55 * plateH + 0.22 * macro + 0.20 * band + 0.08 * rockMask + mountain);

    float craterIn, craterRim;
    h = applyCraters(h, p, seed, kind, 0.75, craterIn, craterRim);

    vec3 sandLo = vec3(0.62, 0.52, 0.28);
    vec3 sandHi = vec3(0.92, 0.86, 0.50);
    if ((paletteIx & 1) != 0) {
      sandLo = vec3(0.56, 0.34, 0.22);
      sandHi = vec3(0.88, 0.62, 0.40);
    }

    vec3 col = mix(sandLo, sandHi, clamp01(0.35 + 0.65 * macro));
    col *= (0.90 + 0.18 * band);

    col = mix(col, vec3(0.34, 0.28, 0.22), rockMask * 0.65);
    col = mix(col, vec3(0.70, 0.64, 0.52), ps.boundary01 * 0.16);

    col = mix(col, col * 0.62, craterIn * 0.50);
    col = mix(col, vec3(0.98, 0.94, 0.86), craterRim * 0.18);

    col = mix(col, col * 0.92, (1.0 - h) * 0.25);

    return col;
  }

  if (kind == 2) { // Ocean
    PlateSample ps = samplePlates(pWarp, seed, kind);

    float plateH = plateHeightBias(kind, seed, ps.idx);
    float macro = fbm01(pWarp * 2.1, seed ^ 0x0CE4u, 6, 2.0, 0.5);
    float detail = fbm01(pWarp * 6.5, seed ^ 0xE1E7u, 5, 2.1, 0.55);
    float mountain = pow(ps.boundary01, 1.85) * (0.18 + 0.26 * detail);
    float elevRaw = clamp01(0.68 * plateH + 0.32 * macro + mountain);

    float tSea = float((uint(uSeedLo) >> 11u) & 255u) / 255.0;
    float sea = 0.50 + tSea * 0.09;

    float land = smoothstep01(sea - 0.02, sea + 0.02, elevRaw);
    float landElev = clamp01((elevRaw - sea) / max(1.0e-6, (1.0 - sea)));
    float coast = clamp01(1.0 - abs(elevRaw - sea) * 22.0);

    vec3 waterDeep = vec3(0.02, 0.07, 0.20);
    vec3 waterShallow = vec3(0.07, 0.26, 0.52);
    if ((paletteIx & 2) != 0) {
      waterDeep = vec3(0.03, 0.09, 0.18);
      waterShallow = vec3(0.08, 0.30, 0.40);
    }
    vec3 water = mix(waterDeep, waterShallow, coast);

    float lat01 = abs(sin(lat));

    float temp = 1.0 - pow(lat01, 0.95);
    temp += plateTempBias(kind, seed, ps.idx);
    temp -= landElev * 0.35;
    temp = clamp01(temp);

    float moist = fbm01(pWarp * 3.5, seed ^ 0xC0A57u, 5, 2.0, 0.5);
    moist = clamp01((moist - 0.35) * 1.45);
    moist *= (0.62 + 0.38 * cos(lat * 2.0));
    moist = clamp01(moist + plateMoistBias(kind, seed, ps.idx));

    vec3 landCol = biomeColor(temp, moist);

    // Beaches.
    float beach = coast * (1.0 - smoothstep01(sea + 0.01, sea + 0.04, elevRaw));
    landCol = mix(landCol, vec3(0.92, 0.86, 0.62), beach);

    float mountainMask = max(smoothstep01(0.72, 0.92, landElev), ps.boundary01 * 0.55);
    landCol = mix(landCol, vec3(0.52, 0.50, 0.48), mountainMask * 0.55);

    float snowLat = smoothstep01(0.62, 0.88, lat01);
    float snowElev = smoothstep01(0.80, 0.98, landElev);
    float snow = clamp01(max(snowLat, snowElev) * (0.55 + 0.45 * (1.0 - temp)));
    landCol = mix(landCol, vec3(0.90, 0.95, 1.00), snow);

    vec3 col = mix(water, landCol, land);

    float cloud = fbm01(pWarp * 4.0, seed ^ 0xC10Du, 5, 2.0, 0.5);
    float cloudMask = smoothstep01(0.74, 0.93, cloud) * (1.0 - land * 0.35);
    col = mix(col, vec3(1.0), cloudMask * 0.08);

    return col;
  }

  if (kind == 3) { // Ice
    PlateSample ps = samplePlates(pWarp, seed, kind);
    float plateH = plateHeightBias(kind, seed, ps.idx);

    float n = fbm01(pWarp * 3.3, seed ^ 0x1CEu, 6, 2.0, 0.5);
    float crack = fbm01(pWarp * 20.0, seed ^ 0xC24Cu, 3, 2.4, 0.55);
    float crackMask = smoothstep01(0.78, 0.93, crack);
    float ridge = pow(ps.boundary01, 1.35) * 0.28;

    float h = clamp01(0.52 * plateH + 0.32 * n + 0.16 * crackMask + ridge);

    float craterIn, craterRim;
    h = applyCraters(h, p, seed, kind, 0.85, craterIn, craterRim);

    vec3 iceLo = vec3(0.70, 0.82, 0.95);
    vec3 iceHi = vec3(0.93, 0.98, 1.00);
    vec3 col = mix(iceLo, iceHi, clamp01(0.15 + 0.85 * n));

    col = mix(col, vec3(0.28, 0.32, 0.38), crackMask * 0.35);

    float dust = fbm01(pWarp * 10.0, seed ^ 0x1CEDu, 4, 2.1, 0.55);
    float eq = 1.0 - abs(sin(lat));
    float dustMask = clamp01((dust - 0.45) * 1.45) * eq;
    col = mix(col, vec3(0.62, 0.56, 0.44), dustMask * 0.35);

    col = mix(col, vec3(0.86, 0.90, 0.98), ps.boundary01 * 0.10);

    col = mix(col, col * 0.75, craterIn * 0.55);
    col = mix(col, vec3(1.0), craterRim * 0.18);

    float pole = smoothstep01(0.62, 0.92, abs(sin(lat)));
    col = mix(col, vec3(0.92, 0.98, 1.00), pole * 0.15);

    return col;
  }

  if (kind == 4) { // GasGiant
    vec3 q = warpDir(p, seed ^ 0x6A59u, 2.6, 0.25);
    float warp = fbm01(q * 2.6, seed ^ 0x6A59u, 5, 2.0, 0.55);
    float bandFreq = 22.0 + 12.0 * warp;
    float band = 0.5 + 0.5 * sin((lat * bandFreq) + (warp - 0.5) * 4.2);
    float fine = fbm01(q * 12.0, seed ^ 0x6A5Au, 3, 2.3, 0.55);
    float turb = fbm01(q * 7.0, seed ^ 0x7B3Bu, 4, 2.1, 0.55);

    float t = clamp01(0.60 * band + 0.25 * turb + 0.15 * fine);

    vec3 lo = vec3(0.84, 0.68, 0.42);
    vec3 hi = vec3(0.96, 0.88, 0.66);
    if ((paletteIx & 1) != 0) {
      lo = vec3(0.70, 0.76, 0.86);
      hi = vec3(0.95, 0.94, 0.88);
    }

    vec3 col = mix(lo, hi, t);

    float pole = smoothstep01(0.55, 0.92, abs(sin(lat)));
    col = mix(col, vec3(0.68, 0.76, 0.88), pole * 0.18);

    int count = vortexCount(seed, SALT_VGAS);
    uint vSeed = seed ^ SALT_VGAS;
    for (int i = 0; i < 10; ++i) {
      if (i >= count) break;

      vec3 dir0 = randDir(vSeed, uint(i), SALT_VGAS);
      dir0.y = clamp(dir0.y, -0.92, 0.92);
      vec3 dir = normalize(dir0);

      float radDeg = mix(8.0, 28.0, hash01(uvec3(uint(i), 2u, SALT_VGAS), seed));
      float cosR = cos(radians(radDeg));
      float invSpan = 1.0 / max(1.0e-6, (1.0 - cosR));

      float dd = dot(p, dir);
      if (dd <= cosR) continue;

      float u = clamp01((dd - cosR) * invSpan);
      float w = smoothstep01(0.0, 1.0, u);
      float str = mix(0.35, 0.85, hash01(uvec3(uint(i), 3u, SALT_VGAS), seed));

      float swirl = 0.5 + 0.5 * sin(lon * 2.8 + lat * 7.5 + float(seed & 7u));
      float local = clamp01(t + (swirl - 0.5) * 0.25 * w * str);
      col = mix(col, mix(lo, hi, local), 0.35 * w * str);

      col = mix(col, vortexTint(i, seed), w * str * 0.65);
    }

    return col;
  }

  if (kind == 5) { // Star
    float coarse = fbm01(p * 5.0, seed ^ 0x57A4u, 5, 2.0, 0.55);
    float fine   = fbm01(p * 18.0, seed ^ 0xF1AEu, 3, 2.4, 0.55);
    float t = clamp01(0.65 * coarse + 0.35 * fine);

    vec3 c0 = vec3(1.00, 0.70, 0.18);
    vec3 c1 = vec3(1.00, 0.92, 0.58);
    if ((paletteIx & 1) != 0) {
      c0 = vec3(1.00, 0.84, 0.52);
      c1 = vec3(1.00, 0.98, 0.86);
    }

    vec3 col = mix(c0, c1, t);

    float spot = fbm01(p * 9.0, seed ^ 0xABCDEFu, 4, 2.0, 0.5);
    float spotMask = smoothstep01(0.78, 0.92, spot);
    col = mix(col, col * 0.82, spotMask * 0.25);

    return col;
  }

  if (kind == 6) { // Clouds
    outAlpha = surfaceCloudAlpha(seed, p, lat, lon);
    return vec3(0.94, 0.97, 1.00);
  }

  return vec3(0.6);
}

float kindScaleForNormal(int kind) {
  if (kind == 0) return 1.25; // Rocky
  if (kind == 1) return 1.05; // Desert
  if (kind == 2) return 0.85; // Ocean
  if (kind == 3) return 1.10; // Ice
  if (kind == 4) return 0.22; // Gas
  if (kind == 5) return 0.34; // Star
  if (kind == 6) return 0.60; // Clouds
  return 1.0;
}

void main() {
  // UV is in [0,1] across the viewport-filling part of the triangle.
  float u = vUv.x;
  float v = vUv.y;

  // Rotation so systems don't all align on the seam.
  uint seed = seedMix();
  uint rotBits = (uint(uSeedLo) >> 8u) & 0xFFFFu;
  float rot = float(rotBits) / 65535.0 * (2.0 * PI);

  float lat = (v - 0.5) * PI;
  float lon = u * (2.0 * PI) - PI + rot;
  vec3 p = latLonToSphere(lat, lon);

  if (uMode == 0) {
    float a;
    vec3 col = surfaceAlbedo(uKind, seed, p, lat, lon, a);
    FragColor = vec4(col, a);
    return;
  }

  // Normal map.
  float du = uInvSize.x;
  float dv = uInvSize.y;

  // Physical distances on a unit sphere for one texel step.
  float dLon = du * (2.0 * PI);
  float dLat = dv * PI;
  float cosLat = cos(lat);
  float dx = max(1.0e-6, abs(cosLat) * dLon);
  float dy = max(1.0e-6, dLat);

  // Wrap in U, clamp in V.
  float uL = fract(u - du);
  float uR = fract(u + du);
  float vD = clamp01(v - dv);
  float vU = clamp01(v + dv);

  float latL = (v - 0.5) * PI;
  float latD = (vD - 0.5) * PI;
  float latU = (vU - 0.5) * PI;

  float lonL = uL * (2.0 * PI) - PI + rot;
  float lonR = uR * (2.0 * PI) - PI + rot;
  float lonC = u * (2.0 * PI) - PI + rot;

  vec3 pL = latLonToSphere(latL, lonL);
  vec3 pR = latLonToSphere(latL, lonR);
  vec3 pD = latLonToSphere(latD, lonC);
  vec3 pU = latLonToSphere(latU, lonC);

  float hL = surfaceHeight(uKind, seed, pL, latL, lonL);
  float hR = surfaceHeight(uKind, seed, pR, latL, lonR);
  float hD = surfaceHeight(uKind, seed, pD, latD, lonC);
  float hU = surfaceHeight(uKind, seed, pU, latU, lonC);

  float dHdx = (hR - hL) / (2.0 * dx);
  float dHdy = (hU - hD) / (2.0 * dy);

  float k = kindScaleForNormal(uKind);
  vec3 nrm = normalize(vec3(-dHdx * k, -dHdy * k, 1.0));
  FragColor = vec4(nrm * 0.5 + 0.5, 1.0);
}

)GLSL";

static core::i32 bitcastU32ToI32(core::u32 u) {
  core::i32 s = 0;
  static_assert(sizeof(s) == sizeof(u));
  std::memcpy(&s, &u, sizeof(u));
  return s;
}

} // namespace

GpuSurfaceCache::~GpuSurfaceCache() { destroy(); }

void GpuSurfaceCache::destroy() {
  clear();

  if (vao_) {
    gl::DeleteVertexArrays(1, &vao_);
    vao_ = 0;
  }
  if (fbo_) {
    gl::DeleteFramebuffers(1, &fbo_);
    fbo_ = 0;
  }

  shader_ = ShaderProgram{};
  inited_ = false;
}

bool GpuSurfaceCache::init(std::string* outError) {
  if (inited_) return true;

  std::string err;
  if (!shader_.build(kVS, kFS, &err)) {
    if (outError) *outError = err;
    destroy();
    return false;
  }

  gl::GenFramebuffers(1, &fbo_);
  gl::GenVertexArrays(1, &vao_);

  // Validate that the FBO object exists.
  gl::BindFramebuffer(GL_FRAMEBUFFER, fbo_);
  gl::BindFramebuffer(GL_FRAMEBUFFER, 0);

  inited_ = true;
  return true;
}

void GpuSurfaceCache::clear() {
  albedoCache_.clear();
  normalCache_.clear();
  tick_ = 0;
}

core::u64 GpuSurfaceCache::makeKey(const char* tag, SurfaceKind kind, core::u64 seed, int widthPx) const {
  core::u64 h = core::fnv1a64(tag);
  h = core::hashCombine(h, (core::u64)(core::u8)kind);
  h = core::hashCombine(h, seed);
  h = core::hashCombine(h, (core::u64)(core::i64)widthPx);
  return h;
}

void GpuSurfaceCache::evictIfNeeded(std::unordered_map<core::u64, Entry>& cache) {
  if (maxEntries_ == 0) {
    cache.clear();
    return;
  }
  while (cache.size() > maxEntries_) {
    core::u64 oldestKey = 0;
    core::u64 oldestTick = (core::u64)-1;
    for (const auto& kv : cache) {
      if (kv.second.lastUseTick < oldestTick) {
        oldestTick = kv.second.lastUseTick;
        oldestKey = kv.first;
      }
    }
    if (oldestKey != 0) cache.erase(oldestKey);
    else break;
  }
}

bool GpuSurfaceCache::render(Texture2D& out, SurfaceKind kind, core::u64 seed, int w, int h, int mode) {
  if (!inited_) return false;
  if (w <= 0 || h <= 0) return false;

  // Save a small set of GL state we touch.
  GLint prevFbo = 0;
  GLint prevVao = 0;
  GLint prevProg = 0;
  GLint prevViewport[4] = {0,0,0,0};
  ::glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
  ::glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &prevVao);
  ::glGetIntegerv(GL_CURRENT_PROGRAM, &prevProg);
  ::glGetIntegerv(GL_VIEWPORT, prevViewport);

  const GLboolean prevDepth = glIsEnabled(GL_DEPTH_TEST);
  const GLboolean prevBlend = glIsEnabled(GL_BLEND);
  const GLboolean prevCull  = glIsEnabled(GL_CULL_FACE);
  const GLboolean prevScissor = glIsEnabled(GL_SCISSOR_TEST);

  // Allocate output texture.
  out = Texture2D{};
  out.allocateRGBA(w, h,
                   /*generateMips=*/true,
                   /*nearestFilter=*/false,
                   /*clampToEdge=*/false);

  gl::BindFramebuffer(GL_FRAMEBUFFER, fbo_);
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, out.handle(), 0);
  ::glDrawBuffer(GL_COLOR_ATTACHMENT0);

  const GLenum status = gl::CheckFramebufferStatus(GL_FRAMEBUFFER);
  if (status != GL_FRAMEBUFFER_COMPLETE) {
    // Restore minimal state and bail.
    gl::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    gl::BindVertexArray((GLuint)prevVao);
    gl::UseProgram((GLuint)prevProg);
    glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);
    return false;
  }

  glViewport(0, 0, w, h);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  glDisable(GL_CULL_FACE);
  glDisable(GL_SCISSOR_TEST);

  // Clear (mainly to avoid driver warnings on some setups).
  ::glClearColor(0, 0, 0, 0);
  ::glClear(GL_COLOR_BUFFER_BIT);

  shader_.bind();
  shader_.setUniform1i("uKind", (int)(core::u8)kind);
  shader_.setUniform1i("uMode", mode);

  const core::u32 seedLo = (core::u32)(seed & 0xFFFFFFFFull);
  const core::u32 seedHi = (core::u32)((seed >> 32) & 0xFFFFFFFFull);
  shader_.setUniform1i("uSeedLo", (int)bitcastU32ToI32(seedLo));
  shader_.setUniform1i("uSeedHi", (int)bitcastU32ToI32(seedHi));

  shader_.setUniform2f("uInvSize", 1.0f / (float)w, 1.0f / (float)h);

  gl::BindVertexArray(vao_);
  ::glDrawArrays(GL_TRIANGLES, 0, 3);

  // Build mip chain after the base level is rendered.
  gl::BindTexture(GL_TEXTURE_2D, out.handle());
  gl::GenerateMipmap(GL_TEXTURE_2D);

  // Detach to keep FBO clean.
  gl::FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, 0);

  // Restore previous state.
  gl::BindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
  gl::BindVertexArray((GLuint)prevVao);
  gl::UseProgram((GLuint)prevProg);
  glViewport(prevViewport[0], prevViewport[1], prevViewport[2], prevViewport[3]);

  if (prevDepth) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
  if (prevBlend) glEnable(GL_BLEND); else glDisable(GL_BLEND);
  if (prevCull) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);
  if (prevScissor) glEnable(GL_SCISSOR_TEST); else glDisable(GL_SCISSOR_TEST);

  return true;
}

const Texture2D& GpuSurfaceCache::albedo(SurfaceKind kind, core::u64 seed, int widthPx) {
  ++tick_;
  widthPx = std::clamp(widthPx, 64, 2048);
  const int h = std::max(2, widthPx / 2);

  const core::u64 k = makeKey("gpu_surface_albedo", kind, seed, widthPx);
  if (auto it = albedoCache_.find(k); it != albedoCache_.end()) {
    it->second.lastUseTick = tick_;
    return it->second.tex;
  }

  Entry e{};
  (void)render(e.tex, kind, seed, widthPx, h, /*mode=*/0);
  e.lastUseTick = tick_;

  auto [it, inserted] = albedoCache_.emplace(k, std::move(e));
  (void)inserted;
  evictIfNeeded(albedoCache_);
  return it->second.tex;
}

const Texture2D& GpuSurfaceCache::normal(SurfaceKind kind, core::u64 seed, int widthPx) {
  ++tick_;
  widthPx = std::clamp(widthPx, 64, 2048);
  const int h = std::max(2, widthPx / 2);

  const core::u64 k = makeKey("gpu_surface_normal", kind, seed, widthPx);
  if (auto it = normalCache_.find(k); it != normalCache_.end()) {
    it->second.lastUseTick = tick_;
    return it->second.tex;
  }

  Entry e{};
  (void)render(e.tex, kind, seed, widthPx, h, /*mode=*/1);
  e.lastUseTick = tick_;

  auto [it, inserted] = normalCache_.emplace(k, std::move(e));
  (void)inserted;
  evictIfNeeded(normalCache_);
  return it->second.tex;
}

} // namespace stellar::render
