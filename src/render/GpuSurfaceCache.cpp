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
static const char* kFS_0 = R"GLSL(
#version 330 core

in vec2 vUv;

uniform int uKind;    // render::SurfaceKind (0..)
uniform int uMode;    // 0=albedo, 1=normal
uniform int uSeedLo;  // bitwise seed (lower 32)
uniform int uSeedHi;  // bitwise seed (upper 32)
uniform vec2 uInvSize; // (1/w, 1/h)

out vec4 FragColor;

const float PI = 3.1415926535897932384626433832795;

float clamp01(float x) { return clamp(x, 0.0, 1.0); }

// Quintic fade (Perlin-style): 6t^5 - 15t^4 + 10t^3
float fade(float t) {
  return t*t*t*(t*(t*6.0 - 15.0) + 10.0);
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

float smoothstep01(float e0, float e1, float x) {
  float t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3.0 - 2.0 * t);
}


// -----------------------------------------------------------------------------
// Crater field synthesis (Worley-style nearest feature on p*freq).
//
// This is deliberately tiny and deterministic so the GPU surface cache can
// produce "impact" detail without any lookup textures.

float u01(uint h) { return float(h) / 4294967295.0; }

struct CraterSample {
  float h;
  float bowl;
  float rim;
  float ejecta;
};

void worleyNearest(vec3 x, uint seed, out float outDist, out vec3 outFeature, out uint outId) {
  ivec3 base = ivec3(floor(x));
  float minD = 1.0e9;
  uint bestId = 0u;
  vec3 bestF = vec3(0.0);

  for (int dz = -1; dz <= 1; ++dz) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        ivec3 c = base + ivec3(dx, dy, dz);
        uvec3 uc = uvec3(c);

        uint id = pcg_hash(uc.x + pcg_hash(uc.y + pcg_hash(uc.z + seed)));

        float rx = u01(pcg_hash(id ^ 0xA531u));
        float ry = u01(pcg_hash(id ^ 0xC0DEu));
        float rz = u01(pcg_hash(id ^ 0xBEEFu));

        vec3 f = vec3(c) + vec3(rx, ry, rz);
        float d = length(x - f);

        if (d < minD) {
          minD = d;
          bestId = id;
          bestF = f;
        }
      }
    }
  }

  outDist = minD;
  outId = bestId;
  outFeature = bestF;
}

CraterSample craterLayer(vec3 p, uint seed, float freq, float amp, bool withRays) {
  vec3 x = p * freq;

  float dist;
  vec3 feature;
  uint id;
  worleyNearest(x, seed, dist, feature, id);

  float radius = mix(0.28, 0.58, u01(pcg_hash(id ^ 0xC0DEC0DEu)));
  float depthScale = mix(0.75, 1.40, u01(pcg_hash(id ^ 0xBEEF1234u)));

  float d = dist / max(1.0e-6, radius);

  float t = clamp01(1.0 - d);
  float bowl = t * t;

  float rim = smoothstep01(0.86, 1.0, d) * (1.0 - smoothstep01(1.0, 1.14, d));

  float ejectaRing = smoothstep01(1.0, 1.08, d) * (1.0 - smoothstep01(1.08, 1.70, d));
  float ejecta = ejectaRing;

  if (withRays && ejectaRing > 1.0e-5) {
    vec3 cpos = feature / freq;
    vec3 cdir = (length(cpos) > 1.0e-6) ? normalize(cpos) : p;

    vec3 up = (abs(cdir.y) > 0.85) ? vec3(1.0, 0.0, 0.0) : vec3(0.0, 1.0, 0.0);
    vec3 t0 = normalize(cross(up, cdir));
    vec3 t1 = cross(cdir, t0);

    vec3 v = p - cdir * dot(p, cdir);
    float vLen = length(v);
    if (vLen > 1.0e-6) {
      v /= vLen;
      float ang = atan(dot(v, t1), dot(v, t0));

      float rayCount = float(6 + int(pcg_hash(id ^ 0x51DEu) % 15u));
      float phase = u01(pcg_hash(id ^ 0xA11CE5u)) * (2.0 * PI);

      float rv = 0.5 + 0.5 * cos(ang * rayCount + phase);
      rv = pow(clamp01(rv), 5.0);
      ejecta *= mix(0.35, 1.0, rv);
    }
  }

  CraterSample cs;
  cs.h = (-amp * depthScale * bowl) + (amp * 0.45 * rim) + (amp * 0.10 * ejectaRing);
  cs.bowl = bowl;
  cs.rim = rim;
  cs.ejecta = ejecta;
  return cs;
}

CraterSample craterField(vec3 p, uint seed) {
  CraterSample a = craterLayer(p, seed ^ 0xC7A73u, 2.2, 0.10, true);
  CraterSample b = craterLayer(p, seed ^ 0xF00Du,  7.5, 0.05, false);
  CraterSample c = craterLayer(p, seed ^ 0xB16D0u, 16.0, 0.02, false);

  CraterSample cs;
  cs.h = a.h + b.h + c.h;
  cs.bowl = max(a.bowl, max(b.bowl * 0.75, c.bowl * 0.50));
  cs.rim = max(a.rim, max(b.rim * 0.75, c.rim * 0.50));
  cs.ejecta = max(a.ejecta, max(b.ejecta * 0.65, c.ejecta * 0.45));
  return cs;
}


vec3 latLonToSphere(float lat, float lon) {
  float cosLat = cos(lat);
  float sinLat = sin(lat);
  float cosLon = cos(lon);
  float sinLon = sin(lon);
  return vec3(cosLat * cosLon, sinLat, cosLat * sinLon);
}

float surfaceCloudAlpha(uint seed, vec3 p, float lat, float lon) {
  // Cloud density based on layered FBM.
  float n1 = fbm3D(p * 3.6, seed ^ 0xC10Du, 6, 2.0, 0.5);
  float n2 = fbm3D(p * 10.5 + vec3(0.0, lat * 0.35, lon * 0.15), seed ^ 0xBEEFu, 4, 2.1, 0.55);
  float n = 0.65 * n1 + 0.35 * n2;
  float d = smoothstep01(0.55, 0.80, n);
  // Subtle latitude banding.
  float band = 0.5 + 0.5 * sin(lat * 8.0 + (n2 - 0.5) * 3.0);
  d *= 0.80 + 0.20 * band;
  return clamp01(d);
}

float surfaceHeight(int kind, uint seed, vec3 p, float lat, float lon) {
  // Height fields are intentionally approximate; only relative variation matters.
  if (kind == 0) { // Rocky
    float e0 = fbm3D(p * 2.8, seed ^ 0xA11CE5u, 6, 2.0, 0.5);
    float e  = clamp01((e0 - 0.25) * 1.35);
    float r0 = fbm3D(p * 7.5, seed ^ 0xBADC0FFEu, 4, 2.1, 0.55);
    float ridged = 1.0 - abs(r0 * 2.0 - 1.0);
    float dust = fbm3D(p * 15.0, seed ^ 0xC0FFEEu, 3, 2.4, 0.55);
    float h = clamp01(0.70 * e + 0.30 * ridged);
    h = clamp01(h * (0.92 + 0.10 * dust));

    CraterSample cs = craterField(p, seed ^ 0xC0DEF00Du);
    h = clamp01(h + cs.h);

    return h;
  }
  if (kind == 1) { // Desert
    float e0 = fbm3D(p * 2.2, seed ^ 0xD35E7u, 6, 2.0, 0.5);
    float e  = clamp01((e0 - 0.22) * 1.30);
    float warp = fbm3D(p * 3.4, seed ^ 0x51DEu, 4, 2.1, 0.55);
    float band = 0.5 + 0.5 * sin((lat * 16.0) + (lon * 2.2) + (warp - 0.5) * 3.2);
    float rock = fbm3D(p * 6.0, seed ^ 0x0BADC0DEu, 4, 2.1, 0.52);
    float rockMask = smoothstep01(0.70, 0.88, rock);
    float h = 0.78 * e + 0.22 * band;
    h += rockMask * 0.12;

    CraterSample cs = craterField(p, seed ^ 0xC0DEF00Du);
    h = clamp01(h + cs.h * 0.85);

    return h;
  }
  if (kind == 2) { // Ocean
    float n = fbm3D(p * 2.1, seed ^ 0x0CE4u, 6, 2.0, 0.5);
    float land = smoothstep01(0.48, 0.55, n);
    float e0 = fbm3D(p * 6.5, seed ^ 0xE1E7u, 5, 2.1, 0.55);
    float elev = clamp01((e0 - 0.30) * 1.55);
    float wave = fbm3D(p * 18.0, seed ^ 0xA57EA11Cu, 3, 2.2, 0.55);
    float hOcean = 0.03 * (wave - 0.5);
    float hLand = 0.25 + 0.75 * elev;
    return clamp01(mix(hOcean, hLand, land));
  }
  if (kind == 3) { // Ice
    float n = fbm3D(p * 3.3, seed ^ 0x1CEu, 6, 2.0, 0.5);
    float crack = fbm3D(p * 20.0, seed ^ 0xC24Cu, 3, 2.4, 0.55);
    float crackMask = smoothstep01(0.76, 0.93, crack);
    float h = clamp01((n - 0.22) * 1.35);
    h = clamp01(h + 0.22 * crackMask);

    CraterSample cs = craterField(p, seed ^ 0xC0DEF00Du);
    h = clamp01(h + cs.h * 0.65);

    return h;
  }
  if (kind == 4) { // GasGiant
    float warp = fbm3D(p * 2.6, seed ^ 0x6A59u, 5, 2.0, 0.55);
    float band = 0.5 + 0.5 * sin((lat * 22.0) + (warp - 0.5) * 4.2);
    float turb = fbm3D(p * 7.0, seed ^ 0x7B3Bu, 4, 2.1, 0.55);
    float h = 0.5 + 0.10 * (band - 0.5) + 0.08 * (turb - 0.5);
    return clamp01(h);
  }
  if (kind == 5) { // Star
    float coarse = fbm3D(p * 5.0, seed ^ 0x57A4u, 5, 2.0, 0.55);
    float fine   = fbm3D(p * 18.0, seed ^ 0xF1AEu, 3, 2.4, 0.55);
    float h = 0.55 + 0.55 * coarse + 0.20 * (fine - 0.5);
    return clamp01(h);
  }
  if (kind == 6) { // Clouds
    return surfaceCloudAlpha(seed, p, lat, lon);
  }
  if (kind == 7) { // CityLights
    // Emissive-only map; keep normals essentially flat by returning a constant.
    return 0.5;
  }
  return 0.5;
}

vec3 surfaceAlbedo(int kind, uint seed, vec3 p, float lat, float lon) {
  if (kind == 0) { // Rocky
    float e0 = fbm3D(p * 2.8, seed ^ 0xA11CE5u, 6, 2.0, 0.5);
    float e  = clamp01((e0 - 0.25) * 1.35);
    float r0 = fbm3D(p * 7.5, seed ^ 0xBADC0FFEu, 4, 2.1, 0.55);
    float ridged = 1.0 - abs(r0 * 2.0 - 1.0);
    float dust = fbm3D(p * 15.0, seed ^ 0xC0FFEEu, 3, 2.4, 0.55);
    vec3 baseLo = vec3(0.18, 0.16, 0.15);
    vec3 baseHi = vec3(0.62, 0.56, 0.48);
    vec3 col = mix(baseLo, baseHi, e);
    col = mix(col, col * 0.65, (1.0 - ridged) * 0.35);
    col *= (0.85 + 0.28 * dust);
    float speck = fbm3D(p * 26.0, seed ^ 0xFEEDBEEFu, 2, 2.5, 0.6);
    float speckMask = smoothstep01(0.82, 0.96, speck);
    col = mix(col, vec3(0.88, 0.84, 0.78), speckMask * 0.32);

    CraterSample cs = craterField(p, seed ^ 0xC0DEF00Du);
    col = mix(col, col * 0.55, cs.bowl * 0.25);
    col = mix(col, vec3(0.92, 0.90, 0.86), cs.rim * 0.35);
    col = mix(col, vec3(0.88, 0.84, 0.78), cs.ejecta * 0.18);

    return col;
  }
  if (kind == 1) { // Desert
    float e0 = fbm3D(p * 2.2, seed ^ 0xD35E7u, 6, 2.0, 0.5);
    float e  = clamp01((e0 - 0.22) * 1.30);
    vec3 sandLo = vec3(0.62, 0.52, 0.28);
    vec3 sandHi = vec3(0.92, 0.86, 0.50);
    vec3 col = mix(sandLo, sandHi, e);
    float warp = fbm3D(p * 3.4, seed ^ 0x51DEu, 4, 2.1, 0.55);
    float band = 0.5 + 0.5 * sin((lat * 16.0) + (lon * 2.2) + (warp - 0.5) * 3.2);
    col *= (0.90 + 0.18 * band);
    float rock = fbm3D(p * 6.0, seed ^ 0x0BADC0DEu, 4, 2.1, 0.52);
    float rockMask = smoothstep01(0.70, 0.88, rock);
    col = mix(col, vec3(0.36, 0.30, 0.22), rockMask * 0.55);

    CraterSample cs = craterField(p, seed ^ 0xC0DEF00Du);
    col = mix(col, col * 0.75, cs.bowl * 0.22);
    col = mix(col, vec3(0.97, 0.93, 0.78), cs.rim * 0.28);
    col = mix(col, vec3(0.96, 0.92, 0.84), cs.ejecta * 0.22);

    return col;
  }
  if (kind == 2) { // Ocean
    float n = fbm3D(p * 2.1, seed ^ 0x0CE4u, 6, 2.0, 0.5);
    float land = smoothstep01(0.48, 0.55, n);
    float e0 = fbm3D(p * 6.5, seed ^ 0xE1E7u, 5, 2.1, 0.55);
    float elev = clamp01((e0 - 0.30) * 1.55);
    float coast = clamp01(1.0 - abs(n - 0.52) * 18.0);
    vec3 waterDeep = vec3(0.03, 0.08, 0.22);
    vec3 waterShallow = vec3(0.07, 0.22, 0.45);
    vec3 water = mix(waterDeep, waterShallow, coast);
    vec3 landLo = vec3(0.12, 0.30, 0.16);
    vec3 landHi = vec3(0.52, 0.46, 0.26);
    vec3 landCol = mix(landLo, landHi, elev);
    vec3 col = mix(water, landCol, land);

    float cap = smoothstep01(0.62, 0.90, abs(sin(lat)));
    if (cap > 1.0e-6) {
      float capNoise = fbm3D(p * 9.0, seed ^ 0xC4F5u, 3, 2.2, 0.55);
      float capMask = cap * (0.65 + 0.35 * capNoise);
      col = mix(col, vec3(0.90, 0.95, 1.00), capMask);
    }
    float cloudN = fbm3D(p * 4.0, seed ^ 0xC10Du, 5, 2.0, 0.5);
    float cloud = smoothstep01(0.70, 0.92, cloudN) * (1.0 - land * 0.25);
    col = mix(col, vec3(1.0), cloud * 0.10);
    return col;
  }
  if (kind == 3) { // Ice
    float n = fbm3D(p * 3.3, seed ^ 0x1CEu, 6, 2.0, 0.5);
    vec3 iceLo = vec3(0.70, 0.84, 0.94);
    vec3 iceHi = vec3(0.95, 0.99, 1.00);
    vec3 col = mix(iceLo, iceHi, clamp01((n - 0.22) * 1.35));
    float crack = fbm3D(p * 20.0, seed ^ 0xC24Cu, 3, 2.4, 0.55);
    float crackMask = smoothstep01(0.76, 0.93, crack);
    col = mix(col, vec3(0.25, 0.35, 0.40), crackMask * 0.22);
    col *= (0.92 + 0.16 * fbm3D(p * 12.0, seed ^ 0xF00Du, 3, 2.1, 0.55));

    CraterSample cs = craterField(p, seed ^ 0xC0DEF00Du);
    float bowl = cs.bowl * 0.65;
    float rim = cs.rim * 0.65;
    float ej = cs.ejecta * 0.55;
    col = mix(col, vec3(0.62, 0.74, 0.86), bowl * 0.18);
    col = mix(col, vec3(0.98, 1.00, 1.00), rim * 0.25);
    col = mix(col, vec3(0.90, 0.96, 1.00), ej * 0.16);

    return col;
  }
  if (kind == 4) { // GasGiant
    float warp = fbm3D(p * 2.6, seed ^ 0x6A59u, 5, 2.0, 0.55);
    float band = 0.5 + 0.5 * sin((lat * 22.0) + (warp - 0.5) * 4.2);
    float turb = fbm3D(p * 7.0, seed ^ 0x7B3Bu, 4, 2.1, 0.55);
    vec3 a = vec3(0.86, 0.68, 0.46);
    vec3 b = vec3(0.98, 0.90, 0.74);
    vec3 c = vec3(0.66, 0.44, 0.32);
    vec3 col = mix(a, b, band);
    col = mix(col, c, (turb - 0.5) * 0.35);
    // Occasional storm spots.
    float spot = fbm3D(p * 14.0, seed ^ 0xA11CEu, 3, 2.2, 0.55);
    float spotMask = smoothstep01(0.80, 0.93, spot);
    col = mix(col, vec3(1.0, 0.92, 0.80), spotMask * 0.25);
    return clamp(col, 0.0, 1.0);
  }
  if (kind == 5) { // Star
    float coarse = fbm3D(p * 5.0, seed ^ 0x57A4u, 5, 2.0, 0.55);
    float fine   = fbm3D(p * 18.0, seed ^ 0xF1AEu, 3, 2.4, 0.55);
    float grain  = fbm3D(p * 42.0, seed ^ 0xABCDEFu, 2, 2.4, 0.55);
    float h = clamp01(0.55 + 0.55 * coarse + 0.20 * (fine - 0.5));
    vec3 base = mix(vec3(0.95, 0.62, 0.22), vec3(1.0, 0.92, 0.70), clamp01(h));
    base *= (1.15 + 0.20 * (grain - 0.5));
    // Hot granules.
    float hot = smoothstep01(0.78, 0.95, fine);
    base = mix(base, vec3(1.35, 1.10, 0.65), hot * 0.35);
    return base;
  }
  if (kind == 6) { // Clouds (RGB near-white, alpha density)
    return vec3(0.94, 0.97, 1.00);
  }
  if (kind == 7) { // CityLights
    float cont = fbm3D(p * 2.1, seed ^ 0x0CE4u, 6, 2.0, 0.5);
    float land = smoothstep01(0.48, 0.55, cont);
    float coast = clamp01(1.0 - abs(cont - 0.52) * 18.0);

    float latW = pow(abs(cos(lat)), 0.30);

    float pop0 = fbm3D(p * 6.0, seed ^ 0xC1717u, 5, 2.1, 0.55);
    float pop1 = fbm3D(p * 14.0 + vec3(12.3, -7.1, 3.7), seed ^ 0xB16D0u, 4, 2.2, 0.55);
    float pop = clamp01(0.65 * pop0 + 0.35 * pop1);

    float cluster = smoothstep01(0.72, 0.92, pop);

    float r0 = fbm3D(p * 22.0 + vec3(4.1, -9.7, 2.6), seed ^ 0x70ADu, 3, 2.3, 0.55);
    float rid = 1.0 - abs(r0 * 2.0 - 1.0);
)GLSL";

// NOTE: MSVC has a relatively small per-literal size limit; keep large GLSL sources
// split into multiple literals and concatenate at runtime.
static const char* kFS_1 = R"GLSL(    float roads = smoothstep01(0.78, 0.94, rid) * (0.30 + 0.70 * cluster);

    float mega0 = fbm3D(p * 1.6 + vec3(19.7, -5.2, 11.3), seed ^ 0x31E6Au, 4, 2.0, 0.5);
    float mega = smoothstep01(0.78, 0.95, mega0);
    mega *= mega;

    float d = land * latW * (0.60 * cluster + 0.25 * roads + 0.55 * mega * (0.35 + 0.65 * coast));
    d *= (0.35 + 0.65 * coast);

    float hole = fbm3D(p * 4.5 + vec3(-7.1, 2.2, -3.3), seed ^ 0x4011Eu, 4, 2.0, 0.5);
    d *= (1.0 - 0.55 * smoothstep01(0.55, 0.80, hole));

    d = smoothstep01(0.08, 0.28, d);
    d = pow(clamp01(d), 1.25);

    float spark = fbm3D(p * 55.0, seed ^ 0x5A2F01u, 2, 2.5, 0.6);
    float sparkMask = smoothstep01(0.74, 0.93, spark);
    d = clamp01(d * (1.0 + 0.35 * sparkMask));

    float warmT = clamp01(0.65 + 0.35 * float((seed >> 16u) & 0xFFu) / 255.0);
    vec3 cool = vec3(0.72, 0.88, 1.00);
    vec3 warm = vec3(1.00, 0.84, 0.58);
    vec3 tint = mix(cool, warm, warmT);

    return tint * d;
  }
  return vec3(0.6);
}

float kindScaleForNormal(int kind) {
  if (kind == 0) return 1.35; // Rocky
  if (kind == 1) return 1.05; // Desert
  if (kind == 2) return 0.85; // Ocean
  if (kind == 3) return 1.05; // Ice
  if (kind == 4) return 0.20; // Gas
  if (kind == 5) return 0.30; // Star
  if (kind == 6) return 0.55; // Clouds
  if (kind == 7) return 0.0; // CityLights
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
    vec3 col = surfaceAlbedo(uKind, seed, p, lat, lon);
    float a = 1.0;
    if (uKind == 6) {
      a = surfaceCloudAlpha(seed, p, lat, lon);
    }
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

static const std::string& fullFragmentSource_() {
  static const std::string s = std::string(kFS_0) + std::string(kFS_1);
  return s;
}


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
  if (!shader_.build(kVS, fullFragmentSource_(), &err)) {
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
