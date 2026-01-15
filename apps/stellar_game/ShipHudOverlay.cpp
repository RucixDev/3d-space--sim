#include "ShipHudOverlay.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cctype>

namespace stellar::ui {

namespace {

static float clamp01(float v) {
  return std::clamp(v, 0.0f, 1.0f);
}

static ImU32 toU32(const Color4f& c, float alphaScale = 1.0f) {
  auto b = [](float v) -> int {
    const float cl = std::clamp(v, 0.0f, 1.0f);
    return (int)std::lround(cl * 255.0f);
  };

  const float a = clamp01(c.a * alphaScale);
  return IM_COL32(b(c.r), b(c.g), b(c.b), b(a));
}

static float lerpf(float a, float b, float t) {
  return a + (b - a) * t;
}

// --- Procedural segment-display font (14-seg-ish) --------------------------
//
// HUD value readouts look far more "ship-like" when rendered as a segmented display
// rather than standard UI text. This is not a strict hardware-accurate encoding;
// it's a stylized vector approximation designed for readability at small sizes.
//
// Segment naming follows the common 14-seg convention: a-f + g1/g2 + h-k + l/m.
// We also support a DP (decimal point) segment.

enum SegBits : unsigned int {
  SA  = 1u << 0,
  SB  = 1u << 1,
  SC  = 1u << 2,
  SD  = 1u << 3,
  SE  = 1u << 4,
  SF  = 1u << 5,
  SG1 = 1u << 6,
  SG2 = 1u << 7,
  SH  = 1u << 8,
  SI  = 1u << 9,
  SJ  = 1u << 10,
  SK  = 1u << 11,
  SL  = 1u << 12,
  SM  = 1u << 13,
  SDP = 1u << 14,
};

static float segCharWidth(float h) {
  return h * 0.62f;
}

static float segCharSpacing(float h, float uiScale) {
  return std::max(1.0f * uiScale, h * 0.12f);
}

static ImVec2 segTextSize(const char* s, float h, float spacing) {
  const int n = (s && *s) ? (int)std::strlen(s) : 0;
  const float w = (n > 0) ? (n * segCharWidth(h) + (n - 1) * spacing) : 0.0f;
  return ImVec2(w, h);
}

static unsigned int segMaskForChar(char ch) {
  const unsigned char uc = (unsigned char)ch;
  const char c = (char)std::toupper((int)uc);

  // Digits (classic 7-seg style using a-f + g1/g2).
  switch (c) {
    case '0': return SA | SB | SC | SD | SE | SF;
    case '1': return SB | SC;
    case '2': return SA | SB | SG1 | SG2 | SE | SD;
    case '3': return SA | SB | SC | SG1 | SG2 | SD;
    case '4': return SF | SG1 | SG2 | SB | SC;
    case '5': return SA | SF | SG1 | SG2 | SC | SD;
    case '6': return SA | SF | SE | SD | SC | SG1 | SG2;
    case '7': return SA | SB | SC;
    case '8': return SA | SB | SC | SD | SE | SF | SG1 | SG2;
    case '9': return SA | SB | SC | SD | SF | SG1 | SG2;

    // High-frequency status letters used by the Ship HUD (CHG/JMP/READY/CD).
    case 'A': return SA | SB | SC | SE | SF | SG1 | SG2;
    case 'B': return SA | SB | SC | SD | SE | SF | SG1 | SG2;
    case 'C': return SA | SD | SE | SF;
    case 'D': return SB | SC | SD | SE | SF;
    case 'E': return SA | SD | SE | SF | SG1 | SG2;
    case 'F': return SA | SE | SF | SG1 | SG2;
    case 'G': return SA | SD | SE | SF | SC | SG2;
    case 'H': return SB | SC | SE | SF | SG1 | SG2;
    case 'I': return SA | SD | SL | SM;
    case 'J': return SB | SC | SD | SE;
    case 'K': return SE | SF | SH | SJ | SI | SK;
    case 'L': return SD | SE | SF;
    case 'M': return SA | SB | SC | SE | SF | SH | SI;
    case 'N': return SB | SC | SE | SF | SH | SK;
    case 'O': return SA | SB | SC | SD | SE | SF;
    case 'P': return SA | SB | SE | SF | SG1 | SG2;
    case 'Q': return SA | SB | SC | SD | SE | SF | SK;
    case 'R': return SA | SB | SE | SF | SG1 | SG2 | SK;
    case 'S': return SA | SF | SG1 | SG2 | SC | SD;
    case 'T': return SA | SL | SM;
    case 'U': return SB | SC | SD | SE | SF;
    case 'V': return SE | SF | SJ | SK;
    case 'W': return SB | SC | SD | SE | SF | SJ | SK;
    case 'X': return SH | SI | SJ | SK;
    case 'Y': return SH | SI | SM | SD;
    case 'Z': return SA | SI | SJ | SD;

    // Punctuation handled as special cases in the renderer.
    default: break;
  }

  return 0u;
}

static ImU32 scaleAlpha(ImU32 c, float a) {
  const unsigned int ca = (c >> 24) & 0xFFu;
  const unsigned int na = (unsigned int)std::lround((double)ca * std::clamp(a, 0.0f, 1.0f));
  return (c & 0x00FFFFFFu) | (na << 24);
}

static void segLine(ImDrawList* draw,
                    ImVec2 a,
                    ImVec2 b,
                    ImU32 col,
                    float thickness,
                    float dropoutP,
                    core::SplitMix64* rng) {
  if (dropoutP > 0.0f && rng && rng->chance((double)dropoutP)) return;
  draw->AddLine(a, b, col, thickness);
}

static void segDot(ImDrawList* draw,
                   ImVec2 p,
                   ImU32 col,
                   float r,
                   float dropoutP,
                   core::SplitMix64* rng) {
  if (dropoutP > 0.0f && rng && rng->chance((double)dropoutP)) return;
  draw->AddCircleFilled(p, r, col);
}

static void drawSegChar(ImDrawList* draw,
                        ImVec2 pos,
                        float h,
                        char ch,
                        ImU32 col,
                        float thickness,
                        float dropoutP,
                        core::SplitMix64* rng) {
  const float w = segCharWidth(h);
  const float pad = std::max(1.0f, h * 0.12f);

  const float x0 = pos.x;
  const float y0 = pos.y;
  const float x1 = pos.x + w;
  const float y1 = pos.y + h;
  const float xm = (x0 + x1) * 0.5f;
  const float ym = (y0 + y1) * 0.5f;

  const float xl = x0 + pad;
  const float xr = x1 - pad;
  const float yt = y0 + pad;
  const float yb = y1 - pad;

  // Middle segments
  const float ymid = ym;

  // Special punctuation.
  if (ch == ' ' || ch == '\t') return;
  if (ch == '.') {
    segDot(draw, ImVec2(xr, yb), col, std::max(1.2f, thickness * 0.65f), dropoutP, rng);
    return;
  }
  if (ch == ':') {
    const float r = std::max(1.2f, thickness * 0.65f);
    segDot(draw, ImVec2(xr, y0 + h * 0.35f), col, r, dropoutP, rng);
    segDot(draw, ImVec2(xr, y0 + h * 0.70f), col, r, dropoutP, rng);
    return;
  }
  if (ch == '-') {
    segLine(draw, ImVec2(xl, ymid), ImVec2(xr, ymid), col, thickness, dropoutP, rng);
    return;
  }
  if (ch == '+') {
    segLine(draw, ImVec2(xl, ymid), ImVec2(xr, ymid), col, thickness, dropoutP, rng);
    segLine(draw, ImVec2(xm, yt), ImVec2(xm, yb), col, thickness, dropoutP, rng);
    return;
  }
  if (ch == '/') {
    segLine(draw, ImVec2(xl, yb), ImVec2(xr, yt), col, thickness, dropoutP, rng);
    return;
  }
  if (ch == '%') {
    // Simple percent: diagonal slash + two dots.
    segLine(draw, ImVec2(xl, yb), ImVec2(xr, yt), col, thickness, dropoutP, rng);
    const float r = std::max(1.2f, thickness * 0.60f);
    segDot(draw, ImVec2(xl + (xr - xl) * 0.25f, yt + (yb - yt) * 0.25f), col, r, dropoutP, rng);
    segDot(draw, ImVec2(xl + (xr - xl) * 0.75f, yt + (yb - yt) * 0.75f), col, r, dropoutP, rng);
    return;
  }

  unsigned int m = segMaskForChar(ch);

  // Decimal point bit can be driven by the mapping (rare).
  if ((m & SDP) != 0u) {
    segDot(draw, ImVec2(xr, yb), col, std::max(1.2f, thickness * 0.65f), dropoutP, rng);
    m &= ~SDP;
  }

  // Optional soft glow (very subtle) for readability on dark backgrounds.
  if ((col >> 24) > 0) {
    const ImU32 glow = scaleAlpha(col, 0.22f);
    const float gth = thickness * 2.4f;
    // Only draw glow when not heavily dropping segments.
    if (dropoutP < 0.15f) {
      if (m & SA)  segLine(draw, ImVec2(xl, yt),   ImVec2(xr, yt),   glow, gth, 0.0f, nullptr);
      if (m & SD)  segLine(draw, ImVec2(xl, yb),   ImVec2(xr, yb),   glow, gth, 0.0f, nullptr);
      if (m & SG1) segLine(draw, ImVec2(xl, ymid), ImVec2(xm, ymid), glow, gth, 0.0f, nullptr);
      if (m & SG2) segLine(draw, ImVec2(xm, ymid), ImVec2(xr, ymid), glow, gth, 0.0f, nullptr);
      if (m & SF)  segLine(draw, ImVec2(xl, yt),   ImVec2(xl, ymid), glow, gth, 0.0f, nullptr);
      if (m & SE)  segLine(draw, ImVec2(xl, ymid), ImVec2(xl, yb),   glow, gth, 0.0f, nullptr);
      if (m & SB)  segLine(draw, ImVec2(xr, yt),   ImVec2(xr, ymid), glow, gth, 0.0f, nullptr);
      if (m & SC)  segLine(draw, ImVec2(xr, ymid), ImVec2(xr, yb),   glow, gth, 0.0f, nullptr);
      if (m & SH)  segLine(draw, ImVec2(xl, yt),   ImVec2(xm, ymid), glow, gth, 0.0f, nullptr);
      if (m & SI)  segLine(draw, ImVec2(xr, yt),   ImVec2(xm, ymid), glow, gth, 0.0f, nullptr);
      if (m & SJ)  segLine(draw, ImVec2(xl, yb),   ImVec2(xm, ymid), glow, gth, 0.0f, nullptr);
      if (m & SK)  segLine(draw, ImVec2(xr, yb),   ImVec2(xm, ymid), glow, gth, 0.0f, nullptr);
      if (m & SL)  segLine(draw, ImVec2(xm, yt),   ImVec2(xm, ymid), glow, gth, 0.0f, nullptr);
      if (m & SM)  segLine(draw, ImVec2(xm, ymid), ImVec2(xm, yb),   glow, gth, 0.0f, nullptr);
    }
  }

  // Main strokes (subject to dropout).
  if (m & SA)  segLine(draw, ImVec2(xl, yt),   ImVec2(xr, yt),   col, thickness, dropoutP, rng);
  if (m & SD)  segLine(draw, ImVec2(xl, yb),   ImVec2(xr, yb),   col, thickness, dropoutP, rng);
  if (m & SG1) segLine(draw, ImVec2(xl, ymid), ImVec2(xm, ymid), col, thickness, dropoutP, rng);
  if (m & SG2) segLine(draw, ImVec2(xm, ymid), ImVec2(xr, ymid), col, thickness, dropoutP, rng);
  if (m & SF)  segLine(draw, ImVec2(xl, yt),   ImVec2(xl, ymid), col, thickness, dropoutP, rng);
  if (m & SE)  segLine(draw, ImVec2(xl, ymid), ImVec2(xl, yb),   col, thickness, dropoutP, rng);
  if (m & SB)  segLine(draw, ImVec2(xr, yt),   ImVec2(xr, ymid), col, thickness, dropoutP, rng);
  if (m & SC)  segLine(draw, ImVec2(xr, ymid), ImVec2(xr, yb),   col, thickness, dropoutP, rng);
  if (m & SH)  segLine(draw, ImVec2(xl, yt),   ImVec2(xm, ymid), col, thickness, dropoutP, rng);
  if (m & SI)  segLine(draw, ImVec2(xr, yt),   ImVec2(xm, ymid), col, thickness, dropoutP, rng);
  if (m & SJ)  segLine(draw, ImVec2(xl, yb),   ImVec2(xm, ymid), col, thickness, dropoutP, rng);
  if (m & SK)  segLine(draw, ImVec2(xr, yb),   ImVec2(xm, ymid), col, thickness, dropoutP, rng);
  if (m & SL)  segLine(draw, ImVec2(xm, yt),   ImVec2(xm, ymid), col, thickness, dropoutP, rng);
  if (m & SM)  segLine(draw, ImVec2(xm, ymid), ImVec2(xm, yb),   col, thickness, dropoutP, rng);
}

static void drawSegText(ImDrawList* draw,
                        ImVec2 pos,
                        float h,
                        const char* s,
                        ImU32 col,
                        float thickness,
                        float spacing,
                        float dropoutP,
                        core::SplitMix64* rng) {
  if (!s || !*s) return;
  const float adv = segCharWidth(h) + spacing;

  ImVec2 p = pos;
  for (const char* it = s; *it; ++it) {
    // Add a tiny deterministic jitter under glitch to read as "signal instability".
    ImVec2 j = ImVec2(0, 0);
    if (dropoutP > 0.02f && rng) {
      const float mag = thickness * 0.65f * dropoutP;
      j.x = (float)rng->range<double>(-mag, mag);
      j.y = (float)rng->range<double>(-mag, mag);
    }

    drawSegChar(draw, ImVec2(p.x + j.x, p.y + j.y), h, *it, col, thickness, dropoutP, rng);
    p.x += adv;
  }
}
static void drawSpark(ImDrawList* draw,
                      const ShipHudSparkline& h,
                      ImVec2 a,
                      ImVec2 b,
                      ImU32 c,
                      float uiScale) {
  if (h.count < 2) return;
  const int n = h.count;
  const float w = b.x - a.x;
  const float hgt = b.y - a.y;
  ImVec2 prev = ImVec2(a.x, b.y - h.get(0) * hgt);
  for (int i = 1; i < n; ++i) {
    const float t = (n <= 1) ? 0.0f : (float)i / (float)(n - 1);
    ImVec2 p = ImVec2(a.x + t * w, b.y - h.get(i) * hgt);
    draw->AddLine(prev, p, c, 1.0f * uiScale);
    prev = p;
  }
}

// --- Panel decor -----------------------------------------------------------

enum class DecorKind : int {
  Circuit = 0,
  Hatch,
  Radar,
  Barcode
};

static DecorKind chooseDecorKind(ShipHudInstrument inst, core::SplitMix64& rng) {
  // Slightly bias based on instrument so HUDs "feel" themed.
  switch (inst) {
    case ShipHudInstrument::Target: return DecorKind::Radar;
    case ShipHudInstrument::Shield: return rng.chance(0.55) ? DecorKind::Circuit : DecorKind::Hatch;
    case ShipHudInstrument::Heat: return rng.chance(0.60) ? DecorKind::Hatch : DecorKind::Barcode;
    case ShipHudInstrument::Fuel: return rng.chance(0.55) ? DecorKind::Barcode : DecorKind::Circuit;
    default: break;
  }
  const int k = rng.range<int>(0, 3);
  return (DecorKind)k;
}

static void drawDecorHatch(ImDrawList* draw, ImVec2 pMin, ImVec2 pMax, core::SplitMix64& rng, ImU32 col, float uiScale) {
  const float margin = 6.0f * uiScale;
  pMin.x += margin;
  pMin.y += margin;
  pMax.x -= margin;
  pMax.y -= margin;

  const float w = std::max(1.0f, pMax.x - pMin.x);
  const float h = std::max(1.0f, pMax.y - pMin.y);
  const float spacing = (8.0f + (float)rng.range<int>(0, 6)) * uiScale;

  const bool flip = rng.chance(0.5);
  const float diag = w + h;
  const int n = (int)std::ceil(diag / spacing) + 2;
  for (int i = -1; i < n; ++i) {
    const float t = (float)i * spacing;
    ImVec2 a = flip ? ImVec2(pMin.x + t, pMin.y) : ImVec2(pMin.x + t, pMax.y);
    ImVec2 b = flip ? ImVec2(a.x - h, pMax.y) : ImVec2(a.x - h, pMin.y);
    draw->AddLine(a, b, col, 1.0f * uiScale);
  }
}

static void drawDecorRadar(ImDrawList* draw, ImVec2 pMin, ImVec2 pMax, core::SplitMix64& rng, ImU32 col, float uiScale) {
  const float margin = 7.0f * uiScale;
  pMin.x += margin;
  pMin.y += margin;
  pMax.x -= margin;
  pMax.y -= margin;

  const ImVec2 c = ImVec2((pMin.x + pMax.x) * 0.5f, (pMin.y + pMax.y) * 0.5f);
  const float r = 0.48f * std::min(pMax.x - pMin.x, pMax.y - pMin.y);

  const int rings = 2 + rng.range<int>(0, 2);
  for (int i = 1; i <= rings; ++i) {
    draw->AddCircle(c, r * (float)i / (float)rings, col, 0, 1.0f * uiScale);
  }
  draw->AddLine(ImVec2(c.x - r, c.y), ImVec2(c.x + r, c.y), col, 1.0f * uiScale);
  draw->AddLine(ImVec2(c.x, c.y - r), ImVec2(c.x, c.y + r), col, 1.0f * uiScale);

  const int blips = 3 + rng.range<int>(0, 5);
  for (int i = 0; i < blips; ++i) {
    const float a = (float)rng.range<double>(0.0, 6.28318530718);
    const float rr = r * (float)rng.range<double>(0.1, 0.95);
    ImVec2 p = ImVec2(c.x + std::cos(a) * rr, c.y + std::sin(a) * rr);
    draw->AddCircleFilled(p, (1.5f + (float)rng.range<int>(0, 2)) * uiScale, col, 8);
  }
}

static void drawDecorBarcode(ImDrawList* draw, ImVec2 pMin, ImVec2 pMax, core::SplitMix64& rng, ImU32 col, float uiScale) {
  const float margin = 6.0f * uiScale;
  pMin.x += margin;
  pMin.y += margin;
  pMax.x -= margin;
  pMax.y -= margin;
  const float w = std::max(1.0f, pMax.x - pMin.x);
  const float h = std::max(1.0f, pMax.y - pMin.y);

  const int bars = 10 + rng.range<int>(0, 10);
  float x = pMin.x;
  for (int i = 0; i < bars && x < pMax.x; ++i) {
    const float bw = (1.0f + (float)rng.range<int>(0, 2)) * uiScale;
    const float bh = h * (float)rng.range<double>(0.25, 0.95);
    const float y0 = pMax.y - bh;
    draw->AddLine(ImVec2(x, y0), ImVec2(x, pMax.y), col, bw);
    x += (bw + (2.0f + (float)rng.range<int>(0, 3)) * uiScale);
  }
}

static void drawDecorCircuit(ImDrawList* draw,
                             ImVec2 pMin,
                             ImVec2 pMax,
                             core::SplitMix64& rng,
                             ImU32 col,
                             float uiScale,
                             int detailLevel) {
  const float margin = 6.0f * uiScale;
  pMin.x += margin;
  pMin.y += margin;
  pMax.x -= margin;
  pMax.y -= margin;
  const float w = std::max(1.0f, pMax.x - pMin.x);
  const float h = std::max(1.0f, pMax.y - pMin.y);

  const int cols = std::clamp((int)std::floor(w / (42.0f * uiScale)), 6, 16);
  const int rows = std::clamp((int)std::floor(h / (36.0f * uiScale)), 4, 12);
  if (cols < 2 || rows < 2) return;

  const float dx = w / (float)(cols - 1);
  const float dy = h / (float)(rows - 1);
  auto cell = [&](int cx, int cy) {
    return ImVec2(pMin.x + (float)cx * dx, pMin.y + (float)cy * dy);
  };

  const int traces = 1 + rng.range<int>(0, 1 + detailLevel);
  for (int t = 0; t < traces; ++t) {
    // Pick a border start.
    int x = 0, y = 0;
    const int side = rng.range<int>(0, 3);
    if (side == 0) { // left
      x = 0;
      y = rng.range<int>(0, rows - 1);
    } else if (side == 1) { // right
      x = cols - 1;
      y = rng.range<int>(0, rows - 1);
    } else if (side == 2) { // top
      x = rng.range<int>(0, cols - 1);
      y = 0;
    } else { // bottom
      x = rng.range<int>(0, cols - 1);
      y = rows - 1;
    }

    int dir = rng.range<int>(0, 3); // 0=E 1=W 2=S 3=N
    const int steps = 6 + rng.range<int>(0, 10 + detailLevel * 4);
    ImVec2 prev = cell(x, y);
    draw->AddCircleFilled(prev, 1.2f * uiScale, col, 8);

    for (int i = 0; i < steps; ++i) {
      // Turn a bit.
      const int turn = rng.range<int>(0, 99);
      if (turn < 18) {
        // left
        dir = (dir == 0) ? 3 : (dir == 3) ? 1 : (dir == 1) ? 2 : 0;
      } else if (turn > 82) {
        // right
        dir = (dir == 0) ? 2 : (dir == 2) ? 1 : (dir == 1) ? 3 : 0;
      }

      int nx = x, ny = y;
      if (dir == 0) nx += 1;
      else if (dir == 1) nx -= 1;
      else if (dir == 2) ny += 1;
      else ny -= 1;

      // Clamp; if we hit the wall, bounce.
      if (nx < 0 || nx >= cols || ny < 0 || ny >= rows) {
        dir ^= 1; // flip E<->W or S<->N in a cheap deterministic way
        nx = std::clamp(nx, 0, cols - 1);
        ny = std::clamp(ny, 0, rows - 1);
      }

      x = nx;
      y = ny;
      ImVec2 cur = cell(x, y);
      draw->AddLine(prev, cur, col, 1.0f * uiScale);
      if (rng.chance(0.18)) {
        draw->AddCircleFilled(cur, 1.0f * uiScale, col, 8);
      }
      prev = cur;
    }
  }
}

static void drawPanelDecor(ImDrawList* draw,
                           ImVec2 pMin,
                           ImVec2 pMax,
                           ShipHudInstrument inst,
                           core::SplitMix64& rng,
                           ImU32 col,
                           float uiScale,
                           int detailLevel) {
  const DecorKind kind = chooseDecorKind(inst, rng);
  switch (kind) {
    case DecorKind::Circuit:
      drawDecorCircuit(draw, pMin, pMax, rng, col, uiScale, detailLevel);
      break;
    case DecorKind::Hatch:
      drawDecorHatch(draw, pMin, pMax, rng, col, uiScale);
      break;
    case DecorKind::Radar:
      drawDecorRadar(draw, pMin, pMax, rng, col, uiScale);
      break;
    case DecorKind::Barcode:
      drawDecorBarcode(draw, pMin, pMax, rng, col, uiScale);
      break;
  }
}

// --- Glyphs ---------------------------------------------------------------

static void drawInstrumentGlyph(ImDrawList* draw,
                                ShipHudInstrument inst,
                                ImVec2 c,
                                float sz,
                                ImU32 col,
                                core::SplitMix64& rng,
                                float thickness) {
  const float r = sz * 0.5f;
  const float t = thickness;

  auto line = [&](ImVec2 a, ImVec2 b) {
    draw->AddLine(a, b, col, t);
  };

  switch (inst) {
    case ShipHudInstrument::Speed: {
      // Arrow + motion lines
      ImVec2 tip = ImVec2(c.x + r * 0.75f, c.y);
      ImVec2 base = ImVec2(c.x - r * 0.55f, c.y);
      line(base, tip);
      line(ImVec2(tip.x - r * 0.28f, tip.y - r * 0.28f), tip);
      line(ImVec2(tip.x - r * 0.28f, tip.y + r * 0.28f), tip);
      const int streaks = 2 + rng.range<int>(0, 2);
      for (int i = 0; i < streaks; ++i) {
        const float yy = lerpf(-0.30f, 0.30f, (float)i / (float)std::max(1, streaks - 1));
        const float len = r * (0.25f + (float)rng.range<double>(0.0, 0.25));
        line(ImVec2(c.x - r * 0.75f, c.y + yy * r), ImVec2(c.x - r * 0.75f + len, c.y + yy * r));
      }
    } break;
    case ShipHudInstrument::Shield: {
      // Hex-ish shield
      const float rr = r * 0.85f;
      ImVec2 p0 = ImVec2(c.x, c.y - rr);
      ImVec2 p1 = ImVec2(c.x + rr, c.y - rr * 0.2f);
      ImVec2 p2 = ImVec2(c.x + rr * 0.65f, c.y + rr);
      ImVec2 p3 = ImVec2(c.x - rr * 0.65f, c.y + rr);
      ImVec2 p4 = ImVec2(c.x - rr, c.y - rr * 0.2f);
      line(p0, p1);
      line(p1, p2);
      line(p2, p3);
      line(p3, p4);
      line(p4, p0);
    } break;
    case ShipHudInstrument::Hull: {
      // Box with diagonal brace
      ImVec2 a = ImVec2(c.x - r * 0.75f, c.y - r * 0.60f);
      ImVec2 b = ImVec2(c.x + r * 0.75f, c.y + r * 0.60f);
      draw->AddRect(a, b, col, 0.0f, 0, t);
      line(ImVec2(a.x, b.y), ImVec2(b.x, a.y));
    } break;
    case ShipHudInstrument::Heat: {
      // Wavy line
      const int segs = 6;
      ImVec2 prev = ImVec2(c.x - r, c.y);
      for (int i = 1; i <= segs; ++i) {
        const float tt = (float)i / (float)segs;
        const float x = lerpf(c.x - r, c.x + r, tt);
        const float y = c.y + std::sin(tt * 6.2831853f) * r * 0.35f;
        ImVec2 cur = ImVec2(x, y);
        line(prev, cur);
        prev = cur;
      }
    } break;
    case ShipHudInstrument::Fuel: {
      // Droplet-ish: circle + point
      draw->AddCircle(c, r * 0.55f, col, 0, t);
      line(ImVec2(c.x, c.y + r * 0.10f), ImVec2(c.x, c.y + r * 0.85f));
      line(ImVec2(c.x - r * 0.25f, c.y + r * 0.55f), ImVec2(c.x, c.y + r * 0.85f));
      line(ImVec2(c.x + r * 0.25f, c.y + r * 0.55f), ImVec2(c.x, c.y + r * 0.85f));
    } break;
    case ShipHudInstrument::Pips: {
      // Three stacked bars
      const float w = r * 1.5f;
      const float h = r * 0.25f;
      for (int i = 0; i < 3; ++i) {
        const float yy = c.y - r * 0.55f + i * (h + r * 0.20f);
        draw->AddRect(ImVec2(c.x - w * 0.5f, yy), ImVec2(c.x + w * 0.5f, yy + h), col, 0.0f, 0, t);
      }
    } break;
    case ShipHudInstrument::Fsd: {
      // Star + outward arrow
      draw->AddCircle(c, r * 0.35f, col, 0, t);
      line(ImVec2(c.x, c.y - r), ImVec2(c.x, c.y - r * 0.35f));
      line(ImVec2(c.x, c.y + r), ImVec2(c.x, c.y + r * 0.35f));
      line(ImVec2(c.x - r, c.y), ImVec2(c.x - r * 0.35f, c.y));
      line(ImVec2(c.x + r, c.y), ImVec2(c.x + r * 0.35f, c.y));
      ImVec2 tip = ImVec2(c.x + r * 0.90f, c.y - r * 0.65f);
      line(ImVec2(c.x + r * 0.25f, c.y - r * 0.20f), tip);
      line(ImVec2(tip.x - r * 0.25f, tip.y), tip);
      line(ImVec2(tip.x, tip.y + r * 0.25f), tip);
    } break;
    case ShipHudInstrument::Throttle: {
      // Slider bar
      ImVec2 a = ImVec2(c.x - r * 0.85f, c.y);
      ImVec2 b = ImVec2(c.x + r * 0.85f, c.y);
      line(a, b);
      const float knob = lerpf(a.x, b.x, (float)rng.range<double>(0.15, 0.85));
      draw->AddCircleFilled(ImVec2(knob, c.y), r * 0.18f, col, 12);
    } break;
    case ShipHudInstrument::Target: {
      // Crosshair
      draw->AddCircle(c, r * 0.65f, col, 0, t);
      line(ImVec2(c.x - r, c.y), ImVec2(c.x - r * 0.35f, c.y));
      line(ImVec2(c.x + r * 0.35f, c.y), ImVec2(c.x + r, c.y));
      line(ImVec2(c.x, c.y - r), ImVec2(c.x, c.y - r * 0.35f));
      line(ImVec2(c.x, c.y + r * 0.35f), ImVec2(c.x, c.y + r));
    } break;
    case ShipHudInstrument::Cargo: {
      // Crate
      ImVec2 a = ImVec2(c.x - r * 0.70f, c.y - r * 0.55f);
      ImVec2 b = ImVec2(c.x + r * 0.70f, c.y + r * 0.55f);
      draw->AddRect(a, b, col, 0.0f, 0, t);
      line(ImVec2(a.x, c.y), ImVec2(b.x, c.y));
      line(ImVec2(c.x, a.y), ImVec2(c.x, b.y));
    } break;
    case ShipHudInstrument::GForce: {
      // Stylized "G"
      draw->AddCircle(c, r * 0.70f, col, 0, t);
      line(ImVec2(c.x, c.y), ImVec2(c.x + r * 0.70f, c.y));
      line(ImVec2(c.x + r * 0.70f, c.y), ImVec2(c.x + r * 0.70f, c.y + r * 0.35f));
    } break;
case ShipHudInstrument::Attitude: {
      // Horizon circle + reference line
      draw->AddCircle(c, r * 0.70f, col, 0, t);
      line(ImVec2(c.x - r * 0.70f, c.y), ImVec2(c.x + r * 0.70f, c.y));
      line(ImVec2(c.x, c.y - r * 0.70f), ImVec2(c.x, c.y - r * 0.25f));
      // Roll ticks
      line(ImVec2(c.x - r * 0.40f, c.y - r * 0.55f), ImVec2(c.x - r * 0.25f, c.y - r * 0.70f));
      line(ImVec2(c.x + r * 0.40f, c.y - r * 0.55f), ImVec2(c.x + r * 0.25f, c.y - r * 0.70f));
    } break;
    default: {
      draw->AddCircle(c, r * 0.60f, col, 0, t);
    } break;
  }
}

static void drawMicrotext(ImDrawList* draw,
                          ImVec2 pMin,
                          ImVec2 pMax,
                          core::SplitMix64& rng,
                          ImU32 col,
                          float uiScale,
                          int lines) {
  static const char* kA[] = {
    "SYS", "CAL", "AX", "SIG", "TLM", "BUS", "LO", "HI", "DMP", "FLT"
  };
  static const char* kB[] = {
    "OK", "NOM", "Δ", "PH", "MOD", "LOCK", "SYNC", "BYP", "CRC", "DRIFT"
  };

  const float lineH = 11.0f * uiScale;
  ImVec2 pos = ImVec2(pMin.x + 6.0f * uiScale, pMax.y - (float)lines * lineH - 20.0f * uiScale);

  char buf[48];
  for (int i = 0; i < lines; ++i) {
    const int a = rng.range<int>(0, (int)(sizeof(kA) / sizeof(kA[0])) - 1);
    const int b = rng.range<int>(0, (int)(sizeof(kB) / sizeof(kB[0])) - 1);
    const int n = rng.range<int>(0, 0xFFFF);
    std::snprintf(buf, sizeof(buf), "%s %s %04X", kA[a], kB[b], n);
    draw->AddText(ImVec2(pos.x, pos.y + (float)i * lineH), col, buf);
  }
}

// --- Gauge helpers --------------------------------------------------------

static void drawTicks(ImDrawList* draw,
                      ImVec2 c,
                      float r0,
                      float r1,
                      float a0,
                      float a1,
                      int major,
                      int minor,
                      ImU32 col,
                      float uiScale) {
  major = std::max(1, major);
  minor = std::max(major, minor);

  // Minor ticks
  for (int i = 0; i <= minor; ++i) {
    const float t = (minor <= 0) ? 0.0f : (float)i / (float)minor;
    const float a = lerpf(a0, a1, t);
    const bool isMajor = (i % std::max(1, minor / major)) == 0;
    const float rr0 = isMajor ? r0 : (r0 + (r1 - r0) * 0.55f);
    const float rr1 = r1;
    ImVec2 p0 = ImVec2(c.x + std::cos(a) * rr0, c.y + std::sin(a) * rr0);
    ImVec2 p1 = ImVec2(c.x + std::cos(a) * rr1, c.y + std::sin(a) * rr1);
    draw->AddLine(p0, p1, col, (isMajor ? 1.4f : 1.0f) * uiScale);
  }
}

static void drawBar(ImDrawList* draw,
                    ImVec2 pMin,
                    ImVec2 pMax,
                    float frac,
                    ImU32 colGrid,
                    ImU32 fillCol,
                    float uiScale,
                    int majorTicks) {
  frac = clamp01(frac);
  const float margin = 6.0f * uiScale;
  ImVec2 bMin = ImVec2(pMin.x + margin, pMax.y - margin - 10.0f * uiScale);
  ImVec2 bMax = ImVec2(pMax.x - margin, pMax.y - margin);
  draw->AddRect(bMin, bMax, colGrid, 0.0f, 0, 1.0f * uiScale);
  ImVec2 fMax = ImVec2(bMin.x + (bMax.x - bMin.x) * frac, bMax.y);
  draw->AddRectFilled(bMin, fMax, fillCol, 0.0f);

  // Tick marks
  majorTicks = std::clamp(majorTicks, 3, 10);
  for (int i = 0; i <= majorTicks; ++i) {
    const float t = (float)i / (float)majorTicks;
    const float x = lerpf(bMin.x, bMax.x, t);
    const float len = (i == 0 || i == majorTicks) ? 6.0f : 4.0f;
    draw->AddLine(ImVec2(x, bMin.y - len * uiScale), ImVec2(x, bMin.y), colGrid, 1.0f * uiScale);
  }
}

static void drawArc(ImDrawList* draw,
                    ImVec2 pMin,
                    ImVec2 pMax,
                    float frac,
                    ImU32 colGrid,
                    ImU32 arcCol,
                    float uiScale,
                    int majorTicks,
                    int minorTicks) {
  frac = clamp01(frac);
  ImVec2 c = ImVec2((pMin.x + pMax.x) * 0.5f, (pMin.y + pMax.y) * 0.55f);
  const float r = std::min(pMax.x - pMin.x, pMax.y - pMin.y) * 0.38f;
  const float a0 = 3.14159f * 0.75f;
  const float a1 = 3.14159f * 2.25f;

  drawTicks(draw, c, r * 0.88f, r * 1.02f, a0, a1, majorTicks, minorTicks, colGrid, uiScale);

  draw->PathArcTo(c, r, a0, a1, 28);
  draw->PathStroke(colGrid, false, 2.0f * uiScale);
  draw->PathArcTo(c, r, a0, a0 + (a1 - a0) * frac, 28);
  draw->PathStroke(arcCol, false, 3.0f * uiScale);
}

static void drawDial(ImDrawList* draw,
                     ImVec2 pMin,
                     ImVec2 pMax,
                     float frac,
                     ImU32 colGrid,
                     ImU32 needleCol,
                     float uiScale,
                     int majorTicks,
                     int minorTicks) {
  frac = clamp01(frac);
  ImVec2 c = ImVec2((pMin.x + pMax.x) * 0.5f, (pMin.y + pMax.y) * 0.58f);
  const float r = std::min(pMax.x - pMin.x, pMax.y - pMin.y) * 0.34f;
  draw->AddCircle(c, r, colGrid, 28, 2.0f * uiScale);

  const float a0 = 3.14159f * 0.75f;
  const float a1 = 3.14159f * 2.25f;
  drawTicks(draw, c, r * 0.80f, r * 0.98f, a0, a1, majorTicks, minorTicks, colGrid, uiScale);

  const float a = a0 + (a1 - a0) * frac;
  ImVec2 tip = ImVec2(c.x + std::cos(a) * r, c.y + std::sin(a) * r);
  draw->AddLine(c, tip, needleCol, 2.0f * uiScale);
  draw->AddCircleFilled(c, 2.5f * uiScale, needleCol, 12);
}

static void drawSegmented(ImDrawList* draw,
                          ImVec2 pMin,
                          ImVec2 pMax,
                          float frac,
                          ImU32 colGrid,
                          ImU32 segCol,
                          float uiScale,
                          int segs) {
  frac = clamp01(frac);
  segs = std::clamp(segs, 6, 16);
  const float gap = 2.0f * uiScale;
  const float margin = 6.0f * uiScale;
  const float w = (pMax.x - pMin.x) - margin * 2.0f;
  const float h = 10.0f * uiScale;
  const float segW = (w - gap * (segs - 1)) / (float)segs;
  ImVec2 base = ImVec2(pMin.x + margin, pMax.y - margin - h);
  const int filled = (int)std::floor(frac * segs + 0.0001f);
  for (int i = 0; i < segs; ++i) {
    ImVec2 s0 = ImVec2(base.x + i * (segW + gap), base.y);
    ImVec2 s1 = ImVec2(s0.x + segW, base.y + h);
    draw->AddRect(s0, s1, colGrid, 0.0f, 0, 1.0f * uiScale);
    if (i < filled) draw->AddRectFilled(s0, s1, segCol, 0.0f);
  }
}

static float degToRadF(float deg) {
  return deg * 0.017453292519943295769f;
}

static bool lineCircleIntersection(ImVec2 c, float r, ImVec2 p0, ImVec2 dir, ImVec2& outA, ImVec2& outB) {
  // dir must be normalized.
  const float ox = p0.x - c.x;
  const float oy = p0.y - c.y;

  const float b = ox * dir.x + oy * dir.y;
  const float cterm = ox * ox + oy * oy - r * r;
  const float disc = b * b - cterm;
  if (disc < 0.0f) return false;

  const float s = std::sqrt(disc);
  const float t0 = -b - s;
  const float t1 = -b + s;

  outA = ImVec2(p0.x + dir.x * t0, p0.y + dir.y * t0);
  outB = ImVec2(p0.x + dir.x * t1, p0.y + dir.y * t1);
  return true;
}

static ImVec2 clampPointToCircle(ImVec2 c, float r, ImVec2 p) {
  const float dx = p.x - c.x;
  const float dy = p.y - c.y;
  const float d2 = dx * dx + dy * dy;
  const float rr = r * r;
  if (d2 <= rr) return p;

  const float d = std::sqrt(std::max(1e-12f, d2));
  const float s = r / d;
  return ImVec2(c.x + dx * s, c.y + dy * s);
}

static void drawAttitudeIndicator(ImDrawList* draw,
                                 ImVec2 pMin,
                                 ImVec2 pMax,
                                 const ShipHudOverlayTelemetry& t,
                                 ImU32 colGrid,
                                 ImU32 colText,
                                 ImU32 colPrim,
                                 ImU32 colAcc,
                                 ImU32 colDanger,
                                 float uiScale,
                                 int theme) {
  if (!draw) return;

  const float margin = 8.0f * uiScale;
  const float topPad = 22.0f * uiScale;
  const float botPad = 14.0f * uiScale;

  ImVec2 a = ImVec2(pMin.x + margin, pMin.y + topPad);
  ImVec2 b = ImVec2(pMax.x - margin, pMax.y - botPad);

  const float w = std::max(1.0f, b.x - a.x);
  const float h = std::max(1.0f, b.y - a.y);

  ImVec2 c = ImVec2((a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f);
  const float r = std::min(w, h) * 0.48f;

  const float ringThick = (theme == 2 ? 2.2f : 1.8f) * uiScale;
  const float lineThick = (theme == 3 ? 1.10f : 1.30f) * uiScale;

  // Outer ring + center dot
  draw->AddCircle(c, r, colGrid, 48, ringThick);
  draw->AddCircleFilled(c, std::max(1.0f, r * 0.018f), colGrid, 12);

  if (!t.attitudeValid) {
    const char* msg = "NO REF";
    ImVec2 sz = ImGui::CalcTextSize(msg);
    draw->AddText(ImVec2(c.x - sz.x * 0.5f, c.y - sz.y * 0.5f), colDanger, msg);
    return;
  }

  const float maxPitch = 45.0f;
  const float pitch = std::clamp(t.attitudePitchDeg, -89.0f, 89.0f);
  const float roll = t.attitudeRollDeg;

  const float angle = -degToRadF(roll);
  ImVec2 dir = ImVec2(std::cos(angle), std::sin(angle));
  ImVec2 nrm = ImVec2(-dir.y, dir.x); // perpendicular (down for angle=0)

  const float pitchNorm = std::clamp(pitch / maxPitch, -1.2f, 1.2f);
  const float offset = pitchNorm * (r * 0.65f);
  ImVec2 hPt = ImVec2(c.x + nrm.x * offset, c.y + nrm.y * offset);

  // Horizon line
  {
    ImVec2 p0, p1;
    if (lineCircleIntersection(c, r * 0.985f, hPt, dir, p0, p1)) {
      draw->AddLine(p0, p1, colGrid, lineThick);
    }
  }

  // Pitch ladder (±10/20/30)
  const int steps[] = {10, 20, 30};
  for (int si = 0; si < 3; ++si) {
    const float deg = (float)steps[si];
    for (int sgn = -1; sgn <= 1; sgn += 2) {
      const float pDeg = deg * (float)sgn;
      const float o = std::clamp(pDeg / maxPitch, -2.0f, 2.0f) * (r * 0.65f);
      ImVec2 lPt = ImVec2(hPt.x - nrm.x * o, hPt.y - nrm.y * o);

      ImVec2 a0, b0;
      if (!lineCircleIntersection(c, r * 0.985f, lPt, dir, a0, b0)) continue;

      // Shorten the ladder line.
      const float shrink = (deg == 10.0f) ? 0.60f : (deg == 20.0f) ? 0.52f : 0.44f;
      ImVec2 mid = ImVec2((a0.x + b0.x) * 0.5f, (a0.y + b0.y) * 0.5f);
      ImVec2 va = ImVec2(a0.x - mid.x, a0.y - mid.y);
      ImVec2 vb = ImVec2(b0.x - mid.x, b0.y - mid.y);
      ImVec2 a1 = ImVec2(mid.x + va.x * shrink, mid.y + va.y * shrink);
      ImVec2 b1 = ImVec2(mid.x + vb.x * shrink, mid.y + vb.y * shrink);

      draw->AddLine(a1, b1, colGrid, lineThick);

      // End ticks.
      const float tick = 4.0f * uiScale;
      draw->AddLine(a1, ImVec2(a1.x + nrm.x * tick, a1.y + nrm.y * tick), colGrid, lineThick);
      draw->AddLine(b1, ImVec2(b1.x + nrm.x * tick, b1.y + nrm.y * tick), colGrid, lineThick);

      // Label on right side only for readability.
      if (sgn > 0) {
        char buf[8];
        std::snprintf(buf, sizeof(buf), "%d", (int)deg);
        ImVec2 ts = ImGui::CalcTextSize(buf);
        draw->AddText(ImVec2(b1.x + 6.0f * uiScale, b1.y - ts.y * 0.5f), colText, buf);
      }
    }
  }

  // Fixed aircraft symbol (wings + nose mark)
  {
    const float wing = r * 0.42f;
    const float gap = r * 0.10f;
    draw->AddLine(ImVec2(c.x - wing, c.y), ImVec2(c.x - gap, c.y), colText, lineThick);
    draw->AddLine(ImVec2(c.x + gap, c.y), ImVec2(c.x + wing, c.y), colText, lineThick);
    draw->AddLine(ImVec2(c.x, c.y), ImVec2(c.x, c.y + r * 0.16f), colText, lineThick);
  }

  // Bank index notch at top (purely for orientation)
  {
    const float notchR = r * 0.92f;
    const float notchW = 7.0f * uiScale;
    const float notchH = 5.0f * uiScale;
    ImVec2 p0 = ImVec2(c.x, c.y - notchR);
    ImVec2 p1 = ImVec2(c.x - notchW, c.y - notchR + notchH);
    ImVec2 p2 = ImVec2(c.x + notchW, c.y - notchR + notchH);
    draw->AddTriangleFilled(p0, p1, p2, colGrid);
  }

  // Prograde marker (velocity vector)
  if (t.progradeValid) {
    const float maxA = 60.0f;
    const float ox = std::clamp(t.progradeYawDeg / maxA, -1.0f, 1.0f) * (r * 0.70f);
    const float oy = std::clamp(-t.progradePitchDeg / maxA, -1.0f, 1.0f) * (r * 0.70f);
    ImVec2 m = clampPointToCircle(c, r * 0.80f, ImVec2(c.x + ox, c.y + oy));

    const float rr = 6.0f * uiScale;
    draw->AddCircle(m, rr, colAcc, 18, lineThick);
    draw->AddLine(ImVec2(m.x - rr, m.y), ImVec2(m.x + rr, m.y), colAcc, lineThick);
    draw->AddLine(ImVec2(m.x, m.y - rr), ImVec2(m.x, m.y + rr), colAcc, lineThick);
  }


  // Retrograde marker (opposite velocity vector)
  if (t.retrogradeValid) {
    const float maxA = 60.0f;
    const float ox = std::clamp(t.retrogradeYawDeg / maxA, -1.0f, 1.0f) * (r * 0.70f);
    const float oy = std::clamp(-t.retrogradePitchDeg / maxA, -1.0f, 1.0f) * (r * 0.70f);
    ImVec2 m = clampPointToCircle(c, r * 0.80f, ImVec2(c.x + ox, c.y + oy));

    const float rr = 6.0f * uiScale;
    draw->AddCircle(m, rr, colText, 18, lineThick);
    draw->AddLine(ImVec2(m.x - rr, m.y - rr), ImVec2(m.x + rr, m.y + rr), colText, lineThick);
    draw->AddLine(ImVec2(m.x - rr, m.y + rr), ImVec2(m.x + rr, m.y - rr), colText, lineThick);
  }


  // Gravity marker (down vector)
  if (t.gravityValid) {
    const float maxA = 60.0f;
    const float ox = std::clamp(t.gravityYawDeg / maxA, -1.0f, 1.0f) * (r * 0.74f);
    const float oy = std::clamp(-t.gravityPitchDeg / maxA, -1.0f, 1.0f) * (r * 0.74f);
    ImVec2 m = clampPointToCircle(c, r * 0.83f, ImVec2(c.x + ox, c.y + oy));

    const float rr = 6.5f * uiScale;
    ImVec2 a0 = ImVec2(m.x, m.y + rr);
    ImVec2 b0 = ImVec2(m.x - rr, m.y - rr);
    ImVec2 c0 = ImVec2(m.x + rr, m.y - rr);
    draw->AddTriangleFilled(a0, b0, c0, colPrim);
    draw->AddTriangle(a0, b0, c0, colGrid, lineThick);
  }


  // Target direction marker (nearest contact)
  if (t.targetDirValid) {
    const float maxA = 60.0f;
    const float ox = std::clamp(t.targetYawDeg / maxA, -1.0f, 1.0f) * (r * 0.72f);
    const float oy = std::clamp(-t.targetPitchDeg / maxA, -1.0f, 1.0f) * (r * 0.72f);
    ImVec2 m = clampPointToCircle(c, r * 0.82f, ImVec2(c.x + ox, c.y + oy));

    const float rr = 6.0f * uiScale;
    ImVec2 p0(m.x, m.y - rr);
    ImVec2 p1(m.x + rr, m.y);
    ImVec2 p2(m.x, m.y + rr);
    ImVec2 p3(m.x - rr, m.y);
    draw->AddLine(p0, p1, colPrim, lineThick);
    draw->AddLine(p1, p2, colPrim, lineThick);
    draw->AddLine(p2, p3, colPrim, lineThick);
    draw->AddLine(p3, p0, colPrim, lineThick);
  }


  // Reference label strip
  {
    char buf[40];
    if (t.attitudeUsesGravity) {
      std::snprintf(buf, sizeof(buf), "GRAV %.2fg", (double)t.attitudeRefG);
    } else {
      std::snprintf(buf, sizeof(buf), "ORBIT");
    }
    draw->AddText(ImVec2(pMin.x + 6.0f * uiScale, pMax.y - 14.0f * uiScale), colGrid, buf);

    char prBuf[24];
    std::snprintf(prBuf, sizeof(prBuf), "P%+0.0f R%+0.0f", (double)pitch, (double)roll);
    ImVec2 ts = ImGui::CalcTextSize(prBuf);
    draw->AddText(ImVec2(pMax.x - 6.0f * uiScale - ts.x, pMax.y - 14.0f * uiScale), colText, prBuf);
  }
}


} // namespace

void ShipHudSparkline::push(float v) {
  v = clamp01(v);
  samples[(std::size_t)head] = v;
  head = (head + 1) % (int)samples.size();
  count = std::min(count + 1, (int)samples.size());
}

float ShipHudSparkline::get(int i) const {
  if (count <= 0) return 0.0f;
  i = std::clamp(i, 0, count - 1);
  int idx = head - count + i;
  while (idx < 0) idx += (int)samples.size();
  return samples[(std::size_t)(idx % (int)samples.size())];
}

void shipHudHistoryUpdate(ShipHudHistory& h,
                          double dtRealSec,
                          double speedKmS,
                          float shieldFrac,
                          float hullFrac,
                          float heatFrac,
                          float fuelFrac) {
  h.acc += dtRealSec;
  const double step = 0.10;
  while (h.acc >= step) {
    h.acc -= step;

    const float speedNorm = (float)std::clamp(speedKmS / 2.0, 0.0, 1.0);
    h.speed.push(speedNorm);
    h.shield.push(clamp01(shieldFrac));
    h.hull.push(clamp01(hullFrac));
    h.heat.push(clamp01(heatFrac));
    h.fuel.push(clamp01(fuelFrac));
  }
}

void drawShipHudOverlay(ImDrawList* draw,
                        ImVec2 winPos,
                        ImVec2 winSize,
                        float uiScale,
                        const ProceduralShipHudPlan& plan,
                        const ShipHudHistory& hist,
                        const ShipHudOverlayTelemetry& t,
                        const ShipHudOverlayConfig& cfg,
                        const ShipHudOverlayPalette& pal) {
  if (!draw) return;

  const ImU32 colGrid = toU32(pal.grid, 1.0f);
  const ImU32 colText = toU32(pal.text, 1.0f);
  const ImU32 colPrim = toU32(pal.primary, 1.0f);
  const ImU32 colAcc = toU32(pal.accent, 1.0f);
  const ImU32 colDanger = toU32(pal.danger, 1.0f);

  // Theme: 0=Blueprint, 1=Industrial, 2=Military, 3=Minimal
  const int theme = (int)plan.themeVariant % 4;
  const float frameThick = (theme == 2 ? 1.6f : theme == 3 ? 1.10f : 1.25f) * uiScale;
  const float outerRound = (theme == 3 ? 4.0f : 6.0f) * uiScale;
  const float panelRound = (theme == 2 ? 2.0f : 4.0f) * uiScale;
  const int majorTicks = (theme == 3 ? 4 : theme == 0 ? 7 : 6);
  const int minorTicks = (theme == 3 ? 12 : theme == 0 ? 28 : 22);
  const int segs = (theme == 2 ? 12 : theme == 0 ? 14 : 10);
  const float decorAlpha = cfg.showDecor ? std::clamp(cfg.decorAlpha, 0.0f, 1.0f) : 0.0f;
  const ImU32 colDecor = toU32(pal.grid, 0.25f * decorAlpha);
  const ImU32 colMicro = toU32(pal.grid, (theme == 0 ? 0.45f : 0.35f));

  // Outer frame
  draw->AddRect(winPos, ImVec2(winPos.x + winSize.x, winPos.y + winSize.y), colGrid, outerRound, 0, frameThick);

  // Inner client area (leave room for a small title strip)
  const float pad = (theme == 2 ? 11.0f : 10.0f) * uiScale;
  const float titleH = (theme == 3 ? 16.0f : 18.0f) * uiScale;
  ImVec2 innerMin = ImVec2(winPos.x + pad, winPos.y + pad + titleH);
  ImVec2 innerMax = ImVec2(winPos.x + winSize.x - pad, winPos.y + winSize.y - pad);
  const float innerW = std::max(1.0f, innerMax.x - innerMin.x);
  const float innerH = std::max(1.0f, innerMax.y - innerMin.y);

  // Title strip
  {
    const char* themeName = (theme == 0) ? "BP" : (theme == 1) ? "IND" : (theme == 2) ? "MIL" : "MIN";
    char title[32];
    std::snprintf(title, sizeof(title), "SHIP %s", themeName);
    draw->AddText(ImVec2(winPos.x + pad, winPos.y + pad), colText, title);
    char seedBuf[40];
    std::snprintf(seedBuf, sizeof(seedBuf), "#%04llX", (unsigned long long)(plan.seed & 0xFFFFu));
    ImVec2 sSize = ImGui::CalcTextSize(seedBuf);
    draw->AddText(ImVec2(winPos.x + winSize.x - pad - sSize.x, winPos.y + pad), colGrid, seedBuf);
  }

  const float speedKmS = t.speedKmS;
  const float shieldFrac = clamp01(t.shieldFrac);
  const float hullFrac = clamp01(t.hullFrac);
  const float heatFrac = clamp01(t.heatFrac);
  const float fuelFrac = clamp01(t.fuelFrac);

  float glitchStrength = 0.0f;
  if (cfg.glitchFx) {
    glitchStrength += std::max(0.0f, (heatFrac - 0.85f) / 0.15f);
    glitchStrength += std::max(0.0f, (0.25f - shieldFrac) / 0.25f) * 0.35f;
    glitchStrength += std::max(0.0f, (0.20f - hullFrac) / 0.20f) * 0.55f;
    glitchStrength = std::clamp(glitchStrength, 0.0f, 1.0f);
  }

  // Render panels from the plan
  for (const ShipHudPanel& p : plan.panels) {
    ImVec2 pMin = ImVec2(innerMin.x + p.x0 * innerW, innerMin.y + p.y0 * innerH);
    ImVec2 pMax = ImVec2(innerMin.x + p.x1 * innerW, innerMin.y + p.y1 * innerH);

    const float inset = plan.cellPaddingNorm * std::min(innerW, innerH);
    pMin.x += inset;
    pMin.y += inset;
    pMax.x -= inset;
    pMax.y -= inset;

    const float pw = std::max(1.0f, pMax.x - pMin.x);
    const float ph = std::max(1.0f, pMax.y - pMin.y);

    // Panel-local RNG for decor/glyphs/microtext + glitch dropout.
    core::u64 panelSeed = plan.seed;
    panelSeed = core::hashCombine(panelSeed, (core::u64)p.instrument);
    panelSeed = core::hashCombine(panelSeed, (core::u64)p.style);
    panelSeed = core::hashCombine(panelSeed, (core::u64)p.accentVariant);
    core::SplitMix64 prng(core::hashCombine(panelSeed, core::fnv1a64("ship_hud_panel")));

    // Panel dropout / jitter: quantize time so flicker isn't per-frame noisy.
    bool dropout = false;
    ImVec2 jitter(0.0f, 0.0f);
    if (cfg.panelDropouts && glitchStrength > 0.02f) {
      const core::u64 tSlice = (core::u64)std::floor(std::max(0.0, t.timeRealSec) * 8.0);
      core::SplitMix64 frng(core::hashCombine(panelSeed, tSlice));
      float instBias = 1.0f;
      if (p.instrument == ShipHudInstrument::Target) instBias = 1.25f;
      if (p.instrument == ShipHudInstrument::Fsd) instBias = 1.10f;
      if (p.instrument == ShipHudInstrument::Pips) instBias = 0.85f;

      const float chance = std::clamp(glitchStrength * glitchStrength * 0.28f * instBias, 0.0f, 0.65f);
      dropout = frng.chance(chance);

      const float j = (dropout ? 3.5f : 1.5f) * glitchStrength * uiScale;
      jitter.x = (float)frng.range<double>(-j, j);
      jitter.y = (float)frng.range<double>(-j, j);
    }

    pMin.x += jitter.x;
    pMin.y += jitter.y;
    pMax.x += jitter.x;
    pMax.y += jitter.y;

    // Panel frame
    draw->AddRect(pMin, pMax, colGrid, panelRound, 0, 1.0f * uiScale);

    // Background decor
    if (!dropout && decorAlpha > 0.001f) {
      // For minimal theme, make decor sparser.
      const int decorDetail = (theme == 3) ? std::max(0, plan.detailLevel - 1) : plan.detailLevel;
      core::SplitMix64 drng(prng.nextU64());
      drawPanelDecor(draw, pMin, pMax, p.instrument, drng, colDecor, uiScale, decorDetail);
    }

    // Glyph (top-right)
    if (!dropout && cfg.showGlyphs) {
      const float gSz = std::min(pw, ph) * 0.18f;
      const ImVec2 gc = ImVec2(pMax.x - 14.0f * uiScale - gSz * 0.5f, pMin.y + 14.0f * uiScale + gSz * 0.5f);
      core::SplitMix64 grng(prng.nextU64());
      drawInstrumentGlyph(draw, p.instrument, gc, gSz, colGrid, grng, (theme == 2 ? 1.4f : 1.1f) * uiScale);
    }

    // Instrument data
    const char* label = toString(p.instrument);
    float frac = 0.0f;
    ImU32 gaugeCol = colPrim;
    char valBuf[64]{};
    bool hasVal = false;

    switch (p.instrument) {
      case ShipHudInstrument::Speed: {
        if (speedKmS < 1.0) {
          std::snprintf(valBuf, sizeof(valBuf), "%.0f m/s", speedKmS * 1000.0);
        } else {
          std::snprintf(valBuf, sizeof(valBuf), "%.2f km/s", speedKmS);
        }
        hasVal = true;
        frac = (float)std::clamp(speedKmS / 2.0, 0.0, 1.0);
        gaugeCol = colAcc;
      } break;
      case ShipHudInstrument::Shield: {
        std::snprintf(valBuf, sizeof(valBuf), "%.0f%%", shieldFrac * 100.0f);
        hasVal = true;
        frac = shieldFrac;
        if (shieldFrac < 0.25f) gaugeCol = colDanger;
      } break;
      case ShipHudInstrument::Hull: {
        std::snprintf(valBuf, sizeof(valBuf), "%.0f%%", hullFrac * 100.0f);
        hasVal = true;
        frac = hullFrac;
        if (hullFrac < 0.35f) gaugeCol = colDanger;
      } break;
      case ShipHudInstrument::Heat: {
        std::snprintf(valBuf, sizeof(valBuf), "%.0f%%", heatFrac * 100.0f);
        hasVal = true;
        frac = heatFrac;
        if (heatFrac > 0.85f) gaugeCol = colDanger;
      } break;
      case ShipHudInstrument::Fuel: {
        std::snprintf(valBuf, sizeof(valBuf), "%.0f%%", fuelFrac * 100.0f);
        hasVal = true;
        frac = fuelFrac;
      } break;
      case ShipHudInstrument::Pips: {
        std::snprintf(valBuf, sizeof(valBuf), "S%d E%d W%d", t.pipSys, t.pipEng, t.pipWep);
        hasVal = true;
        frac = (float)(t.pipSys + t.pipEng + t.pipWep) / 12.0f;
      } break;
      case ShipHudInstrument::Fsd: {
        std::snprintf(valBuf, sizeof(valBuf), "%s", t.fsdText[0] ? t.fsdText : "--");
        hasVal = true;
        frac = clamp01(t.fsdFrac);
        if (std::strstr(valBuf, "JMP")) gaugeCol = colAcc;
      } break;
      case ShipHudInstrument::Throttle: {
        const float th = clamp01(t.throttleFrac);
        frac = th;
        std::snprintf(valBuf, sizeof(valBuf), "%.0f%%", th * 100.0f);
        hasVal = true;
      } break;
      case ShipHudInstrument::Target: {
        if (t.hasNearest) {
          if (t.nearestKm < 1.0) {
            std::snprintf(valBuf, sizeof(valBuf), "%.0f m", t.nearestKm * 1000.0);
          } else {
            std::snprintf(valBuf, sizeof(valBuf), "%.1f km", t.nearestKm);
          }
          hasVal = true;
          frac = (float)std::clamp(1.0 - (t.nearestKm / 50.0), 0.0, 1.0);
        } else {
          std::snprintf(valBuf, sizeof(valBuf), "--");
          hasVal = true;
          frac = 0.0f;
        }
      } break;
      case ShipHudInstrument::Cargo: {
        std::snprintf(valBuf, sizeof(valBuf), "%.0f cr", t.cargoCr);
        hasVal = true;
        frac = (float)std::clamp(t.cargoCr / 25000.0, 0.0, 1.0);
      } break;
      case ShipHudInstrument::GForce: {
        std::snprintf(valBuf, sizeof(valBuf), "%.2fg", t.gForce);
        hasVal = true;
        frac = (float)std::clamp(t.gForce / 6.0, 0.0, 1.0);
        if (t.gForce > 4.0) gaugeCol = colDanger;
      } break;
      case ShipHudInstrument::Attitude: {
        // Custom instrument: draws horizon + velocity/gravity vectors.
        hasVal = false;
        frac = 0.0f;
        gaugeCol = colGrid;
      } break;
      default: break;
    }

    // Label (top-left)
    draw->AddText(ImVec2(pMin.x + 6.0f * uiScale, pMin.y + 4.0f * uiScale), dropout ? colDecor : colText, label);

    if (dropout) {
      // No signal: strike + noise
      core::SplitMix64 nrng(core::hashCombine(panelSeed, (core::u64)std::floor(t.timeRealSec * 24.0)));
      const int n = 8 + nrng.range<int>(0, 10);
      for (int i = 0; i < n; ++i) {
        const float x0 = (float)nrng.range<double>(pMin.x, pMax.x);
        const float x1 = (float)nrng.range<double>(pMin.x, pMax.x);
        const float y = (float)nrng.range<double>(pMin.y, pMax.y);
        draw->AddLine(ImVec2(x0, y), ImVec2(x1, y), colDecor, 1.0f * uiScale);
      }
      draw->AddLine(pMin, pMax, colDanger, 1.3f * uiScale);
      draw->AddLine(ImVec2(pMin.x, pMax.y), ImVec2(pMax.x, pMin.y), colDanger, 1.3f * uiScale);
      const char* msg = "NO SIGNAL";
      ImVec2 ms = ImGui::CalcTextSize(msg);
      draw->AddText(ImVec2(pMin.x + (pw - ms.x) * 0.5f, pMin.y + (ph - ms.y) * 0.5f), colDanger, msg);
      continue;
    }

    // Sparkline
    if (p.sparkline && cfg.showSparklines) {
      const float spH = 18.0f * uiScale;
      ImVec2 spMin = ImVec2(pMin.x + 6.0f * uiScale, pMin.y + 18.0f * uiScale);
      ImVec2 spMax = ImVec2(pMax.x - 6.0f * uiScale, spMin.y + spH);
      const ShipHudSparkline* sh = nullptr;
      switch (p.instrument) {
        case ShipHudInstrument::Speed: sh = &hist.speed; break;
        case ShipHudInstrument::Shield: sh = &hist.shield; break;
        case ShipHudInstrument::Hull: sh = &hist.hull; break;
        case ShipHudInstrument::Heat: sh = &hist.heat; break;
        case ShipHudInstrument::Fuel: sh = &hist.fuel; break;
        default: sh = nullptr; break;
      }
      if (sh) drawSpark(draw, *sh, spMin, spMax, colDecor, uiScale);
    }

    // Microtext
    if (cfg.showMicrotext && plan.detailLevel >= 2 && p.instrument != ShipHudInstrument::Attitude) {
      core::SplitMix64 mrng(prng.nextU64());
      const int lines = (plan.detailLevel >= 3 && theme != 3) ? 3 : 2;
      drawMicrotext(draw, pMin, pMax, mrng, colMicro, uiScale, lines);
    }

    // Special-case: Attitude draws a navball-style horizon + vector instrument.
    if (p.instrument == ShipHudInstrument::Attitude) {
      drawAttitudeIndicator(draw, pMin, pMax, t, colGrid, colText, colPrim, colAcc, colDanger, uiScale, theme);
    } else if (p.instrument == ShipHudInstrument::Pips) {
      const float margin = 6.0f * uiScale;
      const float barH = 7.0f * uiScale;
      const float gap = 3.0f * uiScale;
      ImVec2 base = ImVec2(pMin.x + margin, pMax.y - margin - (barH * 3.0f + gap * 2.0f));
      auto mini = [&](int idx, int val, const char* lbl) {
        const float y0 = base.y + idx * (barH + gap);
        ImVec2 b0 = ImVec2(base.x, y0);
        ImVec2 b1 = ImVec2(pMax.x - margin, y0 + barH);
        draw->AddRect(b0, b1, colGrid, 0.0f, 0, 1.0f * uiScale);
        const float f = (float)std::clamp((double)val / 4.0, 0.0, 1.0);
        ImVec2 f1 = ImVec2(b0.x + (b1.x - b0.x) * f, b1.y);
        draw->AddRectFilled(b0, f1, gaugeCol, 0.0f);
        draw->AddText(ImVec2(b0.x, y0 - 12.0f * uiScale), colDecor, lbl);
      };
      mini(0, t.pipSys, "SYS");
      mini(1, t.pipEng, "ENG");
      mini(2, t.pipWep, "WEP");
    } else {
      switch (p.style) {
        case ShipHudGaugeStyle::Bar:
          drawBar(draw, pMin, pMax, frac, colGrid, gaugeCol, uiScale, majorTicks);
          break;
        case ShipHudGaugeStyle::Arc:
          drawArc(draw, pMin, pMax, frac, colGrid, gaugeCol, uiScale, majorTicks, minorTicks);
          break;
        case ShipHudGaugeStyle::Dial:
          drawDial(draw, pMin, pMax, frac, colGrid, gaugeCol, uiScale, majorTicks, minorTicks);
          break;
        case ShipHudGaugeStyle::Segmented:
          drawSegmented(draw, pMin, pMax, frac, colGrid, gaugeCol, uiScale, segs);
          break;
        default:
          drawBar(draw, pMin, pMax, frac, colGrid, gaugeCol, uiScale, majorTicks);
          break;
      }
    }

    // Value centered
    if (hasVal) {
      if (cfg.segmentText) {
        // Procedural segment-display readout (optional)
        float valH = std::min(24.0f * uiScale, ph * 0.30f);
        float spacing = segCharSpacing(valH, uiScale);
        ImVec2 sz = segTextSize(valBuf, valH, spacing);

        const float maxW = std::max(1.0f, pw - 10.0f * uiScale);
        float scale = 1.0f;
        if (sz.x > maxW && sz.x > 1.0f) {
          scale = std::clamp(maxW / sz.x, 0.55f, 1.0f);
          valH *= scale;
          spacing = segCharSpacing(valH, uiScale);
          sz = segTextSize(valBuf, valH, spacing);
        }

        const float thBase = (plan.themeVariant == 0) ? 1.35f : (plan.themeVariant == 2 ? 2.10f : 1.70f);
        const float th = std::max(1.0f, thBase * uiScale * scale);

        // When glitch FX is active, drop a few segments deterministically per time-slice.
        float dropoutP = 0.0f;
        if (cfg.glitchFx) {
          dropoutP = std::clamp(glitchStrength * 0.42f, 0.0f, 0.28f);
        }

        core::SplitMix64 srng(core::hashCombine(panelSeed, (core::u64)std::floor(t.timeRealSec * 16.0)));

        const ImVec2 vPos = ImVec2(pMin.x + (pw - sz.x) * 0.5f, pMin.y + ph * 0.55f - sz.y * 0.5f);
        drawSegText(draw, vPos, valH, valBuf, gaugeCol, th, spacing, dropoutP, dropoutP > 0.0f ? &srng : nullptr);
      } else {
        ImVec2 vSize = ImGui::CalcTextSize(valBuf);
        ImVec2 vPos = ImVec2(pMin.x + (pw - vSize.x) * 0.5f, pMin.y + ph * 0.55f - vSize.y * 0.5f);
        draw->AddText(vPos, gaugeCol, valBuf);
      }
    }
  }

  // Global glitch overlay (scanlines/blocks)
  if (glitchStrength > 0.01f && cfg.glitchFx) {
    core::SplitMix64 grng(core::hashCombine(plan.seed, (core::u64)(t.timeRealSec * 60.0)));
    const int lines = (int)std::round(20.0f * glitchStrength);
    for (int i = 0; i < lines; ++i) {
      const float y = (float)grng.range<double>(winPos.y, winPos.y + winSize.y);
      const float jitter = (float)grng.range<double>(-18.0, 18.0) * uiScale * glitchStrength;
      const float a = 0.05f + 0.18f * glitchStrength;
      draw->AddLine(ImVec2(winPos.x + jitter, y), ImVec2(winPos.x + winSize.x + jitter, y), toU32(pal.grid, a), 1.0f * uiScale);
    }

    const int blocks = (int)std::round(12.0f * glitchStrength);
    for (int i = 0; i < blocks; ++i) {
      const float x = (float)grng.range<double>(winPos.x, winPos.x + winSize.x);
      const float y = (float)grng.range<double>(winPos.y, winPos.y + winSize.y);
      const float w = (float)grng.range<double>(12.0, 64.0) * uiScale;
      const float h = (float)grng.range<double>(2.0, 10.0) * uiScale;
      draw->AddRectFilled(ImVec2(x, y), ImVec2(x + w, y + h), toU32(pal.background, 0.12f * glitchStrength));
    }
  }
}

} // namespace stellar::ui
