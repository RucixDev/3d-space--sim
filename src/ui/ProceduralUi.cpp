#include "stellar/ui/ProceduralUi.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::ui {

// Classic HSV (H in degrees [0..360), S/V in [0..1]).
// (Copied/ported in a lightweight form from common reference implementations.
//  We keep this local to UI to avoid introducing a hard dependency on render code.)
static void hsvToRgb(float hDeg, float s, float v, float& outR, float& outG, float& outB) {
  hDeg = std::fmod(hDeg, 360.0f);
  if (hDeg < 0.0f) hDeg += 360.0f;
  s = std::clamp(s, 0.0f, 1.0f);
  v = std::clamp(v, 0.0f, 1.0f);

  if (s <= 0.0f) {
    outR = v;
    outG = v;
    outB = v;
    return;
  }

  const float h = hDeg / 60.0f;
  const int i = static_cast<int>(std::floor(h)) % 6;
  const float f = h - std::floor(h);
  const float p = v * (1.0f - s);
  const float q = v * (1.0f - s * f);
  const float t = v * (1.0f - s * (1.0f - f));

  switch (i) {
    default:
    case 0: outR = v; outG = t; outB = p; break;
    case 1: outR = q; outG = v; outB = p; break;
    case 2: outR = p; outG = v; outB = t; break;
    case 3: outR = p; outG = q; outB = v; break;
    case 4: outR = t; outG = p; outB = v; break;
    case 5: outR = v; outG = p; outB = q; break;
  }
}

static Color4f rgb(float r, float g, float b, float a = 1.0f) {
  return Color4f{std::clamp(r, 0.0f, 1.0f), std::clamp(g, 0.0f, 1.0f), std::clamp(b, 0.0f, 1.0f), std::clamp(a, 0.0f, 1.0f)};
}

static Color4f mulRgb(const Color4f& c, float k, float a = 1.0f) {
  return rgb(c.r * k, c.g * k, c.b * k, c.a * a);
}

static Color4f lerpRgb(const Color4f& a, const Color4f& b, float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  return rgb(a.r + (b.r - a.r) * t,
             a.g + (b.g - a.g) * t,
             a.b + (b.b - a.b) * t,
             a.a + (b.a - a.a) * t);
}

ProceduralHudSkin makeProceduralHudSkin(core::u64 seed) {
  // Split seed streams so changing the layout gen doesn't reshuffle the palette.
  core::SplitMix64 rng(core::hashCombine(seed, core::fnv1a64("hud_skin")));

  ProceduralHudSkin out{};

  // Choose a "family" first; these are intentionally biased toward sci-fi cockpit palettes.
  // 0: cyan/blue, 1: green, 2: amber, 3: red, 4: monochrome
  const int family = rng.range<int>(0, 4);

  float hue = 200.0f;
  float sat = 0.35f;
  float val = 0.95f;

  switch (family) {
    default:
    case 0: // Cyan / blue
      hue = rng.range<float>(186.0f, 220.0f);
      sat = rng.range<float>(0.25f, 0.55f);
      val = rng.range<float>(0.86f, 1.00f);
      break;
    case 1: // Green phosphor
      hue = rng.range<float>(110.0f, 145.0f);
      sat = rng.range<float>(0.45f, 0.90f);
      val = rng.range<float>(0.82f, 1.00f);
      break;
    case 2: // Amber
      hue = rng.range<float>(30.0f, 62.0f);
      sat = rng.range<float>(0.50f, 0.95f);
      val = rng.range<float>(0.84f, 1.00f);
      break;
    case 3: // Red alert
      // Use two bands to avoid only-pink reds.
      if (rng.chance(0.5)) hue = rng.range<float>(0.0f, 18.0f);
      else hue = rng.range<float>(342.0f, 360.0f);
      sat = rng.range<float>(0.55f, 1.00f);
      val = rng.range<float>(0.80f, 1.00f);
      break;
    case 4: // Monochrome
      hue = 0.0f;
      sat = 0.0f;
      val = rng.range<float>(0.86f, 0.98f);
      break;
  }

  // Primary color.
  {
    float r = 1, g = 1, b = 1;
    hsvToRgb(hue, sat, val, r, g, b);
    out.colorPrimary = rgb(r, g, b, 1.0f);
  }

  // Accent: complementary-ish hue, but keep it in a visually distinct band.
  {
    float ah = hue;
    if (family == 0) ah = rng.range<float>(18.0f, 55.0f);           // cyan -> orange
    else if (family == 1) ah = rng.range<float>(36.0f, 65.0f);      // green -> amber
    else if (family == 2) ah = rng.range<float>(180.0f, 220.0f);    // amber -> cyan
    else if (family == 3) ah = rng.range<float>(36.0f, 65.0f);      // red -> amber
    else if (family == 4) ah = rng.range<float>(0.0f, 360.0f);      // mono -> subtle tint

    const float as = (family == 4) ? rng.range<float>(0.04f, 0.15f) : rng.range<float>(0.55f, 0.95f);
    const float av = rng.range<float>(0.78f, 1.00f);
    float r = 1, g = 1, b = 1;
    hsvToRgb(ah, as, av, r, g, b);
    out.colorAccent = rgb(r, g, b, 1.0f);
  }

  // Danger: red, slightly varied.
  {
    float dh = rng.chance(0.5) ? rng.range<float>(0.0f, 18.0f) : rng.range<float>(342.0f, 360.0f);
    float ds = rng.range<float>(0.65f, 1.00f);
    float dv = rng.range<float>(0.75f, 1.00f);
    float r = 1, g = 1, b = 1;
    hsvToRgb(dh, ds, dv, r, g, b);
    out.colorDanger = rgb(r, g, b, 1.0f);
  }

  // Background always black; we use alpha masks + background alpha for overlays.
  out.colorBackground = rgb(0.0f, 0.0f, 0.0f, 1.0f);

  // Grid & text are derived from primary to keep the palette cohesive.
  out.colorGrid = lerpRgb(out.colorPrimary, out.colorBackground, rng.range<float>(0.35f, 0.55f));
  out.colorText = lerpRgb(rgb(0.94f, 0.94f, 0.97f, 1.0f), out.colorPrimary, rng.range<float>(0.08f, 0.22f));

  // Cosmetics.
  out.overlayBgAlpha = rng.range<float>(0.22f, 0.55f);
  out.overlayBgAlphaEdit = std::clamp(out.overlayBgAlpha + rng.range<float>(0.05f, 0.22f), 0.0f, 1.0f);

  // Tinting reads nicely for phosphor/amber; for monochrome it can be too flat.
  if (family == 4) {
    out.tintRadarIcons = false;
    out.tintTacticalIcons = false;
  } else {
    out.tintRadarIcons = rng.chance(0.70);
    out.tintTacticalIcons = rng.chance(0.60);
  }

  // Ensure palette isn't too dim (helps readability).
  out.colorPrimary = mulRgb(out.colorPrimary, rng.range<float>(0.95f, 1.05f));
  out.colorGrid = mulRgb(out.colorGrid, rng.range<float>(0.90f, 1.05f));

  return out;
}

void applyHudSkinToSettings(const ProceduralHudSkin& skin, HudSettings& inOut) {
  inOut.overlayBgAlpha = skin.overlayBgAlpha;
  inOut.overlayBgAlphaEdit = skin.overlayBgAlphaEdit;
  inOut.tintRadarIcons = skin.tintRadarIcons;
  inOut.tintTacticalIcons = skin.tintTacticalIcons;
  inOut.colorPrimary = skin.colorPrimary;
  inOut.colorAccent = skin.colorAccent;
  inOut.colorDanger = skin.colorDanger;
  inOut.colorGrid = skin.colorGrid;
  inOut.colorText = skin.colorText;
  inOut.colorBackground = skin.colorBackground;
}

void applyProceduralHudLayout(core::u64 seed, HudLayout& inOut) {
  core::SplitMix64 rng(core::hashCombine(seed, core::fnv1a64("hud_layout")));

  // Preserve user tuning.
  const float radarScale = inOut.widget(HudWidgetId::Radar).scale;
  const float objScale = inOut.widget(HudWidgetId::Objective).scale;
  const float threatScale = inOut.widget(HudWidgetId::Threat).scale;
  const float jumpScale = inOut.widget(HudWidgetId::Jump).scale;
  const bool radarEnabled = inOut.widget(HudWidgetId::Radar).enabled;
  const bool objEnabled = inOut.widget(HudWidgetId::Objective).enabled;
  const bool threatEnabled = inOut.widget(HudWidgetId::Threat).enabled;
  const bool jumpEnabled = inOut.widget(HudWidgetId::Jump).enabled;

  // Start from defaults, then procedurally vary positions.
  HudLayout l = makeDefaultHudLayout();
  l.safeMarginPx = inOut.safeMarginPx;

  // A tiny normalized inset (safe across resolutions). Keep it small so users can further drag.
  const float m = rng.range<float>(0.015f, 0.035f);
  const float left = m;
  const float right = 1.0f - m;
  const float top = m;
  const float bottom = 1.0f - m;

  const bool flipX = rng.chance(0.5);
  const int variant = rng.range<int>(0, 2); // 0 classic, 1 stacked right, 2 stacked left

  auto placeCorner = [&](HudWidgetId id, float x, float y, float pivotX, float pivotY) {
    auto& w = l.widget(id);
    w.posNormX = std::clamp(x, 0.0f, 1.0f);
    w.posNormY = std::clamp(y, 0.0f, 1.0f);
    w.pivotX = std::clamp(pivotX, 0.0f, 1.0f);
    w.pivotY = std::clamp(pivotY, 0.0f, 1.0f);
  };

  // Base placements.
  if (variant == 0) {
    // Classic: radar bottom corner, objective + jump at top, threat opposite top.
    if (!flipX) {
      placeCorner(HudWidgetId::Radar, right, bottom, 1.0f, 1.0f);
      placeCorner(HudWidgetId::Objective, right, top, 1.0f, 0.0f);
      placeCorner(HudWidgetId::Threat, left, top, 0.0f, 0.0f);
    } else {
      placeCorner(HudWidgetId::Radar, left, bottom, 0.0f, 1.0f);
      placeCorner(HudWidgetId::Objective, left, top, 0.0f, 0.0f);
      placeCorner(HudWidgetId::Threat, right, top, 1.0f, 0.0f);
    }
    placeCorner(HudWidgetId::Jump, 0.5f + rng.range<float>(-0.03f, 0.03f), top, 0.5f, 0.0f);
  } else if (variant == 1) {
    // Stack on right: objective + jump + threat live on right edge.
    placeCorner(HudWidgetId::Radar, flipX ? left : right, bottom, flipX ? 0.0f : 1.0f, 1.0f);
    const float x = flipX ? left : right;
    const float px = flipX ? 0.0f : 1.0f;
    placeCorner(HudWidgetId::Objective, x, top, px, 0.0f);
    placeCorner(HudWidgetId::Jump, x, top + rng.range<float>(0.10f, 0.16f), px, 0.0f);
    placeCorner(HudWidgetId::Threat, x, top + rng.range<float>(0.22f, 0.30f), px, 0.0f);
  } else {
    // Stack on left: objective + jump + threat live on left edge.
    placeCorner(HudWidgetId::Radar, flipX ? right : left, bottom, flipX ? 1.0f : 0.0f, 1.0f);
    const float x = flipX ? right : left;
    const float px = flipX ? 1.0f : 0.0f;
    placeCorner(HudWidgetId::Objective, x, top, px, 0.0f);
    placeCorner(HudWidgetId::Jump, x, top + rng.range<float>(0.10f, 0.16f), px, 0.0f);
    placeCorner(HudWidgetId::Threat, x, top + rng.range<float>(0.22f, 0.30f), px, 0.0f);
  }

  // Restore user tuning.
  l.widget(HudWidgetId::Radar).scale = radarScale;
  l.widget(HudWidgetId::Objective).scale = objScale;
  l.widget(HudWidgetId::Threat).scale = threatScale;
  l.widget(HudWidgetId::Jump).scale = jumpScale;
  l.widget(HudWidgetId::Radar).enabled = radarEnabled;
  l.widget(HudWidgetId::Objective).enabled = objEnabled;
  l.widget(HudWidgetId::Threat).enabled = threatEnabled;
  l.widget(HudWidgetId::Jump).enabled = jumpEnabled;

  // Clamp positions to the [0..1] range (paranoia; should already be valid).
  for (auto& w : l.widgets) {
    w.posNormX = std::clamp(w.posNormX, 0.0f, 1.0f);
    w.posNormY = std::clamp(w.posNormY, 0.0f, 1.0f);
    w.pivotX = std::clamp(w.pivotX, 0.0f, 1.0f);
    w.pivotY = std::clamp(w.pivotY, 0.0f, 1.0f);
  }

  inOut = l;
}

HudLayout makeProceduralHudLayout(core::u64 seed) {
  HudLayout l = makeDefaultHudLayout();
  applyProceduralHudLayout(seed, l);
  return l;
}

void makeProceduralAccentRgb(core::u64 seed, float outRgb[3]) {
  if (!outRgb) return;
  const ProceduralHudSkin skin = makeProceduralHudSkin(seed);
  outRgb[0] = skin.colorPrimary.r;
  outRgb[1] = skin.colorPrimary.g;
  outRgb[2] = skin.colorPrimary.b;
}

} // namespace stellar::ui
