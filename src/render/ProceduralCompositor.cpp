#include "stellar/render/ProceduralCompositor.h"

#include "stellar/core/Clamp.h"
#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <cmath>
#include <cstdio>
#include <string>

namespace stellar::render {

namespace {

static std::string glslFloat(float v) {
  char buf[64];
  // GLSL doesn't care about trailing zeros; keep it readable.
  std::snprintf(buf, sizeof(buf), "%.6f", (double)v);
  return std::string(buf);
}

static float signedUnit(core::SplitMix64& r) {
  // [-1,1]
  return (float)(r.nextUnit() * 2.0 - 1.0);
}

static float biased(core::SplitMix64& r, float minV, float maxV, float biasExp) {
  // biasExp>1 biases toward minV. biasExp<1 biases toward maxV.
  const double t = std::pow(r.nextUnit(), (double)biasExp);
  return (float)((double)minV + ((double)maxV - (double)minV) * t);
}

} // namespace

const char* compositorKindName(CompositorKind kind) {
  switch (kind) {
    case CompositorKind::Cinematic: return "Cinematic";
    case CompositorKind::Noir: return "Noir";
    case CompositorKind::Vaporwave: return "Vaporwave";
    case CompositorKind::Analog: return "Analog";
    case CompositorKind::Random: return "Random";
  }
  return "Random";
}

CompositorRecipe makeProceduralCompositorRecipe(core::u64 seed, CompositorKind kind) {
  // Make seed-space stable even if enum order changes.
  const core::u64 base = core::hashCombine(seed, core::fnv1a64("procCompositor"));
  core::SplitMix64 r(base);

  // If Random, pick a family first.
  CompositorKind family = kind;
  if (family == CompositorKind::Random) {
    const int k = r.range(0, 3);
    family = (k == 0) ? CompositorKind::Cinematic
            : (k == 1) ? CompositorKind::Noir
            : (k == 2) ? CompositorKind::Vaporwave
                       : CompositorKind::Analog;
  }

  CompositorRecipe out{};

  // Defaults (identity grade).
  out.tonemap = TonemapMode::Aces;
  out.saturation = 1.0f;
  out.contrast = 1.0f;
  out.hueShiftRad = 0.0f;
  out.liftR = out.liftG = out.liftB = 0.0f;
  out.gammaR = out.gammaG = out.gammaB = 1.0f;
  out.gainR = out.gainG = out.gainB = 1.0f;
  out.splitTone = false;
  out.splitStrength = 0.0f;
  out.shadowTintR = out.shadowTintG = out.shadowTintB = 1.0f;
  out.highlightTintR = out.highlightTintG = out.highlightTintB = 1.0f;
  out.vignetteMul = 1.0f;
  out.grainMul = 1.0f;

  // Family-specific recipe.
  if (family == CompositorKind::Cinematic) {
    out.tonemap = r.chance(0.50) ? TonemapMode::Aces : TonemapMode::Uncharted2;
    out.saturation = biased(r, 1.02f, 1.18f, 1.4f);
    out.contrast = biased(r, 1.03f, 1.16f, 1.2f);

    // Very small warm/cool bias via lift/gain.
    const float warm = biased(r, -0.02f, 0.03f, 1.5f);
    out.liftR = 0.00f + warm * 0.35f;
    out.liftG = 0.00f + warm * 0.10f;
    out.liftB = 0.00f - warm * 0.30f;
    out.gainR = 1.00f + warm * 0.25f;
    out.gainG = 1.00f + warm * 0.10f;
    out.gainB = 1.00f - warm * 0.15f;

    out.vignetteMul = biased(r, 0.90f, 1.25f, 1.2f);
    out.grainMul = biased(r, 0.80f, 1.35f, 1.3f);
  } else if (family == CompositorKind::Noir) {
    out.tonemap = TonemapMode::Reinhard;
    out.saturation = biased(r, 0.06f, 0.28f, 1.2f);
    out.contrast = biased(r, 1.10f, 1.45f, 1.1f);

    // Push highlights down a touch (crush), lift shadows slightly.
    const float lift = biased(r, 0.00f, 0.05f, 1.1f);
    out.liftR = out.liftG = out.liftB = lift;
    const float gain = biased(r, 0.88f, 1.00f, 1.6f);
    out.gainR = out.gainG = out.gainB = gain;
    const float gm = biased(r, 0.95f, 1.08f, 1.2f);
    out.gammaR = out.gammaG = out.gammaB = gm;

    out.vignetteMul = biased(r, 1.10f, 1.65f, 1.1f);
    out.grainMul = biased(r, 1.05f, 1.85f, 1.0f);
  } else if (family == CompositorKind::Vaporwave) {
    out.tonemap = TonemapMode::Aces;
    out.saturation = biased(r, 1.20f, 1.55f, 1.0f);
    out.contrast = biased(r, 0.95f, 1.25f, 1.2f);
    out.hueShiftRad = (float)(signedUnit(r) * 0.55); // ~ +/- 31 degrees

    out.splitTone = true;
    out.splitStrength = biased(r, 0.15f, 0.50f, 1.0f);

    // Shadows: magenta-ish. Highlights: cyan-ish.
    out.shadowTintR = biased(r, 1.05f, 1.55f, 1.0f);
    out.shadowTintG = biased(r, 0.65f, 1.05f, 1.1f);
    out.shadowTintB = biased(r, 1.05f, 1.45f, 1.0f);

    out.highlightTintR = biased(r, 0.75f, 1.05f, 1.2f);
    out.highlightTintG = biased(r, 1.05f, 1.55f, 1.0f);
    out.highlightTintB = biased(r, 1.05f, 1.55f, 1.0f);

    out.vignetteMul = biased(r, 0.75f, 1.10f, 1.4f);
    out.grainMul = biased(r, 0.55f, 1.05f, 1.4f);
  } else if (family == CompositorKind::Analog) {
    out.tonemap = r.chance(0.65) ? TonemapMode::Uncharted2 : TonemapMode::Aces;
    out.saturation = biased(r, 0.88f, 1.10f, 1.2f);
    out.contrast = biased(r, 0.95f, 1.20f, 1.2f);

    // Slight hue wander only.
    out.hueShiftRad = (float)(signedUnit(r) * 0.18);

    // Mild split tone with warm shadows.
    out.splitTone = r.chance(0.55);
    out.splitStrength = biased(r, 0.05f, 0.22f, 1.4f);
    out.shadowTintR = biased(r, 1.00f, 1.18f, 1.2f);
    out.shadowTintG = biased(r, 0.92f, 1.05f, 1.4f);
    out.shadowTintB = biased(r, 0.85f, 0.98f, 1.2f);
    out.highlightTintR = biased(r, 0.95f, 1.05f, 1.4f);
    out.highlightTintG = biased(r, 0.95f, 1.10f, 1.2f);
    out.highlightTintB = biased(r, 1.00f, 1.12f, 1.2f);

    out.vignetteMul = biased(r, 0.95f, 1.25f, 1.3f);
    out.grainMul = biased(r, 0.95f, 1.45f, 1.1f);
  }

  // Small safety clamps.
  out.saturation = (float)core::clamp(out.saturation, 0.0f, 2.0f);
  out.contrast = (float)core::clamp(out.contrast, 0.0f, 2.0f);
  out.vignetteMul = (float)core::clamp(out.vignetteMul, 0.0f, 3.0f);
  out.grainMul = (float)core::clamp(out.grainMul, 0.0f, 3.0f);
  out.splitStrength = (float)core::clamp(out.splitStrength, 0.0f, 1.0f);
  return out;
}

std::string buildCompositorShaderDefines(const CompositorRecipe& r) {
  std::string out;
  out.reserve(1024);

  // Tonemap selection.
  out += "#define STELLAR_TONEMAP_MODE ";
  out += std::to_string((int)r.tonemap);
  out += "\n";

  // Grade enabled if any control deviates from identity.
  const bool grade =
      std::abs(r.saturation - 1.0f) > 1e-6f ||
      std::abs(r.contrast - 1.0f) > 1e-6f ||
      std::abs(r.hueShiftRad) > 1e-6f ||
      std::abs(r.liftR) > 1e-6f || std::abs(r.liftG) > 1e-6f || std::abs(r.liftB) > 1e-6f ||
      std::abs(r.gammaR - 1.0f) > 1e-6f || std::abs(r.gammaG - 1.0f) > 1e-6f || std::abs(r.gammaB - 1.0f) > 1e-6f ||
      std::abs(r.gainR - 1.0f) > 1e-6f || std::abs(r.gainG - 1.0f) > 1e-6f || std::abs(r.gainB - 1.0f) > 1e-6f ||
      r.splitTone;

  out += "#define STELLAR_GRADE_ENABLED ";
  out += grade ? "1\n" : "0\n";

  // Grade constants.
  out += "#define STELLAR_GRADE_SAT "; out += glslFloat(r.saturation); out += "\n";
  out += "#define STELLAR_GRADE_CONTRAST "; out += glslFloat(r.contrast); out += "\n";
  out += "#define STELLAR_GRADE_HUE "; out += glslFloat(r.hueShiftRad); out += "\n";

  out += "#define STELLAR_GRADE_LIFT_R "; out += glslFloat(r.liftR); out += "\n";
  out += "#define STELLAR_GRADE_LIFT_G "; out += glslFloat(r.liftG); out += "\n";
  out += "#define STELLAR_GRADE_LIFT_B "; out += glslFloat(r.liftB); out += "\n";

  out += "#define STELLAR_GRADE_GAMMA_R "; out += glslFloat(r.gammaR); out += "\n";
  out += "#define STELLAR_GRADE_GAMMA_G "; out += glslFloat(r.gammaG); out += "\n";
  out += "#define STELLAR_GRADE_GAMMA_B "; out += glslFloat(r.gammaB); out += "\n";

  out += "#define STELLAR_GRADE_GAIN_R "; out += glslFloat(r.gainR); out += "\n";
  out += "#define STELLAR_GRADE_GAIN_G "; out += glslFloat(r.gainG); out += "\n";
  out += "#define STELLAR_GRADE_GAIN_B "; out += glslFloat(r.gainB); out += "\n";

  out += "#define STELLAR_SPLIT_TONE_ENABLED "; out += (r.splitTone ? "1\n" : "0\n");
  out += "#define STELLAR_SPLIT_TONE_STRENGTH "; out += glslFloat(r.splitStrength); out += "\n";
  out += "#define STELLAR_SPLIT_SHADOW_R "; out += glslFloat(r.shadowTintR); out += "\n";
  out += "#define STELLAR_SPLIT_SHADOW_G "; out += glslFloat(r.shadowTintG); out += "\n";
  out += "#define STELLAR_SPLIT_SHADOW_B "; out += glslFloat(r.shadowTintB); out += "\n";
  out += "#define STELLAR_SPLIT_HI_R "; out += glslFloat(r.highlightTintR); out += "\n";
  out += "#define STELLAR_SPLIT_HI_G "; out += glslFloat(r.highlightTintG); out += "\n";
  out += "#define STELLAR_SPLIT_HI_B "; out += glslFloat(r.highlightTintB); out += "\n";

  // Multipliers for existing uniforms.
  out += "#define STELLAR_VIGNETTE_MUL "; out += glslFloat(r.vignetteMul); out += "\n";
  out += "#define STELLAR_GRAIN_MUL "; out += glslFloat(r.grainMul); out += "\n";

  return out;
}

} // namespace stellar::render
