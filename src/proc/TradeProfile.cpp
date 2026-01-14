#include "stellar/proc/TradeProfile.h"

#include "stellar/core/Hash.h"
#include "stellar/core/StableHash.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace stellar::proc {

static inline double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

static inline double hash01(core::u64 h) {
  // Map to [0,1). (Same technique used in Noise.cpp.)
  return (double)((h >> 11) & ((1ull << 53) - 1)) / (double)(1ull << 53);
}

static double starHabitability01(sim::StarClass c) {
  // Very coarse “habitable potential” curve. This is intentionally stylized.
  switch (c) {
    case sim::StarClass::G: return 1.00;
    case sim::StarClass::F: return 0.90;
    case sim::StarClass::K: return 0.85;
    case sim::StarClass::M: return 0.70;
    case sim::StarClass::A: return 0.60;
    case sim::StarClass::B: return 0.40;
    case sim::StarClass::O: return 0.25;
    default: return 0.75;
  }
}

TradeProfile generateTradeProfile(core::u64 universeSeed, const sim::SystemStub& stub) {
  TradeProfile p{};

  // Galaxy-scale coordinates (ly). We feed these into low-frequency noise fields.
  const double x = stub.posLy.x;
  const double y = stub.posLy.y;
  const double z = stub.posLy.z;

  const double r = std::sqrt(x * x + y * y);

  // These constants mirror defaults in GalaxyParams.
  // (We keep this module independent of GalaxyGenerator so UI tools can use it without a generator instance.)
  constexpr double kRadiusLy = 50'000.0;
  constexpr double kRadialScaleLy = 15'000.0;

  // “Galactic core” bias: high near the center, tapers outward.
  const double coreBias = std::exp(-r / kRadialScaleLy);
  const double rimBias = 1.0 - coreBias;

  // Thickness weighting (downweight far-off-plane systems for pop/wealth).
  constexpr double kHalfThicknessLy = 500.0; // ~thicknessLy/2 default
  const double zBias = clamp01(1.0 - (std::abs(z) / kHalfThicknessLy));

  const double hab = starHabitability01(stub.primaryClass);

  // Derive a handful of coherent noise fields. Use only universeSeed for coherence,
  // and use stub.seed for light micro-variation.
  const core::u64 macroSeed = core::hashCombine(universeSeed, core::fnv1a64("trade_profile_macro"));
  const core::u64 microSeed = core::hashCombine(stub.seed, core::fnv1a64("trade_profile_micro"));

  // Macro frequency: ~1000 ly features.
  constexpr double kF = 0.0010;
  constexpr double kF2 = 0.00065;

  const double nRes = ridgedFbmPerlin2D(macroSeed + 11, x * kF, y * kF, 5);
  const double nAg  = fbmPerlin2D(macroSeed + 23, x * kF, y * kF, 5);
  const double nHub = fbmPerlin2D(macroSeed + 37, x * kF2, y * kF2, 4);
  const double nInd = fbmPerlin2D(macroSeed + 41, x * kF, y * kF, 4);
  const double nTec = fbmPerlin2D(macroSeed + 53, x * kF, y * kF, 4);
  const double nWea = fbmPerlin2D(macroSeed + 59, x * kF2, y * kF2, 4);
  const double nLaw = fbmPerlin2D(macroSeed + 67, x * kF, y * kF, 4);
  const double nPop = fbmPerlin2D(macroSeed + 71, x * kF2, y * kF2, 4);

  // Micro frequency: ~50-100 ly “local flavor”.
  constexpr double kMF = 0.015;
  const double micro = perlin3D(microSeed + 3, x * kMF, y * kMF, z * kMF);

  // ---- Macro factors ----
  // Hub: grows near core + in high-nHub corridors.
  p.hub = clamp01(0.25 * coreBias + 0.75 * nHub);

  // Population: habitable + hub + near-plane.
  p.population = clamp01((0.25 * hab + 0.55 * p.hub + 0.20 * nPop) * (0.55 + 0.45 * zBias));

  // Technology & wealth: core + hub-driven.
  p.technology = clamp01(0.55 * coreBias + 0.25 * p.hub + 0.20 * nTec);
  p.wealth     = clamp01(0.35 * p.technology + 0.40 * p.hub + 0.25 * nWea);

  // Resources: ridged noise + rim bias.
  p.resources = clamp01(0.65 * nRes + 0.25 * rimBias + 0.10 * (1.0 - hab));

  // Agriculture: strongly gated by habitability.
  p.agriculture = clamp01((0.55 * nAg + 0.45 * hab) * (0.60 + 0.40 * zBias));

  // Industry: supply-chain: resources + hub + noise.
  p.industry = clamp01(0.30 * p.resources + 0.45 * p.hub + 0.25 * nInd);

  // Faction stability influences lawlessness slightly (purely procedural; does not consult faction data).
  const core::u64 fH = core::hashCombine(universeSeed, core::hashCombine(0xBADC0FFEEull, (core::u64)stub.factionId));
  const double factionStability = 0.35 + 0.65 * hash01(fH);

  // Lawlessness: higher at the rim, lower in wealthy/hub systems.
  // nLaw provides regional “pirate belts”.
  const double orderBias = (0.60 * p.hub + 0.40 * p.wealth) * factionStability;
  p.lawlessness = clamp01(0.55 * rimBias + 0.25 * nLaw + 0.20 * (1.0 - orderBias));

  // Inject a tiny bit of micro variation so adjacent systems don't look identical.
  const double microSigned = (micro - 0.5) * 2.0;
  p.resources    = clamp01(p.resources    + 0.05 * microSigned);
  p.agriculture  = clamp01(p.agriculture  + 0.04 * microSigned);
  p.industry     = clamp01(p.industry     + 0.04 * microSigned);
  p.technology   = clamp01(p.technology   + 0.03 * microSigned);
  p.hub          = clamp01(p.hub          + 0.03 * microSigned);
  p.wealth       = clamp01(p.wealth       + 0.03 * microSigned);
  p.population   = clamp01(p.population   + 0.03 * microSigned);
  p.lawlessness  = clamp01(p.lawlessness  - 0.03 * microSigned);

  // ---- Commodity mapping ----
  const auto addScores = [&](econ::CommodityId id, double expBase, double impBase) {
    const std::size_t i = static_cast<std::size_t>(id);

    // Small commodity-specific modulation (keeps profiles from looking “too templated”).
    constexpr double kCF = 0.0035;
    const double eN = perlin2D(macroSeed + 900 + (core::u64)i * 19ull, x * kCF, y * kCF);
    const double iN = perlin2D(macroSeed + 1200 + (core::u64)i * 23ull, x * kCF, y * kCF);

    const double exp = clamp01(expBase + 0.22 * (eN - 0.5));
    const double imp = clamp01(impBase + 0.22 * (iN - 0.5));

    p.exportScore[i] = exp;
    p.importScore[i] = imp;
  };

  const double res = p.resources;
  const double ag  = p.agriculture;
  const double ind = p.industry;
  const double tec = p.technology;
  const double hub = p.hub;
  const double law = p.lawlessness;
  const double pop = p.population;
  const double wea = p.wealth;

  // Food / water.
  addScores(econ::CommodityId::Food,
            /*export*/ clamp01(ag * (1.0 - 0.35 * pop) + 0.10 * (1.0 - law)),
            /*import*/ clamp01(pop * (1.0 - ag)));

  addScores(econ::CommodityId::Water,
            /*export*/ clamp01(0.45 * ag + 0.25 * res + 0.10 * (1.0 - law) + 0.20 * zBias),
            /*import*/ clamp01(pop * (1.0 - (0.55 * ag + 0.25 * res))));

  // Raw materials.
  addScores(econ::CommodityId::Ore,
            /*export*/ clamp01(res * (1.0 - 0.45 * ind)),
            /*import*/ clamp01(ind * (1.0 - res)));

  addScores(econ::CommodityId::Metals,
            /*export*/ clamp01(0.70 * res + 0.20 * ind),
            /*import*/ clamp01(ind * (1.0 - res) + 0.15 * tec));

  addScores(econ::CommodityId::Fuel,
            /*export*/ clamp01(0.50 * res + 0.20 * tec + 0.10 * law + 0.20 * rimBias),
            /*import*/ clamp01(0.45 * hub + 0.30 * ind + 0.25 * pop));

  // Manufactured goods.
  addScores(econ::CommodityId::Machinery,
            /*export*/ clamp01(0.75 * ind + 0.15 * tec),
            /*import*/ clamp01(0.30 * ag + 0.25 * res + 0.25 * pop * (1.0 - ind)));

  addScores(econ::CommodityId::Medicine,
            /*export*/ clamp01(0.55 * tec + 0.25 * wea + 0.10 * hub),
            /*import*/ clamp01(0.60 * pop * (1.0 - tec) + 0.20 * law));

  addScores(econ::CommodityId::Electronics,
            /*export*/ clamp01(0.65 * tec + 0.20 * ind + 0.10 * wea),
            /*import*/ clamp01(0.45 * hub * (1.0 - tec) + 0.35 * wea * (1.0 - tec)));

  addScores(econ::CommodityId::Luxury,
            /*export*/ clamp01(0.65 * wea + 0.25 * hub + 0.10 * tec),
            /*import*/ clamp01(0.60 * wea + 0.25 * hub + 0.15 * pop));

  // Specialty / illicit-ish.
  addScores(econ::CommodityId::Weapons,
            /*export*/ clamp01(0.45 * ind + 0.20 * tec + 0.35 * law),
            /*import*/ clamp01(0.55 * law + 0.25 * hub + 0.20 * pop * (1.0 - ind)));

  addScores(econ::CommodityId::Stimulants,
            /*export*/ clamp01(0.50 * law + 0.30 * hub + 0.20 * pop),
            /*import*/ clamp01(0.55 * law + 0.35 * pop + 0.10 * hub));

  // ---- Signature ----
  // Quantize to reduce floating noise and make signatures stable under minor math changes.
  auto q16 = [](double v) -> core::u32 {
    const double c = clamp01(v);
    return static_cast<core::u32>(std::llround(c * 65535.0));
  };

  core::StableHash64 h;
  h.addU64(q16(p.resources));
  h.addU64(q16(p.agriculture));
  h.addU64(q16(p.industry));
  h.addU64(q16(p.technology));
  h.addU64(q16(p.hub));
  h.addU64(q16(p.lawlessness));
  h.addU64(q16(p.population));
  h.addU64(q16(p.wealth));

  for (std::size_t i = 0; i < econ::kCommodityCount; ++i) {
    h.addU64(q16(p.exportScore[i]));
    h.addU64(q16(p.importScore[i]));
  }

  // Include some discrete identifiers so two different stubs in the same spot
  // (unlikely, but possible in tools) don't collide.
  h.addU64(stub.seed);
  h.addU64(stub.id);
  p.signature = h.value();

  return p;
}

} // namespace stellar::proc
