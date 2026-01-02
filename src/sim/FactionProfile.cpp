#include "stellar/sim/FactionProfile.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

FactionProfile factionProfile(core::u64 universeSeed, core::u32 factionId) {
  FactionProfile p{};
  p.factionId = factionId;

  if (factionId == 0) {
    // Independent / neutral baseline.
    return p;
  }

  // IMPORTANT: authority/corruption are derived from the same seed as law_profile_v1,
  // so any gameplay code using LawProfile remains consistent with these traits.
  {
    core::u64 s = core::hashCombine(universeSeed, static_cast<core::u64>(factionId));
    s = core::hashCombine(s, core::seedFromText("law_profile_v1"));
    core::SplitMix64 rng(s);

    const double authority = clamp01(rng.nextDouble());
    const double corruptionRaw = clamp01(rng.nextDouble());

    p.authority = authority;
    p.corruption = (0.25 + 0.75 * corruptionRaw) * (1.05 - 0.70 * authority);
    p.corruption = std::clamp(p.corruption, 0.05, 1.0);
  }

  // Extra flavor traits. These use a separate seed so they can evolve independently
  // without disturbing legacy determinism in other proc systems.
  {
    core::u64 s = core::hashCombine(universeSeed, static_cast<core::u64>(factionId));
    s = core::hashCombine(s, core::seedFromText("faction_profile_v1"));
    core::SplitMix64 rng(s);

    const double w = clamp01(rng.nextDouble());
    const double tNoise = clamp01(rng.nextDouble());
    const double mNoise = clamp01(rng.nextDouble());
    const double sNoise = clamp01(rng.nextDouble());

    p.wealth = w;

    // Tech correlates with wealth, but isn't identical.
    p.tech = clamp01(0.20 + 0.80 * (0.65 * w + 0.35 * tNoise));

    // Militarism loosely correlates with authority.
    p.militarism = clamp01(0.15 + 0.85 * (0.55 * p.authority + 0.45 * mNoise));

    // Stability correlates with wealth and low corruption, with some authority influence.
    const double stabilityCore = 0.45 * p.authority + 0.35 * w + 0.20 * (1.0 - p.corruption);
    p.stability = clamp01(0.10 + 0.90 * stabilityCore * (0.85 + 0.30 * sNoise));
  }

  return p;
}

FactionRelationKind classifyFactionRelation(double relation) {
  if (relation <= -0.35) return FactionRelationKind::Hostile;
  if (relation >= 0.35) return FactionRelationKind::Allied;
  return FactionRelationKind::Neutral;
}

const char* factionRelationKindName(FactionRelationKind k) {
  switch (k) {
    case FactionRelationKind::Hostile: return "Hostile";
    case FactionRelationKind::Neutral: return "Neutral";
    case FactionRelationKind::Allied:  return "Allied";
    default: return "Neutral";
  }
}

// A deterministic, symmetric relation score between two factions.
//
// This is intentionally lightweight and can be used during proc/gen, mission scoring,
// or economy heuristics without needing a precomputed matrix.
double factionRelation(core::u64 universeSeed, const Faction& a, const Faction& b) {
  if (a.id == b.id) return 1.0;

  const FactionProfile pa = factionProfile(universeSeed, a.id);
  const FactionProfile pb = factionProfile(universeSeed, b.id);

  const auto absd = [](double x, double y) { return std::abs(x - y); };

  const double traitDist =
      (absd(pa.authority, pb.authority) +
       absd(pa.corruption, pb.corruption) +
       absd(pa.wealth, pb.wealth) +
       absd(pa.stability, pb.stability) +
       absd(pa.tech, pb.tech) +
       absd(pa.militarism, pb.militarism)) / 6.0;

  const double traitSim = clamp01(1.0 - traitDist);

  // Border pressure: factions with overlapping spheres of influence tend to compete.
  const math::Vec3d dv = a.homePosLy - b.homePosLy;
  const double distLy = std::sqrt(dv.lengthSq());
  const double influenceLy = std::max(1e-6, 0.5 * (a.influenceRadiusLy + b.influenceRadiusLy));
  const double border = std::exp(-distLy / influenceLy); // close => ~1, far => ~0

  // Base relation: similarity helps, overlap hurts.
  double rel = traitSim - border; // [-1, +1]

  // Small deterministic noise per *unordered* pair so the galaxy isn't too uniform.
  core::u64 s = core::hashCombine(universeSeed, core::seedFromText("faction_relation_v1"));
  const core::u32 lo = std::min(a.id, b.id);
  const core::u32 hi = std::max(a.id, b.id);
  const core::u64 key = (static_cast<core::u64>(lo) << 32) | static_cast<core::u64>(hi);
  s = core::hashCombine(s, key);
  core::SplitMix64 rng(s);

  const double noise = (rng.nextDouble() - 0.5) * 0.20; // [-0.10, +0.10]
  rel = std::clamp(rel + noise, -1.0, 1.0);

  return rel;
}

} // namespace stellar::sim
