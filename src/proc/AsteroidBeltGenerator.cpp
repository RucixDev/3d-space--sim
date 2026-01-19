#include "stellar/proc/AsteroidBeltGenerator.h"

#include "stellar/core/Hash.h"
#include "stellar/core/LowDiscrepancy.h"
#include "stellar/core/Random.h"
#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace stellar::proc {

namespace {

constexpr double kTwoPi = 2.0 * stellar::math::kPi;

static inline double clamp01(double v) { return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v); }

static inline double lerp(double a, double b, double t) { return a + (b - a) * t; }

static inline double smoothstep(double e0, double e1, double x) {
  if (e0 == e1) return (x < e0) ? 0.0 : 1.0;
  const double t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3.0 - 2.0 * t);
}

static inline double wrapAnglePi(double a) {
  a = std::fmod(a + stellar::math::kPi, kTwoPi);
  if (a < 0.0) a += kTwoPi;
  return a - stellar::math::kPi;
}

static inline double gauss01(double x, double sigma) {
  const double s = std::max(1e-9, sigma);
  const double t = x / s;
  return std::exp(-0.5 * t * t);
}

static math::Vec3d orbitNormal(const sim::OrbitElements& o) {
  const double i = o.inclinationRad;
  const double O = o.ascendingNodeRad;
  const double si = std::sin(i);
  const double ci = std::cos(i);
  return {std::sin(O) * si, -std::cos(O) * si, ci};
}

static void buildOrthonormalBasis(const math::Vec3d& nIn,
                                  math::Vec3d& outX,
                                  math::Vec3d& outY,
                                  math::Vec3d& outZ) {
  math::Vec3d n = nIn.normalized();
  if (n.lengthSq() < 1e-12) n = {0.0, 1.0, 0.0};

  // Pick an axis not parallel to n.
  const math::Vec3d a = (std::abs(n.y) < 0.90) ? math::Vec3d{0.0, 1.0, 0.0} : math::Vec3d{1.0, 0.0, 0.0};
  math::Vec3d x = math::cross(a, n).normalized();
  if (x.lengthSq() < 1e-12) x = math::cross(math::Vec3d{0.0, 0.0, 1.0}, n).normalized();
  math::Vec3d z = math::cross(n, x).normalized();

  outX = x;
  outY = n;
  outZ = z;
}

static math::Vec3d meanPlanetPlaneNormal(const sim::StarSystem& sys) {
  if (sys.planets.empty()) {
    return {0.0, 1.0, 0.0};
  }

  math::Vec3d acc{0.0, 0.0, 0.0};
  double wsum = 0.0;

  // Weight inner planets slightly higher (they tend to define the "ecliptic" feel),
  // while still including the giant planets.
  for (const auto& p : sys.planets) {
    const double a = std::max(0.05, p.orbit.semiMajorAxisAU);
    const double w = 1.0 / std::sqrt(a);
    const math::Vec3d n = orbitNormal(p.orbit);
    acc += n * w;
    wsum += w;
  }

  if (wsum <= 0.0 || acc.lengthSq() < 1e-12) {
    return {0.0, 1.0, 0.0};
  }
  return (acc * (1.0 / wsum)).normalized();
}

static int pickDominantPlanet(const sim::StarSystem& sys) {
  // Prefer the most massive gas giant; else most massive planet.
  int best = -1;
  double bestMass = -1.0;

  for (std::size_t i = 0; i < sys.planets.size(); ++i) {
    const auto& p = sys.planets[i];
    if (p.type != sim::PlanetType::GasGiant) continue;
    if (p.massEarth > bestMass) {
      bestMass = p.massEarth;
      best = static_cast<int>(i);
    }
  }

  if (best >= 0) return best;

  for (std::size_t i = 0; i < sys.planets.size(); ++i) {
    const auto& p = sys.planets[i];
    if (p.massEarth > bestMass) {
      bestMass = p.massEarth;
      best = static_cast<int>(i);
    }
  }
  return best;
}

static int pickInnerGiant(const sim::StarSystem& sys) {
  int best = -1;
  double bestA = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < sys.planets.size(); ++i) {
    const auto& p = sys.planets[i];
    if (p.type != sim::PlanetType::GasGiant) continue;
    if (p.orbit.semiMajorAxisAU < bestA) {
      bestA = p.orbit.semiMajorAxisAU;
      best = static_cast<int>(i);
    }
  }
  return best;
}

static int pickOuterMostPlanet(const sim::StarSystem& sys) {
  int best = -1;
  double bestA = -1.0;
  for (std::size_t i = 0; i < sys.planets.size(); ++i) {
    const auto& p = sys.planets[i];
    if (p.orbit.semiMajorAxisAU > bestA) {
      bestA = p.orbit.semiMajorAxisAU;
      best = static_cast<int>(i);
    }
  }
  return best;
}

static double resonanceRadiusAU(double aPlanetAU, int m, int n, bool interior) {
  if (m <= 0 || n <= 0) return aPlanetAU;
  const double ratio = interior ? (static_cast<double>(n) / static_cast<double>(m))
                                : (static_cast<double>(m) / static_cast<double>(n));
  // Kepler: P ~ a^(3/2) => a ~ P^(2/3)
  return aPlanetAU * std::pow(std::max(1e-9, ratio), 2.0 / 3.0);
}

static void addKirkwoodGaps(AsteroidBelt& belt,
                            const sim::Planet& planet,
                            core::SplitMix64& rng) {
  // A small curated set of classic interior resonances that tend to carve visible gaps.
  // (planet:body period ratio m:n)
  struct R { int m; int n; double base; };
  const R rs[] = {
      {2, 1, 0.55},
      {3, 1, 0.80},
      {5, 2, 0.45},
      {7, 3, 0.32},
      {4, 1, 0.22},
  };

  const double aP = std::max(0.05, planet.orbit.semiMajorAxisAU);
  const double wB = std::max(1e-6, belt.aOuterAU - belt.aInnerAU);

  for (const auto& r : rs) {
    const double aR = resonanceRadiusAU(aP, r.m, r.n, /*interior=*/true);
    if (aR < belt.aInnerAU || aR > belt.aOuterAU) continue;

    BeltResonanceFeature f{};
    f.m = r.m;
    f.n = r.n;
    f.aAU = aR;

    // Gaps get narrower at higher-order resonances.
    const double order = std::max(1, r.m - r.n);
    const double baseWidth = wB * rng.range(0.010, 0.028);
    f.halfWidthAU = std::clamp(baseWidth / std::sqrt(order), wB * 0.004, wB * 0.05);

    // Strength shaped by resonance order + a small deterministic jitter.
    const double s = std::clamp(r.base + rng.range(-0.08, 0.08), 0.08, 0.95);
    f.strength01 = s;

    belt.resonances.push_back(f);
  }

  // Sort by radius for nicer UI.
  std::sort(belt.resonances.begin(), belt.resonances.end(),
            [](const BeltResonanceFeature& a, const BeltResonanceFeature& b) { return a.aAU < b.aAU; });
}

static void addOuterResonantRidges(AsteroidBelt& belt,
                                  const sim::Planet& planet,
                                  core::SplitMix64& rng) {
  // A few exterior resonances that create enhanced populations ("Plutinos" etc).
  struct R { int m; int n; double base; };
  const R rs[] = {
      {3, 2, 0.30},
      {2, 1, 0.22},
      {5, 2, 0.16},
  };

  const double aP = std::max(0.05, planet.orbit.semiMajorAxisAU);
  const double wB = std::max(1e-6, belt.aOuterAU - belt.aInnerAU);

  for (const auto& r : rs) {
    const double aR = resonanceRadiusAU(aP, r.m, r.n, /*interior=*/false);
    if (aR < belt.aInnerAU || aR > belt.aOuterAU) continue;

    BeltResonanceFeature f{};
    f.m = r.m;
    f.n = r.n;
    f.aAU = aR;

    const double order = std::max(1, r.m - r.n);
    const double baseWidth = wB * rng.range(0.012, 0.040);
    f.halfWidthAU = std::clamp(baseWidth / std::sqrt(order), wB * 0.006, wB * 0.08);

    // Negative strength -> ridge/boost.
    f.strength01 = -std::clamp(r.base + rng.range(-0.06, 0.06), 0.06, 0.55);
    belt.resonances.push_back(f);
  }

  std::sort(belt.resonances.begin(), belt.resonances.end(),
            [](const BeltResonanceFeature& a, const BeltResonanceFeature& b) { return a.aAU < b.aAU; });
}

} // namespace

const char* asteroidBeltKindName(AsteroidBeltKind k) {
  switch (k) {
    case AsteroidBeltKind::MainBelt: return "Main Belt";
    case AsteroidBeltKind::OuterBelt: return "Outer Belt";
    case AsteroidBeltKind::TrojanSwarm: return "Trojan Swarm";
    case AsteroidBeltKind::DebrisDisk: return "Debris Disk";
    default: return "Belt";
  }
}

AsteroidBeltPlan generateAsteroidBelts(core::u64 universeSeed, const sim::StarSystem& system) {
  AsteroidBeltPlan plan;

  const core::u64 stream = core::seedFromText("proc.asteroid_belts.v1");
  const core::u64 seed = core::hashCombine(universeSeed, core::hashCombine(system.stub.seed, stream));
  core::SplitMix64 rng(seed);

  // Shared plane basis derived from the current planet orbits.
  math::Vec3d bx, by, bz;
  buildOrthonormalBasis(meanPlanetPlaneNormal(system), bx, by, bz);

  auto makeId = [&](int idx) -> core::u64 {
    core::u64 h = system.stub.id;
    h = core::hashCombine(h, stream);
    h = core::hashCombine(h, static_cast<core::u64>(idx));
    return h;
  };

  const int innerGiant = pickInnerGiant(system);
  const int outerMost = pickOuterMostPlanet(system);
  const int dominant = pickDominantPlanet(system);

  double aNonGiantMax = 0.0;
  double aOuterMax = 0.0;
  for (const auto& p : system.planets) {
    aOuterMax = std::max(aOuterMax, p.orbit.semiMajorAxisAU);
    if (p.type != sim::PlanetType::GasGiant) {
      aNonGiantMax = std::max(aNonGiantMax, p.orbit.semiMajorAxisAU);
    }
  }

  int beltIndex = 0;

  // ---------------------------------------------------------------------------
  // Main asteroid belt (between inner non-giants and the first giant).
  // ---------------------------------------------------------------------------
  if (innerGiant >= 0 && aNonGiantMax > 0.0) {
    const double aG = system.planets[(std::size_t)innerGiant].orbit.semiMajorAxisAU;

    AsteroidBelt b{};
    b.id = makeId(beltIndex++);
    b.kind = AsteroidBeltKind::MainBelt;

    // Ensure a visible gap.
    const double inner = aNonGiantMax * rng.range(1.04, 1.14);
    const double outer = aG * rng.range(0.80, 0.93);

    if (outer - inner > 0.35 && outer > inner * 1.12) {
      b.aInnerAU = std::max(0.12, inner);
      b.aOuterAU = std::max(b.aInnerAU + 0.10, outer);

      const double w = b.aOuterAU - b.aInnerAU;
      b.thicknessAU = std::clamp(w * rng.range(0.010, 0.030), 0.002, 0.12);

      b.basisX = bx;
      b.basisY = by;
      b.basisZ = bz;

      b.controllingPlanetIndex = innerGiant;

      // Subtle azimuthal modulation (belt clumps / spokes).
      b.mMode = rng.range(3, 9);
      b.mStrength01 = rng.range(0.05, 0.18);
      b.mPhaseRad = rng.range(0.0, kTwoPi);

      // Resonance gaps with the giant planet (Kirkwood-style).
      {
        core::SplitMix64 rrng(core::hashCombine(seed, core::hashCombine(b.id, core::seedFromText("gaps"))));
        addKirkwoodGaps(b, system.planets[(std::size_t)innerGiant], rrng);
      }

      plan.belts.push_back(std::move(b));
    }
  }

  // ---------------------------------------------------------------------------
  // Outer belt (Kuiper-like) beyond the outermost planet.
  // ---------------------------------------------------------------------------
  if (outerMost >= 0 && aOuterMax > 0.25) {
    const bool spawn = rng.chance(system.planets.size() >= 5 ? 0.85 : 0.65);
    if (spawn) {
      AsteroidBelt b{};
      b.id = makeId(beltIndex++);
      b.kind = AsteroidBeltKind::OuterBelt;

      const double aP = system.planets[(std::size_t)outerMost].orbit.semiMajorAxisAU;
      b.aInnerAU = std::max(0.20, aP * rng.range(1.18, 1.55));
      b.aOuterAU = std::max(b.aInnerAU + 0.25, b.aInnerAU * rng.range(1.55, 2.85));

      const double w = b.aOuterAU - b.aInnerAU;
      b.thicknessAU = std::clamp(w * rng.range(0.020, 0.060), 0.004, 0.40);

      b.basisX = bx;
      b.basisY = by;
      b.basisZ = bz;

      b.controllingPlanetIndex = outerMost;

      // Outer belts are a bit less "spokey".
      b.mMode = rng.range(0, 6);
      b.mStrength01 = (b.mMode > 0) ? rng.range(0.02, 0.10) : 0.0;
      b.mPhaseRad = rng.range(0.0, kTwoPi);

      // Add resonance ridges (enhanced populations) with the outermost planet.
      {
        core::SplitMix64 rrng(core::hashCombine(seed, core::hashCombine(b.id, core::seedFromText("ridges"))));
        addOuterResonantRidges(b, system.planets[(std::size_t)outerMost], rrng);
      }

      plan.belts.push_back(std::move(b));
    }
  }

  // ---------------------------------------------------------------------------
  // Trojan swarms near the dominant planet.
  // ---------------------------------------------------------------------------
  if (dominant >= 0 && rng.chance(0.70)) {
    const auto& p = system.planets[(std::size_t)dominant];
    const double aP = std::max(0.08, p.orbit.semiMajorAxisAU);

    // Avoid trojans when the planet is extremely close to the star.
    if (aP > 0.35) {
      AsteroidBelt b{};
      b.id = makeId(beltIndex++);
      b.kind = AsteroidBeltKind::TrojanSwarm;

      b.aInnerAU = aP * rng.range(0.86, 0.93);
      b.aOuterAU = aP * rng.range(1.07, 1.14);

      b.basisX = bx;
      b.basisY = by;
      b.basisZ = bz;

      b.controllingPlanetIndex = dominant;

      b.trojanCenterThetaRad = wrapAnglePi(p.orbit.meanAnomalyAtEpochRad);
      b.trojanWidthRad = rng.range(stellar::math::degToRad(12.0), stellar::math::degToRad(35.0));

      // Spread scales with orbital radius.
      b.trojanRadialSigmaAU = std::clamp(aP * rng.range(0.020, 0.070), 0.01, 0.80);
      b.trojanVerticalSigmaAU = std::clamp(aP * rng.range(0.006, 0.022), 0.002, 0.25);

      // Trojan swarms are not "spoked".
      b.mMode = 0;
      b.mStrength01 = 0.0;
      b.mPhaseRad = 0.0;

      // Optionally add a weak 1:1 resonance ridge at the planet orbit (just for UI).
      BeltResonanceFeature f{};
      f.m = 1;
      f.n = 1;
      f.aAU = aP;
      f.halfWidthAU = std::max(0.004, (b.aOuterAU - b.aInnerAU) * 0.20);
      f.strength01 = -rng.range(0.10, 0.22);
      b.resonances.push_back(f);

      plan.belts.push_back(std::move(b));
    }
  }

  // ---------------------------------------------------------------------------
  // Fallback: a thin debris disk if the system has few planets and no belts.
  // ---------------------------------------------------------------------------
  if (plan.belts.empty()) {
    AsteroidBelt b{};
    b.id = makeId(beltIndex++);
    b.kind = AsteroidBeltKind::DebrisDisk;

    const double a0 = std::max(0.15, aOuterMax > 0.0 ? aOuterMax * 0.85 : rng.range(0.6, 3.0));
    b.aInnerAU = std::max(0.10, a0 * rng.range(0.65, 0.92));
    b.aOuterAU = std::max(b.aInnerAU + 0.25, a0 * rng.range(1.25, 1.90));
    b.thicknessAU = std::clamp((b.aOuterAU - b.aInnerAU) * rng.range(0.010, 0.030), 0.002, 0.12);

    b.basisX = bx;
    b.basisY = by;
    b.basisZ = bz;

    b.controllingPlanetIndex = -1;
    b.mMode = rng.range(0, 8);
    b.mStrength01 = (b.mMode > 0) ? rng.range(0.03, 0.16) : 0.0;
    b.mPhaseRad = rng.range(0.0, kTwoPi);

    plan.belts.push_back(std::move(b));
  }

  return plan;
}

double asteroidBeltDensity01(const AsteroidBelt& belt, double aAU, double thetaRad) {
  const double a0 = std::max(0.0, belt.aInnerAU);
  const double a1 = std::max(a0 + 1e-9, belt.aOuterAU);

  // Normalize radial coordinate.
  const double t = (aAU - a0) / (a1 - a0);
  if (t <= -0.10 || t >= 1.10) return 0.0;

  // Soft edges (avoid hard cut-offs).
  const double edge = 0.055;
  double d = smoothstep(0.0, edge, t) * (1.0 - smoothstep(1.0 - edge, 1.0, t));

  // Gentle radial shape: main belts slightly peak toward the middle; outer belts are broader.
  const double mid = 1.0 - std::abs(t - 0.5) * 2.0;
  const double midPow = (belt.kind == AsteroidBeltKind::OuterBelt) ? 0.55 : 1.10;
  d *= std::pow(std::clamp(mid, 0.0, 1.0), midPow) * 0.55 + 0.45;

  // Resonance features (gaps + ridges).
  for (const auto& f : belt.resonances) {
    if (f.halfWidthAU <= 1e-12) continue;
    const double x = (aAU - f.aAU) / f.halfWidthAU;
    const double g = gauss01(x, 1.0);

    // strength01 > 0 -> dip, < 0 -> ridge.
    const double s = std::clamp(f.strength01, -0.95, 0.95);
    if (s >= 0.0) {
      d *= (1.0 - s * g);
    } else {
      d *= (1.0 + (-s) * g);
    }
  }

  // Trojan swarms: L4/L5 at +/- 60 degrees relative to trojanCenterThetaRad.
  if (belt.kind == AsteroidBeltKind::TrojanSwarm) {
    const double aCenter = 0.5 * (a0 + a1);
    const double ra = gauss01(aAU - aCenter, std::max(1e-6, belt.trojanRadialSigmaAU));

    const double w = std::max(1e-6, belt.trojanWidthRad);
    const double th0 = belt.trojanCenterThetaRad;
    const double thL4 = th0 + stellar::math::degToRad(60.0);
    const double thL5 = th0 - stellar::math::degToRad(60.0);

    const double d4 = gauss01(wrapAnglePi(thetaRad - thL4), w);
    const double d5 = gauss01(wrapAnglePi(thetaRad - thL5), w);

    // Two lobes; keep background low so it reads as two swarms.
    d *= (0.05 + 0.95 * std::max(d4, d5)) * ra;
  }

  // Azimuthal modulation (subtle clumping).
  if (belt.mMode > 0 && belt.mStrength01 > 0.0) {
    const double amp = std::clamp(belt.mStrength01, 0.0, 0.95);
    const double m = static_cast<double>(belt.mMode);
    const double wave = 0.5 + 0.5 * std::cos(m * thetaRad + belt.mPhaseRad);
    d *= (1.0 - amp) + amp * (0.25 + 0.75 * wave);
  }

  return clamp01(d);
}

std::vector<AsteroidBeltPoint> sampleAsteroidBeltPoints(core::u64 universeSeed,
                                                        const sim::StarSystem& system,
                                                        const AsteroidBelt& belt,
                                                        int pointCount,
                                                        int candidatesPerPoint) {
  pointCount = std::clamp(pointCount, 0, 120000);
  candidatesPerPoint = std::clamp(candidatesPerPoint, 1, 96);

  std::vector<AsteroidBeltPoint> out;
  out.reserve(static_cast<std::size_t>(pointCount));

  if (pointCount <= 0) return out;

  const core::u64 stream = core::seedFromText("proc.asteroid_belt_points.v1");
  const core::u64 seed = core::hashCombine(universeSeed,
                                           core::hashCombine(system.stub.seed,
                                                             core::hashCombine(stream, belt.id)));
  core::SplitMix64 rng(seed);

  // A loose target spacing (in AU) used for the clearance term.
  const double w = std::max(1e-6, belt.aOuterAU - belt.aInnerAU);
  const double target = std::clamp(w / std::sqrt((double)pointCount), w * 0.002, w * 0.08);
  const double target2 = target * target;

  auto candidatePos = [&](std::uint32_t index1Based) -> std::pair<math::Vec3d, double> {
    // Use Halton for deterministic low-discrepancy sampling.
    const core::Halton3 h = core::halton3(index1Based);
    const double u0 = h.x;
    const double u1 = h.y;
    const double u2 = h.z;

    // Slightly bias toward mid-radii for aesthetics.
    const double t = std::pow(std::clamp(u0, 0.0, 1.0), 0.88);
    const double a = lerp(belt.aInnerAU, belt.aOuterAU, t);

    const double theta = u1 * kTwoPi;

    // Vertical distribution: gaussian-ish via a smooth remap.
    const double th = std::max(1e-9, (belt.kind == AsteroidBeltKind::TrojanSwarm)
                                     ? belt.trojanVerticalSigmaAU
                                     : belt.thicknessAU);
    const double y = (u2 * 2.0 - 1.0) * th;

    const double x = a * std::cos(theta);
    const double z = a * std::sin(theta);

    const math::Vec3d pos = belt.basisX * x + belt.basisZ * z + belt.basisY * y;
    const double dens2D = asteroidBeltDensity01(belt, a, theta);

    // Apply vertical falloff (gaussian).
    const double dv = gauss01(y, th);
    return {pos, clamp01(dens2D * dv)};
  };

  // Mitchell-style best-candidate sampling with an importance weight.
  //
  // Score blends:
  //  - density (importance)
  //  - clearance (blue-noise-ish)
  //  - a tiny deterministic jitter (breaks ties)
  for (int i = 0; i < pointCount; ++i) {
    double bestScore = -1.0e30;
    AsteroidBeltPoint best{};

    for (int j = 0; j < candidatesPerPoint; ++j) {
      const std::uint32_t idx = static_cast<std::uint32_t>(i * candidatesPerPoint + j + 1);
      auto [pos, dens] = candidatePos(idx);

      // Clearance term.
      double minD2 = std::numeric_limits<double>::infinity();
      for (const auto& p : out) {
        const math::Vec3d d = pos - p.posAU;
        const double dd = d.lengthSq();
        if (dd < minD2) minD2 = dd;
      }

      const double clear = (out.empty() || !std::isfinite(minD2)) ? 1.0 : std::sqrt(minD2 / target2);
      const double clear01 = std::clamp(clear, 0.0, 3.0) / 3.0;

      // Combine terms.
      const double jitter = rng.nextDouble() * 1e-3;
      const double score = (dens * 0.72 + clear01 * 0.28) + jitter;

      if (score > bestScore) {
        bestScore = score;
        best.posAU = pos;
        best.density01 = dens;
      }
    }

    out.push_back(best);
  }

  return out;
}

} // namespace stellar::proc
