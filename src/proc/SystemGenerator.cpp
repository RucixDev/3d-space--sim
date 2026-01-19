#include "stellar/proc/SystemGenerator.h"

#include "stellar/core/Random.h"
#include "stellar/econ/Economy.h"
#include "stellar/math/Math.h"
#include "stellar/proc/NameGenerator.h"
#include "stellar/proc/TradeEconomy.h"
#include "stellar/proc/TradeProfile.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

// NOTE:
// This file intentionally stays deterministic from (SystemStub.seed + galaxy factions).
// Be careful when introducing new RNG draws; it will perturb downstream procedural
// content and may require updating signature-based regression tests.

namespace stellar::proc {

static sim::Star makeStar(sim::StarClass cls, core::SplitMix64& rng) {
  sim::Star s{};
  s.cls = cls;

  auto rr = [&](double a, double b) { return rng.range(a, b); };

  switch (cls) {
    case sim::StarClass::O:
      s.massSol = rr(16.0, 60.0);
      s.radiusSol = rr(6.0, 15.0);
      s.luminositySol = rr(30000.0, 500000.0);
      s.temperatureK = rr(30000.0, 50000.0);
      break;
    case sim::StarClass::B:
      s.massSol = rr(2.1, 16.0);
      s.radiusSol = rr(2.0, 6.0);
      s.luminositySol = rr(25.0, 30000.0);
      s.temperatureK = rr(10000.0, 30000.0);
      break;
    case sim::StarClass::A:
      s.massSol = rr(1.4, 2.1);
      s.radiusSol = rr(1.4, 2.5);
      s.luminositySol = rr(5.0, 25.0);
      s.temperatureK = rr(7500.0, 10000.0);
      break;
    case sim::StarClass::F:
      s.massSol = rr(1.04, 1.4);
      s.radiusSol = rr(1.15, 1.6);
      s.luminositySol = rr(1.5, 5.0);
      s.temperatureK = rr(6000.0, 7500.0);
      break;
    case sim::StarClass::G:
      s.massSol = rr(0.8, 1.04);
      s.radiusSol = rr(0.9, 1.2);
      s.luminositySol = rr(0.6, 1.6);
      s.temperatureK = rr(5200.0, 6000.0);
      break;
    case sim::StarClass::K:
      s.massSol = rr(0.45, 0.8);
      s.radiusSol = rr(0.7, 0.95);
      s.luminositySol = rr(0.08, 0.6);
      s.temperatureK = rr(3700.0, 5200.0);
      break;
    case sim::StarClass::M:
    default:
      s.massSol = rr(0.08, 0.45);
      s.radiusSol = rr(0.1, 0.7);
      s.luminositySol = rr(0.0001, 0.08);
      s.temperatureK = rr(2400.0, 3700.0);
      break;
  }
  return s;
}

static sim::PlanetType pickPlanetType(double aAU, core::SplitMix64& rng) {
  // crude zone heuristic + randomness
  const double r = rng.nextDouble();
  if (aAU < 0.6) {
    return (r < 0.7) ? sim::PlanetType::Rocky : sim::PlanetType::Desert;
  }
  if (aAU < 2.0) {
    if (r < 0.5) return sim::PlanetType::Rocky;
    if (r < 0.8) return sim::PlanetType::Ocean;
    return sim::PlanetType::Ice;
  }
  if (aAU < 6.0) {
    return (r < 0.3) ? sim::PlanetType::GasGiant : sim::PlanetType::Ice;
  }
  return (r < 0.6) ? sim::PlanetType::GasGiant : sim::PlanetType::Ice;
}

static void setPlanetMassRadius(sim::Planet& p, core::SplitMix64& rng) {
  auto rr = [&](double a, double b) { return rng.range(a, b); };
  switch (p.type) {
    case sim::PlanetType::GasGiant:
      p.radiusEarth = rr(3.0, 11.0);
      p.massEarth = rr(20.0, 320.0);
      break;
    case sim::PlanetType::Ocean:
      p.radiusEarth = rr(0.8, 2.2);
      p.massEarth = rr(0.6, 8.0);
      break;
    case sim::PlanetType::Ice:
      p.radiusEarth = rr(0.4, 1.6);
      p.massEarth = rr(0.2, 5.0);
      break;
    case sim::PlanetType::Desert:
      p.radiusEarth = rr(0.6, 2.0);
      p.massEarth = rr(0.4, 7.0);
      break;
    case sim::PlanetType::Rocky:
    default:
      p.radiusEarth = rr(0.3, 1.8);
      p.massEarth = rr(0.1, 6.0);
      break;
  }
}

// -------------------------
// Moons (procedural only)
// -------------------------
//
// IMPORTANT: These helpers intentionally use their own RNG stream(s) seeded
// from SystemStub.seed, so adding moons does NOT perturb downstream draws for
// stations, names, or economy state.

static constexpr double kEarthMassSol = 3.0034146856628466e-6;    // M_earth / M_sun
static constexpr double kEarthRadiusAU = 4.258750455597227e-5;    // R_earth in AU

static double planetMassSol(const sim::Planet& p) {
  return std::max(0.0, p.massEarth) * kEarthMassSol;
}

static double planetRadiusAU(const sim::Planet& p) {
  return std::max(0.0, p.radiusEarth) * kEarthRadiusAU;
}

// Hill radius of a planet around its star (AU).
static double hillRadiusAU(const sim::Planet& p, const sim::Star& star) {
  const double mStar = std::max(0.08, star.massSol);
  const double mPlanet = planetMassSol(p);
  if (mPlanet <= 0.0 || p.orbit.semiMajorAxisAU <= 0.0) return 0.0;
  return p.orbit.semiMajorAxisAU * std::cbrt(mPlanet / (3.0 * mStar));
}

static sim::PlanetType pickMoonType(const sim::Planet& host,
                                    double hostSemiMajorAU,
                                    const sim::Star& star,
                                    core::SplitMix64& rng) {
  // Rough temperature proxy: use star luminosity to compute an Earth-normalized
  // "distance" scale. Colder systems bias moons toward icy.
  const double L = std::max(1e-6, star.luminositySol);
  const double aEq = std::sqrt(L);
  const double aRel = (aEq > 1e-9) ? (hostSemiMajorAU / aEq) : hostSemiMajorAU;

  const double r = rng.nextDouble();

  // Moons should never be gas giants.
  if (host.type == sim::PlanetType::GasGiant) {
    // Beyond ~2 AU (scaled by luminosity), most large moons are icy.
    if (aRel > 2.0) {
      if (r < 0.75) return sim::PlanetType::Ice;
      if (r < 0.93) return sim::PlanetType::Rocky;
      return sim::PlanetType::Ocean;
    }
    // Inner systems: mostly rocky/volcanic (we approximate with Rocky/Desert).
    if (r < 0.65) return sim::PlanetType::Rocky;
    if (r < 0.90) return sim::PlanetType::Desert;
    return sim::PlanetType::Ocean;
  }

  // Terrestrial planets: rare moons, typically rocky/icy depending on distance.
  if (aRel > 1.5) {
    return (r < 0.75) ? sim::PlanetType::Ice : sim::PlanetType::Rocky;
  }
  return (r < 0.85) ? sim::PlanetType::Rocky : sim::PlanetType::Desert;
}

static void setMoonMassRadius(sim::Moon& m, core::SplitMix64& rng) {
  auto rr = [&](double a, double b) { return rng.range(a, b); };

  // Radius ranges (Earth radii). Keep them modest so moons read visually as moons.
  switch (m.type) {
    case sim::PlanetType::Ocean:
      m.radiusEarth = rr(0.07, 0.40);
      break;
    case sim::PlanetType::Ice:
      m.radiusEarth = rr(0.05, 0.35);
      break;
    case sim::PlanetType::Desert:
      m.radiusEarth = rr(0.05, 0.38);
      break;
    case sim::PlanetType::Rocky:
    default:
      m.radiusEarth = rr(0.05, 0.42);
      break;
  }

  // Simple density proxy relative to Earth.
  double density = 1.0;
  switch (m.type) {
    case sim::PlanetType::Ice: density = rr(0.55, 0.75); break;
    case sim::PlanetType::Ocean: density = rr(0.65, 0.90); break;
    case sim::PlanetType::Desert: density = rr(0.85, 1.05); break;
    case sim::PlanetType::Rocky:
    default: density = rr(0.90, 1.20); break;
  }

  // Mass ~ density * radius^3 (Earth units).
  const double r = std::max(0.01, m.radiusEarth);
  m.massEarth = std::clamp(density * (r * r * r), 1e-5, 0.25);
}

static std::string moonSuffix(std::size_t idx) {
  // a, b, c... (wrap after z)
  const char c = static_cast<char>('a' + (idx % 26));
  return std::string(1, c);
}

static void appendMoonsForPlanet(std::vector<sim::Moon>& out,
                                 const sim::SystemStub& stub,
                                 const sim::Star& star,
                                 const sim::Planet& host,
                                 std::size_t hostIndex) {
  // Separate stream per host planet.
  //
  // NOTE: This stream is intentionally versioned so future moon generator tweaks
  // do not perturb other procedural content (stations/economy/etc.).
  const core::u64 stream = core::seedFromText("proc.moons.v2");
  core::SplitMix64 rng(core::hashCombine(stub.seed, core::hashCombine(stream, (core::u64)hostIndex)));

  // Quick gate: tiny bodies rarely hold moons.
  if (host.massEarth < 0.12 || host.radiusEarth < 0.25) return;

  const std::size_t baseCount = out.size();
  const std::size_t hardCap = 32; // total moons per system cap (safety)
  if (baseCount >= hardCap) return;

  // Stable satellite band.
  const double hill = hillRadiusAU(host, star);
  if (hill <= 0.0) return;

  // Prograde stable region is ~0.5 R_H; we stay comfortably inside.
  const double aMaxAU = hill * 0.45;

  // --- Physics-inspired inner cutoff (Roche + size buffer) -----------------
  auto densityRelPlanet = [&](const sim::Planet& p) {
    const double r = std::max(0.05, p.radiusEarth);
    const double m = std::max(1.0e-6, p.massEarth);
    return std::clamp(m / (r * r * r), 0.05, 10.0);
  };

  auto densityRelMoon = [&](sim::PlanetType t) {
    switch (t) {
      case sim::PlanetType::Ice:   return 0.65;
      case sim::PlanetType::Ocean: return 0.80;
      case sim::PlanetType::Desert:return 0.95;
      case sim::PlanetType::Rocky:
      default:                     return 1.05;
    }
  };

  auto rocheLimitAU = [&](sim::PlanetType moonType) {
    // Fluid Roche limit approximation: 2.44 R_p * (rho_p / rho_s)^(1/3)
    const double rhoP = densityRelPlanet(host);
    const double rhoS = std::max(0.25, densityRelMoon(moonType));
    const double k = 2.44 * std::cbrt(rhoP / rhoS);
    return planetRadiusAU(host) * k;
  };

  // Inner cutoff: stay outside Roche and also a conservative multiple of
  // the host radius (covers atmosphere + “close orbit” clutter).
  // NOTE: tests assume aAU > 4*R_host.
  const double hostRAU = planetRadiusAU(host);
  const double aHardMinAU = std::max(hostRAU * 4.0, 1.0e-6);

  // We don't know moon type yet, so estimate with a typical density.
  const double aRocheEst = std::max(rocheLimitAU(sim::PlanetType::Rocky), rocheLimitAU(sim::PlanetType::Ice));
  const double aMinAU = std::max(aHardMinAU, aRocheEst * 1.15);

  if (!(aMaxAU > aMinAU * 1.25)) return;

  // --- Target moon counts --------------------------------------------------
  int targetRegular = 0;
  switch (host.type) {
    case sim::PlanetType::GasGiant:
      targetRegular = rng.range(3, 9);
      break;
    case sim::PlanetType::Ice:
      targetRegular = rng.chance(0.70) ? rng.range(1, 4) : rng.range(0, 2);
      break;
    case sim::PlanetType::Ocean:
      targetRegular = rng.chance(0.28) ? rng.range(1, 2) : 0;
      break;
    case sim::PlanetType::Desert:
    case sim::PlanetType::Rocky:
    default:
      targetRegular = (host.massEarth > 0.8 && rng.chance(0.34)) ? rng.range(1, 2) : 0;
      break;
  }

  // If the host is a gas giant and has room for satellites, guarantee at least one.
  if (host.type == sim::PlanetType::GasGiant) targetRegular = std::max(1, targetRegular);

  // Room check against hard cap.
  targetRegular = std::min<int>(targetRegular, (int)(hardCap - baseCount));
  if (targetRegular <= 0) return;

  // --- Utility: mutual Hill separation -------------------------------------
  auto mutualHillAU = [&](double a1, double a2, double m1Earth, double m2Earth) {
    const double M = std::max(1.0e-6, host.massEarth);
    const double mu = std::max(1.0e-12, (m1Earth + m2Earth) / (3.0 * M));
    const double aBar = 0.5 * (a1 + a2);
    return std::cbrt(mu) * std::max(1.0e-12, aBar);
  };

  struct Placed {
    double aAU{0.0};
    double mEarth{0.0};
    std::size_t idx{0};
  };

  auto minSepInMutualHill = [&](double aAU, double mEarth, const std::vector<Placed>& placed) {
    double best = 1.0e9;
    for (const auto& p : placed) {
      const double rH = mutualHillAU(aAU, p.aAU, mEarth, p.mEarth);
      const double sep = std::abs(aAU - p.aAU) / std::max(1.0e-12, rH);
      best = std::min(best, sep);
    }
    return best;
  };

  auto sampleLogUniform = [&](double lo, double hi, double biasPow) {
    lo = std::max(lo, 1.0e-12);
    hi = std::max(hi, lo * 1.000001);
    const double l0 = std::log(lo);
    const double l1 = std::log(hi);
    double t = rng.nextDouble();
    // biasPow > 1 => bias toward inner orbits (lo)
    if (biasPow > 1.0) t = std::pow(t, biasPow);
    return std::exp(math::lerp(l0, l1, t));
  };

  auto pickBestA = [&](double lo, double hi, double mEarth, const std::vector<Placed>& placed,
                       int candidates, double biasPow, double sepMin, double& outA) {
    double bestScore = -1.0;
    double bestA = lo;
    for (int c = 0; c < candidates; ++c) {
      const double a = sampleLogUniform(lo, hi, biasPow);
      const double score = placed.empty() ? 1.0e9 : minSepInMutualHill(a, mEarth, placed);
      if (score > bestScore) {
        bestScore = score;
        bestA = a;
      }
    }
    if (!placed.empty() && bestScore < sepMin) return false;
    outA = bestA;
    return true;
  };

  // --- Generate regular moons (low-i, mostly prograde) ----------------------
  std::vector<sim::Moon> moons;
  moons.reserve((std::size_t)targetRegular + 4);

  for (int i = 0; i < targetRegular; ++i) {
    if (baseCount + moons.size() >= hardCap) break;

    sim::Moon m{};
    m.parentPlanetIndex = static_cast<core::u32>(hostIndex);
    m.id = core::hashCombine(stub.id, core::hashCombine((core::u64)hostIndex, (core::u64)i));
    m.name = host.name + " " + moonSuffix((std::size_t)i);

    m.type = pickMoonType(host, host.orbit.semiMajorAxisAU, star, rng);
    setMoonMassRadius(m, rng);

    moons.push_back(std::move(m));
  }

  if (moons.empty()) return;

  // Place larger moons first (more restrictive) for higher success rate.
  std::vector<std::size_t> order(moons.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    if (moons[a].massEarth != moons[b].massEarth) return moons[a].massEarth > moons[b].massEarth;
    return moons[a].id < moons[b].id;
  });

  // Regular satellite disk occupies a small inner band (gas giants) to create a
  // visually familiar "Galilean-style" structure, while leaving room for
  // irregular moons further out.
  double aRegMaxAU = aMaxAU;
  if (host.type == sim::PlanetType::GasGiant) {
    aRegMaxAU = std::min(aMaxAU, std::max(aMinAU * 4.0, hill * rng.range(0.06, 0.11)));
  }

  const int candidates = (host.type == sim::PlanetType::GasGiant) ? 96 : 64;
  const double biasPow = (host.type == sim::PlanetType::GasGiant) ? 2.2 : 1.4;
  const double sepMin = (host.type == sim::PlanetType::GasGiant) ? 4.5 : 5.0;

  std::vector<Placed> placed;
  placed.reserve(moons.size());

  for (std::size_t oi = 0; oi < order.size(); ++oi) {
    const std::size_t mi = order[oi];
    double aAU = 0.0;

    if (!pickBestA(aMinAU, aRegMaxAU, moons[mi].massEarth, placed, candidates, biasPow, sepMin, aAU)) {
      // If we can't place the smallest moons, just stop.
      continue;
    }

    moons[mi].orbit.semiMajorAxisAU = aAU;
    placed.push_back(Placed{aAU, moons[mi].massEarth, mi});
  }

  // Remove any moons we failed to place.
  moons.erase(std::remove_if(moons.begin(), moons.end(), [](const sim::Moon& m) {
               return !(m.orbit.semiMajorAxisAU > 0.0);
             }),
             moons.end());

  if (moons.empty()) return;

  // Resonance snapping pass (gentle): if two neighbors are already close to a
  // simple rational period ratio, nudge the outer one toward that resonance.
  // This creates more interesting, less "random" satellite systems.
  struct Res {
    int num;
    int den;
    double r; // period ratio
  };
  const std::array<Res, 6> kRes = {Res{4,3,4.0/3.0}, Res{3,2,1.5}, Res{5,3,5.0/3.0},
                                  Res{2,1,2.0}, Res{5,2,2.5}, Res{3,1,3.0}};

  auto nearestResonance = [&](double periodRatio) {
    Res best = kRes[0];
    double bestErr = 1.0e9;
    for (const auto& rr : kRes) {
      const double err = std::abs(periodRatio - rr.r) / rr.r;
      if (err < bestErr) {
        bestErr = err;
        best = rr;
      }
    }
    return std::pair<Res, double>(best, bestErr);
  };

  std::sort(moons.begin(), moons.end(), [](const sim::Moon& a, const sim::Moon& b) {
    return a.orbit.semiMajorAxisAU < b.orbit.semiMajorAxisAU;
  });

  const double snapTol = (host.type == sim::PlanetType::GasGiant) ? 0.030 : 0.022;
  const double snapProb = (host.type == sim::PlanetType::GasGiant) ? 0.78 : 0.40;

  for (std::size_t i = 1; i < moons.size(); ++i) {
    if (!rng.chance(snapProb)) continue;

    const double a0 = moons[i - 1].orbit.semiMajorAxisAU;
    const double a1 = moons[i].orbit.semiMajorAxisAU;
    if (!(a1 > a0 * 1.05)) continue;

    const double pr = std::pow(a1 / a0, 1.5);
    const auto [res, relErr] = nearestResonance(pr);
    if (relErr > snapTol) continue;

    const double aRes = a0 * std::pow(res.r, 2.0 / 3.0);
    const double jitter = 1.0 + rng.range(-0.0025, 0.0025);
    const double aNew = aRes * jitter;

    if (!(aNew > aMinAU && aNew < aRegMaxAU)) continue;

    // Check mutual Hill spacing vs immediate neighbors only.
    bool ok = true;
    {
      const double sepPrev = std::abs(aNew - a0) / std::max(1.0e-12, mutualHillAU(aNew, a0, moons[i].massEarth, moons[i - 1].massEarth));
      if (sepPrev < sepMin) ok = false;
    }
    if (ok && i + 1 < moons.size()) {
      const double a2 = moons[i + 1].orbit.semiMajorAxisAU;
      const double sepNext = std::abs(a2 - aNew) / std::max(1.0e-12, mutualHillAU(a2, aNew, moons[i + 1].massEarth, moons[i].massEarth));
      if (sepNext < sepMin) ok = false;
    }

    if (ok) moons[i].orbit.semiMajorAxisAU = aNew;
  }

  // Fill orbital elements for regular moons (low eccentricity, near-host plane).
  for (auto& m : moons) {
    const double aAU = m.orbit.semiMajorAxisAU;
    const double rel = std::clamp((aAU - aMinAU) / std::max(1.0e-9, (aRegMaxAU - aMinAU)), 0.0, 1.0);

    // Eccentricities: close-in moons tend to be more circular.
    const double eMax = math::lerp(0.015, 0.060, rel);
    m.orbit.eccentricity = rng.range(0.0, eMax);

    // Small mutual inclination offsets (degrees).
    const double dInc = math::degToRad(rng.range(-1.2, 1.2));
    m.orbit.inclinationRad = host.orbit.inclinationRad + dInc;
    m.orbit.ascendingNodeRad = host.orbit.ascendingNodeRad + rng.range(-0.25, 0.25);
    m.orbit.argPeriapsisRad = rng.range(0.0, math::kPi * 2.0);
    m.orbit.meanAnomalyAtEpochRad = rng.range(0.0, math::kPi * 2.0);
    m.orbit.epochDays = 0.0;

    // Period: Kepler (a in AU, M in solar masses -> P in years).
    const double mHostSol = std::max(1e-12, planetMassSol(host));
    const double years = std::sqrt((aAU * aAU * aAU) / mHostSol);
    m.orbit.periodDays = years * 365.25;
  }

  // --- Optional irregular moons for gas giants --------------------------------
  // These are high-inclination, more eccentric satellites near the edge of the
  // Hill stability band. They add flavor and variety without needing any new
  // simulation systems.
  if (host.type == sim::PlanetType::GasGiant && (baseCount + moons.size()) < hardCap) {
    const int maxIr = std::min<int>(3, (int)(hardCap - (baseCount + moons.size())));
    const int irregularTarget = rng.chance(0.55) ? rng.range(0, maxIr) : 0;

    if (irregularTarget > 0) {
      // Define an outer band separated from the regular disk.
      const double aIrMinAU = std::max(aRegMaxAU * rng.range(1.35, 2.10), hill * rng.range(0.18, 0.28));
      const double aIrMaxAU = aMaxAU;

      if (aIrMaxAU > aIrMinAU * 1.10) {
        // Build a combined placed-list for spacing checks.
        std::vector<Placed> placedAll;
        placedAll.reserve(moons.size() + (std::size_t)irregularTarget);
        for (std::size_t i = 0; i < moons.size(); ++i) {
          placedAll.push_back(Placed{moons[i].orbit.semiMajorAxisAU, moons[i].massEarth, i});
        }

        const int irCandidates = 72;
        const double irBiasPow = 0.75; // bias toward outer edge
        const double irSepMin = 3.5;

        for (int j = 0; j < irregularTarget; ++j) {
          const int idx = targetRegular + j;

          sim::Moon m{};
          m.parentPlanetIndex = static_cast<core::u32>(hostIndex);
          m.id = core::hashCombine(stub.id, core::hashCombine((core::u64)hostIndex, (core::u64)idx));
          m.name = host.name + " " + moonSuffix((std::size_t)idx);

          // Irregulars are usually small; bias toward icy/rocky.
          const double r = rng.nextDouble();
          m.type = (r < 0.70) ? sim::PlanetType::Ice : sim::PlanetType::Rocky;
          setMoonMassRadius(m, rng);
          // Make them smaller on average.
          m.radiusEarth *= rng.range(0.35, 0.75);
          m.massEarth = std::clamp(m.massEarth * rng.range(0.15, 0.55), 1e-6, 0.02);

          double aAU = 0.0;
          if (!pickBestA(aIrMinAU, aIrMaxAU, m.massEarth, placedAll, irCandidates, irBiasPow, irSepMin, aAU)) {
            continue;
          }

          // Eccentric, high inclination, often retrograde.
          double e = rng.range(0.08, 0.35);
          // Avoid plunging into the regular satellite disk.
          const double q = aAU * (1.0 - e);
          if (q < aRegMaxAU * 1.10) {
            e = std::max(0.05, 1.0 - (aRegMaxAU * 1.10) / std::max(1.0e-9, aAU));
            e = std::clamp(e, 0.05, 0.35);
          }

          const bool retro = rng.chance(0.62);
          const double incDeg = retro ? rng.range(110.0, 170.0) : rng.range(30.0, 85.0);

          m.orbit.semiMajorAxisAU = aAU;
          m.orbit.eccentricity = e;
          m.orbit.inclinationRad = std::clamp(host.orbit.inclinationRad + math::degToRad(incDeg), 0.0, math::kPi);
          m.orbit.ascendingNodeRad = rng.range(0.0, math::kPi * 2.0);
          m.orbit.argPeriapsisRad = rng.range(0.0, math::kPi * 2.0);
          m.orbit.meanAnomalyAtEpochRad = rng.range(0.0, math::kPi * 2.0);
          m.orbit.epochDays = 0.0;

          const double mHostSol = std::max(1e-12, planetMassSol(host));
          const double years = std::sqrt((aAU * aAU * aAU) / mHostSol);
          m.orbit.periodDays = years * 365.25;

          // Accept.
          placedAll.push_back(Placed{aAU, m.massEarth, moons.size()});
          moons.push_back(std::move(m));
        }
      }
    }
  }

  // Append to output.
  for (auto& m : moons) {
    if (out.size() >= hardCap) break;
    out.push_back(std::move(m));
  }

  // Keep moons for the host sorted by semi-major axis to make lab/debug views stable.
  if (out.size() > baseCount + 1) {
    std::sort(out.begin() + static_cast<std::ptrdiff_t>(baseCount), out.end(),
              [](const sim::Moon& a, const sim::Moon& b) {
                if (a.parentPlanetIndex != b.parentPlanetIndex) return a.parentPlanetIndex < b.parentPlanetIndex;
                return a.orbit.semiMajorAxisAU < b.orbit.semiMajorAxisAU;
              });
  }
}

static double clamp01(double x) {
  return std::clamp(x, 0.0, 1.0);
}

// Convert a TradeProfile into a single -1..+1 axis for quick station-type biasing.
//
// Negative => agricultural leaning, Positive => industrial leaning.
static double localIndustryBias(const TradeProfile& p) {
  const double x = clamp01(p.industry) - clamp01(p.agriculture);
  return std::clamp(x * 1.25, -1.0, 1.0);
}

static std::array<double, static_cast<std::size_t>(econ::StationType::Count)>
stationTypeWeights(const TradeProfile& p, double bias) {
  using econ::StationType;
  std::array<double, static_cast<std::size_t>(StationType::Count)> w{};

  const double pop = clamp01(p.population);
  const double hub = clamp01(p.hub);
  const double res = clamp01(p.resources);
  const double ag  = clamp01(p.agriculture);
  const double ind = clamp01(p.industry);
  const double tech = clamp01(p.technology);
  const double wealth = clamp01(p.wealth);
  const double law = clamp01(p.lawlessness);

  // A gameplay-first model:
  // - population/hub => more hubs + larger secondary ecosystem
  // - resources/agriculture => upstream producers
  // - industry/tech => downstream producers
  // - lawlessness => more outposty, less mega-structures
  w[(std::size_t)StationType::Outpost]      = 0.08 + 0.22 * (1.0 - pop) + 0.12 * law;
  w[(std::size_t)StationType::Mining]       = 0.08 + 1.15 * res * (0.65 + 0.35 * (1.0 - ag));
  w[(std::size_t)StationType::Agricultural] = 0.08 + 1.15 * ag  * (0.65 + 0.35 * (1.0 - ind));
  w[(std::size_t)StationType::Refinery]     = 0.06 + 0.95 * res * ind;
  w[(std::size_t)StationType::Industrial]   = 0.08 + 1.30 * ind * (0.40 + 0.60 * pop);
  w[(std::size_t)StationType::Research]     = 0.05 + 1.10 * tech * (0.30 + 0.70 * pop);
  w[(std::size_t)StationType::TradeHub]     = 0.05 + 1.80 * hub * (0.30 + 0.70 * pop) + 0.25 * wealth;
  w[(std::size_t)StationType::Shipyard]     = 0.03 + 0.85 * ind * hub * (0.40 + 0.60 * wealth);

  // Apply the industry-vs-agri axis.
  w[(std::size_t)StationType::Agricultural] *= (bias < 0.0) ? (1.0 + (-bias) * 0.70) : (1.0 - bias * 0.10);
  w[(std::size_t)StationType::Industrial]   *= (bias > 0.0) ? (1.0 + ( bias) * 0.70) : (1.0 + (-bias) * 0.10);

  // Suppress the largest installations in anarchic space.
  const double bigMul = (1.0 - 0.45 * law);
  w[(std::size_t)StationType::TradeHub]   *= bigMul;
  w[(std::size_t)StationType::Shipyard]   *= bigMul;
  w[(std::size_t)StationType::Research]   *= (1.0 - 0.25 * law);
  w[(std::size_t)StationType::Industrial] *= (1.0 - 0.25 * law);

  // Ensure non-zero weights.
  for (double& x : w) {
    if (!std::isfinite(x) || x < 0.0) x = 0.0;
    x = std::max(x, 1e-4);
  }
  return w;
}

static econ::StationType pickStationType(core::SplitMix64& rng, const TradeProfile& p, double bias) {
  using econ::StationType;

  auto w = stationTypeWeights(p, bias);

  double sum = 0.0;
  for (double x : w) sum += std::max(0.0, x);
  if (sum <= 1e-12) return StationType::Outpost;

  double r = rng.nextDouble() * sum;
  for (std::size_t i = 0; i < w.size(); ++i) {
    r -= std::max(0.0, w[i]);
    if (r <= 0.0) return static_cast<StationType>(i);
  }
  return StationType::Outpost;
}

static econ::StationType primaryStationType(const TradeProfile& p, double bias) {
  using econ::StationType;
  auto w = stationTypeWeights(p, bias);

  std::size_t bestI = 0;
  double bestW = -1e30;
  for (std::size_t i = 0; i < w.size(); ++i) {
    const double ww = w[i];
    if (ww > bestW + 1e-12) {
      bestW = ww;
      bestI = i;
    }
  }
  return static_cast<StationType>(bestI);
}

static const sim::Faction* findFaction(core::u32 id, const std::vector<sim::Faction>& factions) {
  for (const auto& f : factions) if (f.id == id) return &f;
  return nullptr;
}

static double stationControlWeight(econ::StationType t) {
  // Mirrors the lightweight weighting in sim::computeSystemControl().
  // Keep this local to proc so we don't add a proc->sim dependency.
  using econ::StationType;
  switch (t) {
    case StationType::TradeHub: return 1.35;
    case StationType::Shipyard: return 1.30;
    case StationType::Research: return 1.20;
    case StationType::Industrial: return 1.15;
    case StationType::Refinery: return 1.10;
    case StationType::Mining: return 1.05;
    case StationType::Agricultural: return 1.05;
    case StationType::Outpost: return 1.0;
    default: return 1.0;
  }
}

static double combinedIndustryBias(double factionBias, const TradeProfile& profile) {
  // Blend faction-wide ideology with local macro conditions.
  // Keep in [-1, +1] where negative is agri-leaning and positive is industry-leaning.
  const double local = localIndustryBias(profile);
  return std::clamp(factionBias + 0.65 * local, -1.0, 1.0);
}

static void applyFactionToStation(sim::Station& st,
                                  core::u32 factionId,
                                  const std::vector<sim::Faction>& factions,
                                  const TradeProfile& profile,
                                  core::u64 stationSeed) {
  st.factionId = factionId;
  const sim::Faction* fac = findFaction(factionId, factions);
  const double baseFee = fac ? fac->taxRate : 0.02;
  const double facBias = fac ? fac->industryBias : 0.0;
  const double bias = combinedIndustryBias(facBias, profile);

  st.feeRate = tuneStationFeeRateForTradeProfile(stationSeed, baseFee, st.type, profile);

  const econ::StationEconomyModel baseM = econ::makeEconomyModel(st.type, bias);
  st.economyModel = tuneEconomyModelForTradeProfile(stationSeed, baseM, profile);
}

sim::StarSystem generateSystem(core::u64 universeSeed,
                               const sim::SystemStub& stub,
                               const std::vector<sim::Faction>& factions) {
  core::SplitMix64 rng(stub.seed);

  sim::StarSystem sys{};
  sys.stub = stub;

  sys.star = makeStar(stub.primaryClass, rng);

  NameGenerator ng(stub.seed);

  // Macro trade profile used to shape station types + local economy tuning.
  const TradeProfile tp = generateTradeProfile(universeSeed, stub);

  // Station-type selection bias (agri <-> industrial). Blend the faction's
  // ideology with local conditions, leaning a bit harder on local conditions
  // so the map feels regionally diverse.
  const sim::Faction* ctrlFac = findFaction(stub.factionId, factions);
  const double facBias = ctrlFac ? ctrlFac->industryBias : 0.0;
  const double typeBias = std::clamp(facBias * 0.40 + localIndustryBias(tp) * 0.90, -1.0, 1.0);

  // Planets
  double a = rng.range(0.25, 0.6); // start AU
  const int nPlanets = std::max(0, stub.planetCount);

  sys.planets.reserve(static_cast<std::size_t>(nPlanets));

  for (int i = 0; i < nPlanets; ++i) {
    sim::Planet p{};
    p.name = ng.planetName(stub.name, i);
    a *= rng.range(1.35, 1.9);
    a += rng.range(0.05, 0.25);

    p.orbit.semiMajorAxisAU = a;
    p.orbit.eccentricity = rng.range(0.0, 0.18);
    p.orbit.inclinationRad = rng.range(0.0, stellar::math::degToRad(6.0));
    p.orbit.ascendingNodeRad = rng.range(0.0, 2.0*stellar::math::kPi);
    p.orbit.argPeriapsisRad = rng.range(0.0, 2.0*stellar::math::kPi);
    p.orbit.meanAnomalyAtEpochRad = rng.range(0.0, 2.0*stellar::math::kPi);
    p.orbit.epochDays = 0.0;

    // Kepler-ish: P(years)^2 = a(AU)^3 / M(star)
    const double years = std::sqrt((a*a*a) / std::max(0.08, sys.star.massSol));
    p.orbit.periodDays = years * 365.25;

    p.type = pickPlanetType(a, rng);
    setPlanetMassRadius(p, rng);

    sys.planets.push_back(std::move(p));
  }

  // Moons (procedural-only): generated from a separate RNG stream per planet so
  // we don't perturb downstream RNG draws (stations, etc.).
  sys.moons.clear();
  for (std::size_t i = 0; i < sys.planets.size(); ++i) {
    appendMoonsForPlanet(sys.moons, stub, sys.star, sys.planets[i], i);
  }

  // ---------------------------------------------------------------------------
  // Stations
  // ---------------------------------------------------------------------------
  const int nStations = std::max(0, stub.stationCount);
  sys.stations.reserve(static_cast<std::size_t>(nStations));

  // Helper: physical/docking parameters.
  auto setStationPhysicals = [&](sim::Station& st) {
    auto rr = [&](double x, double y) { return rng.range(x, y); };

    double baseRadius = 6000.0;
    double speed = 0.20;

    switch (st.type) {
      case econ::StationType::Outpost:
        baseRadius = 4500.0;
        speed = 0.18;
        break;
      case econ::StationType::Mining:
        baseRadius = 6500.0;
        speed = 0.20;
        break;
      case econ::StationType::Refinery:
        baseRadius = 7000.0;
        speed = 0.20;
        break;
      case econ::StationType::Agricultural:
        baseRadius = 6500.0;
        speed = 0.22;
        break;
      case econ::StationType::Industrial:
        baseRadius = 8000.0;
        speed = 0.22;
        break;
      case econ::StationType::Research:
        baseRadius = 6000.0;
        speed = 0.20;
        break;
      case econ::StationType::TradeHub:
        baseRadius = 11000.0;
        speed = 0.25;
        break;
      case econ::StationType::Shipyard:
        baseRadius = 13000.0;
        speed = 0.28;
        break;
      default:
        break;
    }

    st.radiusKm = baseRadius * rr(0.85, 1.20);

    // Slot scaled off radius.
    st.slotWidthKm  = st.radiusKm * rr(0.75, 0.95);
    st.slotHeightKm = st.radiusKm * rr(0.30, 0.45);
    st.slotDepthKm  = st.radiusKm * rr(0.90, 1.35);

    // Approach corridor
    st.approachLengthKm = st.radiusKm * rr(8.0, 14.0);
    st.approachRadiusKm = st.radiusKm * rr(1.2, 2.2);
    st.maxApproachSpeedKmS = speed * rr(0.85, 1.15);

    // Comms range for clearance
    st.commsRangeKm = st.radiusKm * rr(10.0, 16.0);
  };

  // --- Station type plan ---
  // We intentionally "plan" the station mix before placing orbits so we can
  // guarantee the system expresses its macro identity (hub, mining, agri...).
  std::vector<econ::StationType> stationTypes;
  stationTypes.reserve(static_cast<std::size_t>(nStations));

  if (nStations > 0) {
    const econ::StationType primary = primaryStationType(tp, typeBias);
    stationTypes.push_back(primary);

    auto pushUnique = [&](econ::StationType t) {
      if ((int)stationTypes.size() >= nStations) return;
      if (std::find(stationTypes.begin(), stationTypes.end(), t) == stationTypes.end()) {
        stationTypes.push_back(t);
      }
    };

    // If the profile is strongly peaked, guarantee at least one matching station.
    if (tp.hub > 0.65) pushUnique(econ::StationType::TradeHub);
    if (tp.resources > 0.70) pushUnique(econ::StationType::Mining);
    if (tp.agriculture > 0.70) pushUnique(econ::StationType::Agricultural);
    if (tp.industry > 0.70) pushUnique(econ::StationType::Industrial);
    if (tp.technology > 0.75) pushUnique(econ::StationType::Research);
    if (tp.industry > 0.62 && tp.hub > 0.55) pushUnique(econ::StationType::Shipyard);
    if (tp.industry > 0.58 && tp.resources > 0.58) pushUnique(econ::StationType::Refinery);

    while ((int)stationTypes.size() < nStations) {
      stationTypes.push_back(pickStationType(rng, tp, typeBias));
    }
  }

  for (int i = 0; i < nStations; ++i) {
    sim::Station st{};
    st.id = core::hashCombine(static_cast<core::u64>(stub.id), static_cast<core::u64>(i + 1));
    st.name = ng.stationName(stub.name, i);
    st.type = stationTypes[static_cast<std::size_t>(i)];

    const core::u64 stationSeed = core::hashCombine(stub.seed, static_cast<core::u64>(st.id));
    applyFactionToStation(st, stub.factionId, factions, tp, stationSeed);

    // Place stations on simple orbits around the primary.
    // Prefer to place near an existing planet orbit (feels like inhabited space), else pick a random AU.
    double aAU = rng.range(0.8, 5.5);
    if (!sys.planets.empty()) {
      const int idx = rng.range<int>(0, (int)sys.planets.size() - 1);
      aAU = sys.planets[(std::size_t)idx].orbit.semiMajorAxisAU * rng.range(0.92, 1.08);
    }

    st.orbit.semiMajorAxisAU = std::max(0.35, aAU);
    st.orbit.eccentricity = rng.range(0.0, 0.06);
    st.orbit.inclinationRad = rng.range(0.0, stellar::math::degToRad(3.0));
    st.orbit.ascendingNodeRad = rng.range(0.0, 2.0 * stellar::math::kPi);
    st.orbit.argPeriapsisRad = rng.range(0.0, 2.0 * stellar::math::kPi);
    st.orbit.meanAnomalyAtEpochRad = rng.range(0.0, 2.0 * stellar::math::kPi);
    st.orbit.epochDays = 0.0;

    const double years = std::sqrt((st.orbit.semiMajorAxisAU * st.orbit.semiMajorAxisAU * st.orbit.semiMajorAxisAU) /
                                   std::max(0.08, sys.star.massSol));
    st.orbit.periodDays = years * 365.25;

    setStationPhysicals(st);

    sys.stations.push_back(std::move(st));
  }

  // Border outposts:
  // If a system sits inside the overlapping influence radii of multiple factions
  // (and it has enough stations), spawn a minority-owned outpost for the strongest
  // nearby rival.
  //
  // This makes contested systems naturally appear at faction borders, which then
  // feeds into the SecurityModel (contest01) and any systems that consume it
  // (encounters, mission rewards, etc.).
  //
  // IMPORTANT: This logic does NOT consume RNG; it should only change station
  // ownership/economy parameters, not orbital layouts.
  if (stub.factionId != 0 && nStations >= 2 && sys.stations.size() >= 2) {
    const sim::Faction* ctrlFac = findFaction(stub.factionId, factions);
    if (ctrlFac && ctrlFac->influenceRadiusLy > 1e-6) {
      const double dCtrl = (stub.posLy - ctrlFac->homePosLy).length();
      const double wCtrl = std::max(0.0, 1.0 - dCtrl / ctrlFac->influenceRadiusLy);

      core::u32 bestOtherId = 0;
      double bestOtherW = 0.0;
      double bestOtherD = 1e30;

      for (const auto& f : factions) {
        if (f.id == 0) continue;
        if (f.id == stub.factionId) continue;
        if (f.influenceRadiusLy <= 1e-6) continue;
        const double d = (stub.posLy - f.homePosLy).length();
        if (d >= f.influenceRadiusLy) continue;
        const double w = std::max(0.0, 1.0 - d / f.influenceRadiusLy);
        // Prefer the strongest overlapping influence; tie-break by closer distance then faction id.
        if (w > bestOtherW + 1e-12 ||
            (std::abs(w - bestOtherW) <= 1e-12 && (d < bestOtherD - 1e-9 || (std::abs(d - bestOtherD) <= 1e-9 && f.id < bestOtherId)))) {
          bestOtherW = w;
          bestOtherD = d;
          bestOtherId = f.id;
        }
      }

      if (bestOtherId != 0 && bestOtherW > 0.18) {
        const double border01 = bestOtherW / std::max(1e-9, wCtrl + bestOtherW);

        // Threshold chosen so deep-core systems remain uniform, while overlap regions
        // gain contested presence.
        if (border01 > 0.28) {
          // Keep the "primary" station (Alpha, index 0) under the controlling faction.
          //
          // For 2-station systems we can only do this safely if Alpha is strictly heavier
          // than the outpost we would flip; otherwise the computed system-control could
          // end up ambiguous (tie) or dominated by the minority.
          if (sys.stations.size() == 2) {
            const double w0 = stationControlWeight(sys.stations[0].type);
            const double w1 = stationControlWeight(sys.stations[1].type);
            if (w0 > w1 + 1e-12) {
              const core::u64 stSeed = core::hashCombine(stub.seed, static_cast<core::u64>(sys.stations[1].id));
              applyFactionToStation(sys.stations[1], bestOtherId, factions, tp, stSeed);
            }
          } else {
            // Choose a "minor" station (lowest weight) among indices [1..].
            std::size_t pickIdx = 1;
            double bestW = 1e30;
            for (std::size_t i = 1; i < sys.stations.size(); ++i) {
              const double sw = stationControlWeight(sys.stations[i].type);
              if (sw < bestW - 1e-12 || (std::abs(sw - bestW) <= 1e-12 && i > pickIdx)) {
                bestW = sw;
                pickIdx = i;
              }
            }

            const core::u64 stSeed = core::hashCombine(stub.seed, static_cast<core::u64>(sys.stations[pickIdx].id));
            applyFactionToStation(sys.stations[pickIdx], bestOtherId, factions, tp, stSeed);
          }
        }
      }
    }
  }

  return sys;
}

sim::StarSystem generateSystem(const sim::SystemStub& stub, const std::vector<sim::Faction>& factions) {
  // Legacy overload for external callers.
  // NOTE: stubs already encode the universe seed into stub.seed, so this is a
  // reasonable fallback.
  return generateSystem(stub.seed, stub, factions);
}

} // namespace stellar::proc
