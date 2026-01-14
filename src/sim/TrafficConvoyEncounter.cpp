#include "stellar/sim/TrafficConvoyEncounter.h"

#include "stellar/core/Hash.h"
#include "stellar/math/Math.h"
#include "stellar/math/Quat.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {
namespace {

static double clamp01(double x) {
  return std::clamp(x, 0.0, 1.0);
}

static math::Quatd randomUnitQuat(core::SplitMix64& rng) {
  // Uniform random quaternion (see: Shoemake 1992).
  const double u1 = rng.nextUnit();
  const double u2 = rng.nextUnit();
  const double u3 = rng.nextUnit();

  const double a = std::sqrt(std::max(0.0, 1.0 - u1));
  const double b = std::sqrt(std::max(0.0, u1));
  const double t1 = 2.0 * math::kPi * u2;
  const double t2 = 2.0 * math::kPi * u3;

  const double x = a * std::sin(t1);
  const double y = a * std::cos(t1);
  const double z = b * std::sin(t2);
  const double w = b * std::cos(t2);
  return math::Quatd{w, x, y, z}.normalized();
}

static std::vector<math::Vec3d> fibonacciSphereDirs(int n) {
  std::vector<math::Vec3d> dirs;
  if (n <= 0) return dirs;
  dirs.reserve((std::size_t)n);

  // Golden angle (radians)
  const double ga = math::kPi * (3.0 - std::sqrt(5.0));

  for (int i = 0; i < n; ++i) {
    const double t = (double)i + 0.5;
    const double y = 1.0 - 2.0 * (t / (double)n);
    const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
    const double phi = ga * (double)i;
    const double x = std::cos(phi) * r;
    const double z = std::sin(phi) * r;
    dirs.push_back({x, y, z});
  }
  return dirs;
}

static std::vector<math::Vec3d> rotatedDirs(int n, core::SplitMix64& rng) {
  std::vector<math::Vec3d> dirs = fibonacciSphereDirs(n);
  if (dirs.empty()) return dirs;

  const math::Quatd q = randomUnitQuat(rng);
  for (auto& d : dirs) d = q.rotate(d);
  return dirs;
}

static core::u64 trafficEncounterSeed(core::u64 universeSeed, core::u64 convoyId) {
  // Keep the legacy tag so existing saves/debug captures retain roughly
  // similar encounter patterns even after refactors.
  return core::hashCombine(universeSeed ^ core::fnv1a64("traffic_convoy_spawn"), convoyId);
}

} // namespace

TrafficConvoyEncounterPlan planTrafficConvoyEncounter(core::u64 universeSeed,
                                                      core::u64 convoyId,
                                                      double cargoValueCr,
                                                      const SystemSecurityProfile& sec,
                                                      bool lawfulSystem,
                                                      int alivePiratesNow,
                                                      const TrafficConvoyEncounterParams& params) {
  TrafficConvoyEncounterPlan out{};
  if (convoyId == 0) return out;

  out.seed = trafficEncounterSeed(universeSeed, convoyId);

  // Use sub-stream RNGs so that adding/removing cosmetic randomness (e.g.
  // formation placement) does not change the *decision* of whether an ambush
  // occurs. This makes the system more stable for gameplay tuning.
  core::SplitMix64 rngDec(core::hashCombine(out.seed, core::fnv1a64("dec")));
  core::SplitMix64 rngLoad(core::hashCombine(out.seed, core::fnv1a64("load")));
  core::SplitMix64 rngPos(core::hashCombine(out.seed, core::fnv1a64("pos")));

  out.cargoValueCr = std::max(0.0, cargoValueCr);
  const double denom = std::max(1e-6, params.valueScaleCr);
  out.value01 = clamp01(out.cargoValueCr / denom);

  out.risk01 = std::clamp(params.pirateRiskBase +
                            params.pirateRiskPiracyWeight * sec.piracy01 +
                            params.pirateRiskValueWeight * out.value01 -
                            params.pirateRiskSecurityWeight * sec.security01,
                          0.0,
                          params.pirateRiskMax);

  // ---- Trader (convoy leader) --------------------------------------------
  {
    TrafficEncounterNpcSpec t{};
    t.role = TrafficEncounterNpcRole::Trader;
    t.leader = true;

    // Position: small offset from the signal so the convoy isn't perfectly on
    // the marker (feels less game-y).
    const auto dir = rotatedDirs(1, rngPos);
    const double r = rngPos.range(params.convoyOffsetMinKm, params.convoyOffsetMaxKm);
    t.posOffsetKm = dir.empty() ? math::Vec3d{r, 0, 0} : dir[0] * r;
    t.velOffsetKmS = {0, 0, 0};

    // Loadout: a lightly armed hauler.
    const int thrMk = 1 + (out.value01 > 0.65);
    const int shieldMk = 1 + (sec.security01 > 0.55);
    const int distMk = 1 + (sec.security01 > 0.75);
    t.hullClass = ShipHullClass::Hauler;
    t.thrusterMk = thrMk;
    t.shieldMk = shieldMk;
    t.distributorMk = distMk;
    t.weapon = WeaponType::MiningLaser;

    t.hullMul = 0.70 + 0.10 * out.value01;
    t.shieldMul = 0.70 + 0.08 * out.value01;
    t.regenMul = 0.0;
    t.accelMul = 0.78;
    t.aiSkill = 0.42;

    t.pips = {4, 1, 1};
    normalizePips(t.pips);

    out.npcs.push_back(std::move(t));
  }

  // ---- Escorts (lawful systems) ------------------------------------------
  out.escortCount = 0;
  if (lawfulSystem) {
    if (sec.security01 > params.escortSecurity2) out.escortCount = 2;
    else if (sec.security01 > params.escortSecurity1) out.escortCount = 1;

    if (out.value01 > params.escortExtraValueThresh && out.escortCount < params.escortMax) {
      if (rngDec.chance(params.escortExtraChance)) ++out.escortCount;
    }
  }
  out.escortCount = std::clamp(out.escortCount, 0, params.escortMax);

  if (out.escortCount > 0) {
    const auto dirs = rotatedDirs(out.escortCount, rngPos);

    for (int i = 0; i < out.escortCount; ++i) {
      TrafficEncounterNpcSpec p{};
      p.role = TrafficEncounterNpcRole::Police;
      p.leader = (i == 0);

      const math::Vec3d d = dirs.empty() ? math::Vec3d{1, 0, 0} : dirs[(std::size_t)i % dirs.size()];
      const double r = rngPos.range(params.escortOffsetMinKm, params.escortOffsetMaxKm);

      // Offset around the convoy leader.
      p.posOffsetKm = out.npcs[0].posOffsetKm + d * r;

      // Small velocity jitter so the formation isn't frozen.
      const auto vd = rotatedDirs(1, rngPos);
      const double vr = rngPos.range(0.0, params.escortVelJitterMaxKmS);
      p.velOffsetKmS = vd.empty() ? math::Vec3d{0, 0, 0} : vd[0] * vr;

      const int tier = std::clamp(1 + (out.value01 > 0.65) + (sec.security01 > 0.75), 1, 3);
      const WeaponType w = (tier >= 3 && p.leader) ? WeaponType::Railgun : WeaponType::PulseLaser;

      const int tMk = 2;
      const int sMk = 1;
      const int dMk = tier;

      // Compute regenMul so the absolute regen is stable across tiers.
      // (Mirrors the prototype tuning in apps/stellar_game/main.cpp.)
      const ShipDerivedStats dsTmp = computeShipDerivedStats(ShipHullClass::Scout, tMk, sMk, dMk);
      const double regenPerSimMin = std::max(0.0, params.npcShieldRegenPerSec) * 60.0;
      const double targetRegenPerSimMin = regenPerSimMin * (1.0 + 0.12 * (double)(tier - 1));
      const double regenMul = (dsTmp.shieldRegenPerSimMin > 1e-9) ? (targetRegenPerSimMin / dsTmp.shieldRegenPerSimMin) : 0.0;

      p.hullClass = ShipHullClass::Scout;
      p.thrusterMk = tMk;
      p.shieldMk = sMk;
      p.distributorMk = dMk;
      p.weapon = w;

      p.hullMul = 0.86 * (1.0 + 0.10 * (double)(tier - 1));
      p.shieldMul = 0.88 * (1.0 + 0.10 * (double)(tier - 1));
      p.regenMul = regenMul;
      p.accelMul = 0.80 * (1.0 + 0.10 * (double)(tier - 1));
      p.aiSkill = p.leader ? 0.70 : 0.60;

      p.pips = {2, 2, 2};
      normalizePips(p.pips);

      out.npcs.push_back(std::move(p));
    }
  }

  // ---- Raiders ------------------------------------------------------------
  out.ambush = false;
  out.pirateCount = 0;
  if (alivePiratesNow < params.pirateAliveCap && rngDec.chance(out.risk01)) {
    out.ambush = true;

    int count = params.pirateMin;
    // At least 2, with a mild chance for 3.
    if (rngDec.chance(0.55)) ++count;

    if (out.risk01 > params.pirateExtra1Risk && rngDec.chance(params.pirateExtra1Chance)) ++count;
    if (out.risk01 > params.pirateExtra2Risk && rngDec.chance(params.pirateExtra2Chance)) ++count;
    out.pirateCount = std::clamp(count, params.pirateMin, params.pirateMax);

    const auto dirs = rotatedDirs(out.pirateCount, rngPos);

    for (int i = 0; i < out.pirateCount; ++i) {
      const bool leader = (i == 0);
      TrafficEncounterNpcSpec p{};
      p.role = TrafficEncounterNpcRole::Pirate;
      p.leader = leader;

      const math::Vec3d d = dirs.empty() ? math::Vec3d{1, 0, 0} : dirs[(std::size_t)i % dirs.size()];
      const double r = rngPos.range(params.pirateOffsetMinKm, params.pirateOffsetMaxKm);
      p.posOffsetKm = out.npcs[0].posOffsetKm + d * r;

      const auto vd = rotatedDirs(1, rngPos);
      const double vr = rngPos.range(0.0, params.pirateVelJitterMaxKmS);
      p.velOffsetKmS = vd.empty() ? math::Vec3d{0, 0, 0} : vd[0] * vr;

      const double statMul = leader ? 1.22 : (0.88 + 0.22 * rngLoad.nextUnit());
      const int tMk = leader ? 2 : 1;
      const int dMk = leader ? 2 : 1;

      WeaponType w = WeaponType::Cannon;
      const double rW = rngLoad.nextUnit();
      if (leader) w = (rW < 0.55) ? WeaponType::Railgun : WeaponType::Cannon;
      else w = (rW < 0.55) ? WeaponType::Cannon : WeaponType::BeamLaser;

      p.hullClass = ShipHullClass::Fighter;
      p.thrusterMk = tMk;
      p.shieldMk = 1;
      p.distributorMk = dMk;
      p.weapon = w;

      p.hullMul = 0.95 * statMul;
      p.shieldMul = 0.67 * statMul;
      p.regenMul = leader ? 0.70 : 0.63;
      p.accelMul = leader ? 0.72 : 0.70;
      p.aiSkill = leader ? 0.68 : 0.55;

      p.pips = {1, 3, 2};
      normalizePips(p.pips);

      out.npcs.push_back(std::move(p));
    }
  }

  out.ok = !out.npcs.empty();
  return out;
}

} // namespace stellar::sim
