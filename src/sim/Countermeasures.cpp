#include "stellar/sim/Countermeasures.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static math::Vec3d safeNormalized(const math::Vec3d& v, const math::Vec3d& fallback) {
  if (v.lengthSq() < 1e-12) return fallback;
  return v.normalized();
}

void spawnCountermeasureBurst(std::vector<Countermeasure>& ioCountermeasures,
                              core::u64& ioNextId,
                              CountermeasureType type,
                              const math::Vec3d& originKm,
                              const math::Vec3d& baseVelKmS,
                              const math::Vec3d& shipRight,
                              const math::Vec3d& shipUp,
                              const math::Vec3d& shipForward,
                              const CountermeasureBurstParams& params) {
  const int n = std::max(0, params.count);
  if (n == 0) return;

  const math::Vec3d r = safeNormalized(shipRight, {1, 0, 0});
  const math::Vec3d u = safeNormalized(shipUp, {0, 1, 0});
  const math::Vec3d f = safeNormalized(shipForward, {0, 0, 1});

  const double spread = std::max(0.0, params.spread);
  const double ttl = std::max(0.0, params.ttlSimSec);

  // Deterministic seed derived from the next id and burst params.
  core::u64 seed = ioNextId;
  seed ^= (core::u64)n * 0x9E3779B97F4A7C15ull;
  seed ^= (core::u64)type * 0xD6E8FEB86659FD93ull;
  core::SplitMix64 rng(seed);

  const math::Vec3d baseDir = (-f).normalized();

  ioCountermeasures.reserve(ioCountermeasures.size() + (std::size_t)n);

  for (int i = 0; i < n; ++i) {
    Countermeasure c{};
    c.id = ioNextId++;
    c.type = type;

    c.posKm = originKm;
    c.radiusKm = std::max(0.0, params.radiusKm);

    c.ttlSimSec = ttl;
    c.ttlMaxSimSec = ttl;

    // Scatter in the ship's local right/up plane around baseDir.
    const double sx = rng.range(-spread, spread);
    const double sy = rng.range(-spread, spread);
    math::Vec3d dir = baseDir + r * sx + u * sy;
    if (dir.lengthSq() < 1e-12) dir = baseDir;
    dir = dir.normalized();

    const double eject = std::max(0.0, params.ejectSpeedKmS);
    c.velKmS = baseVelKmS + dir * eject;

    // Strength defaults are per-burst; type picks which channel is active.
    const double heat = std::max(0.0, params.heatStrength);
    const double radar = std::max(0.0, params.radarStrength);
    if (type == CountermeasureType::Flare) {
      c.heatStrength = heat;
      c.radarStrength = 0.0;
    } else {
      c.heatStrength = 0.0;
      c.radarStrength = radar;
    }

    ioCountermeasures.push_back(c);
  }
}

void stepCountermeasures(std::vector<Countermeasure>& ioCountermeasures, double dtSim) {
  if (dtSim <= 0.0) return;
  if (ioCountermeasures.empty()) return;

  for (auto& c : ioCountermeasures) {
    if (c.ttlSimSec <= 0.0) continue;
    c.posKm = c.posKm + c.velKmS * dtSim;
    c.ttlSimSec -= dtSim;
  }

  ioCountermeasures.erase(
    std::remove_if(ioCountermeasures.begin(), ioCountermeasures.end(),
                   [](const Countermeasure& c) { return c.ttlSimSec <= 0.0; }),
    ioCountermeasures.end());
}

void appendCountermeasureTargets(const std::vector<Countermeasure>& countermeasures,
                                 std::vector<SphereTarget>& ioTargets) {
  if (countermeasures.empty()) return;

  const std::size_t start = ioTargets.size();
  ioTargets.reserve(start + countermeasures.size());

  for (std::size_t i = 0; i < countermeasures.size(); ++i) {
    const Countermeasure& c = countermeasures[i];

    SphereTarget t{};
    t.kind = CombatTargetKind::Decoy;
    t.index = start + i;
    t.id = c.id;

    t.centerKm = c.posKm;
    t.velKmS = c.velKmS;
    t.radiusKm = std::max(0.0, c.radiusKm);

    const double life = (c.ttlMaxSimSec > 1e-9) ? std::clamp(c.ttlSimSec / c.ttlMaxSimSec, 0.0, 1.0) : 0.0;
    t.decoyHeat = std::max(0.0, c.heatStrength) * life;
    t.decoyRadar = std::max(0.0, c.radarStrength) * life;

    ioTargets.push_back(t);
  }
}

} // namespace stellar::sim
