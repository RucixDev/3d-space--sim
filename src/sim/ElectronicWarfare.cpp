#include "stellar/sim/ElectronicWarfare.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include <cmath>
#include <numbers>

namespace stellar::sim {

static inline double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

double computeJammingSnr(double distKm,
                         double jammerPower,
                         const EwJammerParams& params) {
  distKm = std::max(distKm, params.minDistKm);
  jammerPower = std::max(0.0, jammerPower);

  const double r = std::max(params.halfRangeKm, 1e-6);
  const double snr = jammerPower * (r * r) / (distKm * distKm);
  if (!std::isfinite(snr) || snr <= 0.0) return 0.0;
  return snr;
}

double jamming01FromSnr(double snr) {
  if (!std::isfinite(snr) || snr <= 0.0) return 0.0;
  return clamp01(snr / (1.0 + snr));
}

double applyJammingToSensorPower(double baseSensorPower,
                                 double jamming01,
                                 bool pingActive,
                                 const EwJammerParams& params) {
  baseSensorPower = std::max(0.0, baseSensorPower);
  jamming01 = clamp01(jamming01);

  double effJ = jamming01;
  if (pingActive) {
    effJ *= clamp01(params.pingJammingMult);
  }

  const double g = std::max(0.0, params.suppressionGain);
  const double denom = 1.0 + g * effJ;
  if (!std::isfinite(denom) || denom <= 0.0) return baseSensorPower;
  return baseSensorPower / denom;
}

std::vector<EwGhostBlip> generateGhostBlips(core::u64 seed,
                                           double timeSec,
                                           double rangeKm,
                                           double jamming01,
                                           const EwGhostParams& params) {
  std::vector<EwGhostBlip> out;

  jamming01 = clamp01(jamming01);
  if (jamming01 <= 0.0) return out;

  rangeKm = std::max(1.0, rangeKm);

  const int maxN = std::max(0, params.maxBlips);
  const int n = std::clamp((int)std::llround((double)maxN * (0.25 + 0.75 * jamming01)), 1, maxN);
  if (n <= 0) return out;

  const double hz = std::max(1e-6, params.reseedHz);
  const double bucketDur = 1.0 / hz;
  const double t = std::isfinite(timeSec) ? std::max(0.0, timeSec) : 0.0;
  const core::u64 bucket = (core::u64)std::floor(t / bucketDur);
  const double frac = std::clamp((t - (double)bucket * bucketDur) / bucketDur, 0.0, 1.0);

  const core::u64 s = core::hashCombine(seed, core::hashCombine(core::fnv1a64("ew_ghost_bucket"), bucket + 1u));
  core::SplitMix64 rng(s);

  out.reserve((std::size_t)n);

  const double drift = std::max(0.0, params.driftKmPerSec) * bucketDur;
  const double minS = clamp01(params.minStrength01);
  const double maxS = std::max(minS, clamp01(params.maxStrength01));

  for (int i = 0; i < n; ++i) {
    // Base position: uniform in a disc (sqrt radial distribution).
    const double u = rng.nextUnit();
    const double v = rng.nextUnit();
    const double ang = 2.0 * std::numbers::pi * u;
    const double rad = std::sqrt(v) * rangeKm;

    double x = std::cos(ang) * rad;
    double z = std::sin(ang) * rad;

    // Drift direction is independent of base position.
    const double u2 = rng.nextUnit();
    const double ang2 = 2.0 * std::numbers::pi * u2;
    const double dx = std::cos(ang2);
    const double dz = std::sin(ang2);

    // Centered drift: [-0.5..0.5] so patterns don't "walk" too far.
    const double driftK = (frac - 0.5) * drift;
    x += dx * driftK;
    z += dz * driftK;

    // Clamp inside circle.
    const double d2 = x * x + z * z;
    const double r2 = rangeKm * rangeKm;
    if (d2 > r2 && d2 > 1e-9) {
      const double d = std::sqrt(d2);
      const double k = rangeKm / d;
      x *= k;
      z *= k;
    }

    const double s01 = (minS + (maxS - minS) * rng.nextUnit()) * (0.35 + 0.65 * jamming01);

    EwGhostBlip b{};
    b.xKm = x;
    b.zKm = z;
    b.strength01 = clamp01(s01);
    b.id = core::hashCombine(s, (core::u64)(i + 1));
    out.push_back(b);
  }

  return out;
}

} // namespace stellar::sim
