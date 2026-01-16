#include "SensorContacts.h"

#include <algorithm>
#include <cmath>

namespace stellar::game {
namespace {

constexpr double kDaysToSeconds = 86400.0;

inline double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

inline double safeExpDecay(double dtSec, double tauSec) {
  if (tauSec <= 1.0e-6) return 0.0;
  if (dtSec <= 0.0) return 1.0;
  return std::exp(-dtSec / tauSec);
}

inline bool segmentIntersectsSphere(const math::Vec3d& a,
                                    const math::Vec3d& b,
                                    const math::Vec3d& c,
                                    double r) {
  const math::Vec3d ab = b - a;
  const double abLenSq = ab.lengthSq();
  if (abLenSq < 1.0e-12) return false;

  const math::Vec3d ac = c - a;
  double t = math::dot(ac, ab) / abLenSq;
  if (t < 0.0) t = 0.0;
  if (t > 1.0) t = 1.0;
  const math::Vec3d p = a + ab * t;
  const double dSq = (p - c).lengthSq();
  return dSq <= r * r;
}

} // namespace

void SensorContacts::reset() {
  tracks_.clear();
  views_.clear();
  pingUntilDays_ = 0.0;
  nextPingAllowedDays_ = 0.0;
}

bool SensorContacts::tryStartPing(double nowDays, const SensorContactsSettings& settings) {
  if (nowDays < nextPingAllowedDays_) return false;
  const double durDays = settings.pingDurationSec / kDaysToSeconds;
  const double cdDays = settings.pingCooldownSec / kDaysToSeconds;
  pingUntilDays_ = nowDays + durDays;
  nextPingAllowedDays_ = nowDays + cdDays;
  return true;
}

void SensorContacts::boostIdentify(core::u64 id, double amount01) {
  amount01 = clamp01(amount01);
  SensorTrack& t = tracks_[id];
  t.id = id;
  t.identify01 = std::max(t.identify01, amount01);
}

double SensorContacts::computeDetection01(const math::Vec3d& sensorPosKm,
                                          const SensorTarget& target,
                                          const std::vector<SensorOccluder>& occluders,
                                          double nowDays,
                                          const SensorContactsSettings& settings) const {
  if (settings.baseRangeKm <= 0.0) return 0.0;
  if (target.id == 0) return 0.0;

  const math::Vec3d delta = target.positionKm - sensorPosKm;
  const double rangeKm = std::sqrt(delta.lengthSq());

  const double sig = std::max(0.05, target.signature);
  const double effectiveRangeKm = settings.baseRangeKm * sig;
  if (effectiveRangeKm <= 1.0e-6) return 0.0;

  const double r01 = rangeKm / effectiveRangeKm;
  if (r01 >= 1.0) return 0.0;

  // Smooth-ish falloff: strong near, fades toward 0 at the edge.
  double d01 = (1.0 - r01);
  d01 = d01 * d01;

  // Occlusion attenuation (simple LOS sphere test against large bodies).
  if (settings.occlusionEnabled && settings.occlusionPenalty01 > 0.0) {
    const double pen = clamp01(settings.occlusionPenalty01);
    for (const auto& occ : occluders) {
      if (occ.radiusKm <= 0.0) continue;
      if (segmentIntersectsSphere(sensorPosKm, target.positionKm, occ.positionKm, occ.radiusKm)) {
        d01 *= (1.0 - pen);
        break;
      }
    }
  }

  // Active ping boosts detection for a short window.
  if (pingActive(nowDays) && settings.pingBoost01 > 0.0) {
    d01 = clamp01(d01 + settings.pingBoost01);
  }

  return clamp01(d01);
}

std::vector<SensorView> SensorContacts::update(double nowDays,
                                               double dtSec,
                                               const math::Vec3d& sensorPosKm,
                                               const std::vector<SensorTarget>& targets,
                                               const std::vector<SensorOccluder>& occluders,
                                               const SensorContactsSettings& settings) {
  // Decay all tracks first.
  const double memDecay = safeExpDecay(dtSec, settings.trackMemorySec);
  const double idDecay = safeExpDecay(dtSec, settings.identifyDecaySec);
  for (auto& kv : tracks_) {
    SensorTrack& t = kv.second;
    t.detectedNow = false;
    t.confidence01 *= memDecay;
    t.identify01 *= idDecay;
    if (t.confidence01 < 1.0e-6) t.confidence01 = 0.0;
    if (t.identify01 < 1.0e-6) t.identify01 = 0.0;
  }

  // Refresh tracks from targets.
  const double idGain = (settings.identifyTimeSec <= 1.0e-6) ? 1.0 : (dtSec / settings.identifyTimeSec);
  for (const auto& tar : targets) {
    if (tar.id == 0) continue;

    const math::Vec3d delta = tar.positionKm - sensorPosKm;
    const double rangeKm = std::sqrt(delta.lengthSq());
    const double d01 = computeDetection01(sensorPosKm, tar, occluders, nowDays, settings);

    SensorTrack& t = tracks_[tar.id];
    t.id = tar.id;
    t.lastPosKm = tar.positionKm;
    t.lastRangeKm = rangeKm;

    if (d01 > 0.0) {
      t.detectedNow = true;
      t.lastSeenDays = nowDays;
      // Keep the "best" confidence this frame; confidence otherwise decays.
      t.confidence01 = std::max(t.confidence01, d01);
      t.identify01 = clamp01(t.identify01 + d01 * idGain);
    }

    t.identified = (t.identify01 >= settings.identifyThreshold01);
  }

  // Cull stale tracks (keep small memory footprint).
  {
    const double maxAgeSec = std::max(5.0, settings.trackMemorySec * 3.0);
    for (auto it = tracks_.begin(); it != tracks_.end();) {
      const SensorTrack& t = it->second;
      const double ageSec = (nowDays - t.lastSeenDays) * kDaysToSeconds;
      const bool stale = (ageSec > maxAgeSec) && (t.confidence01 <= 0.01);
      if (stale) {
        it = tracks_.erase(it);
      } else {
        ++it;
      }
    }
  }

  // Build views.
  views_.clear();
  views_.reserve(tracks_.size());
  for (const auto& kv : tracks_) {
    const SensorTrack& t = kv.second;
    if (t.confidence01 <= 0.0) continue;
    SensorView v{};
    v.id = t.id;
    v.positionKm = t.lastPosKm;
    v.rangeKm = t.lastRangeKm;
    v.confidence01 = t.confidence01;
    v.identify01 = t.identify01;
    v.identified = t.identified;
    v.detectedNow = t.detectedNow;
    v.ageSec = (nowDays - t.lastSeenDays) * kDaysToSeconds;
    views_.push_back(v);
  }

  std::sort(views_.begin(), views_.end(), [](const SensorView& a, const SensorView& b) {
    if (a.detectedNow != b.detectedNow) return a.detectedNow > b.detectedNow;
    return a.confidence01 > b.confidence01;
  });

  return views_;
}

double SensorContacts::pingCooldownRemainingSec(double nowDays) const {
  const double remDays = nextPingAllowedDays_ - nowDays;
  if (remDays <= 0.0) return 0.0;
  return remDays * kDaysToSeconds;
}

} // namespace stellar::game
