#include "stellar/sim/SalvageSweep.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

SalvageSweepComputer::SalvageSweepComputer() = default;

void SalvageSweepComputer::reset() {
  targetId_ = 0;
  targetScore_ = 0.0;
  lastSwitchTimeSec_ = -1.0e9;
  skip_.clear();
}

void SalvageSweepComputer::skip(core::u64 id, double timeSec) {
  if (id == 0) id = targetId_;
  if (id == 0) return;

  const double until = timeSec + std::max(0.0, params_.skipSeconds);

  // Replace any existing entry for this id.
  for (auto& e : skip_) {
    if (e.id == id) {
      e.untilTimeSec = std::max(e.untilTimeSec, until);
      goto done;
    }
  }
  skip_.push_back(SkipEntry{id, until});

done:
  // If we skipped the active target, clear it so the next update picks a new one.
  if (id == targetId_) {
    targetId_ = 0;
    targetScore_ = 0.0;
  }
}


template <typename SkipEntryT>
static bool isSkipped(core::u64 id, const std::vector<SkipEntryT>& skip, double timeSec) {
  if (id == 0) return true;
  for (const auto& e : skip) {
    if (e.id == id && e.untilTimeSec > timeSec) return true;
  }
  return false;
}

static double valueTerm(double valueCr, double valueScaleCr) {
  // Map value (credits) into a [0..1) curve with diminishing returns.
  // value=0 => 0
  // value=scale => ~0.63
  // value>>scale => -> 1
  const double s = std::max(0.001, valueScaleCr);
  const double x = std::max(0.0, valueCr) / s;
  return 1.0 - std::exp(-x);
}

static double distTerm(double distKm, double distScaleKm) {
  // Map distance into [0..1) (0 near, ->1 far).
  const double s = std::max(1.0, distScaleKm);
  const double d = std::max(0.0, distKm);
  return d / (d + s);
}

static double scoreCandidate(const SalvagePodCandidate& c,
                             const math::Vec3d& shipPosKm,
                             const SalvageSweepParams& p) {
  const double dKm = (c.posKm - shipPosKm).length();
  const double vCr = c.estimatedValueCr;

  const double mission = (c.missionId != 0) ? p.missionBonus : 0.0;
  const double v = p.valueWeight * valueTerm(vCr, p.valueScaleCr);
  const double d = p.distWeight * distTerm(dKm, p.distScaleKm);

  return mission + v - d;
}

SalvageSweepStatus SalvageSweepComputer::update(double timeSec,
                                                const math::Vec3d& shipPosKm,
                                                const std::vector<SalvagePodCandidate>& candidates) {
  // Prune expired skip entries.
  skip_.erase(std::remove_if(skip_.begin(), skip_.end(),
                             [&](const SkipEntry& e) {
                               return e.id == 0 || e.untilTimeSec <= timeSec;
                             }),
              skip_.end());

  // Find best candidate.
  int bestIdx = -1;
  double bestScore = -1.0e99;
  core::u64 bestId = 0;

  for (int i = 0; i < (int)candidates.size(); ++i) {
    const auto& c = candidates[(std::size_t)i];
    if (c.id == 0) continue;
    if (c.units <= 1e-9) continue;
    if (isSkipped(c.id, skip_, timeSec)) continue;

    const double s = scoreCandidate(c, shipPosKm, params_);
    if (s > bestScore) {
      bestScore = s;
      bestIdx = i;
      bestId = c.id;
    }
  }

  // Update current target score if it still exists.
  int curIdx = -1;
  double curScore = -1.0e99;
  if (targetId_ != 0) {
    for (int i = 0; i < (int)candidates.size(); ++i) {
      const auto& c = candidates[(std::size_t)i];
      if (c.id != targetId_) continue;
      if (c.units <= 1e-9) break;
      if (isSkipped(c.id, skip_, timeSec)) break;

      curIdx = i;
      curScore = scoreCandidate(c, shipPosKm, params_);
      break;
    }
    if (curIdx < 0) {
      targetId_ = 0;
      targetScore_ = 0.0;
    } else {
      targetScore_ = curScore;
    }
  }

  // Acquire if needed.
  if (targetId_ == 0) {
    if (bestIdx >= 0) {
      targetId_ = bestId;
      targetScore_ = bestScore;
      lastSwitchTimeSec_ = timeSec;
      return SalvageSweepStatus{targetId_, bestIdx, bestScore, true};
    }
    return SalvageSweepStatus{};
  }

  // Consider switching.
  if (bestIdx >= 0 && bestId != targetId_) {
    const double dt = timeSec - lastSwitchTimeSec_;
    const bool canSwitch = dt >= std::max(0.0, params_.minHoldSeconds);
    const bool isMuchBetter = bestScore > (targetScore_ + std::max(0.0, params_.switchHysteresis));

    if (canSwitch && isMuchBetter) {
      targetId_ = bestId;
      targetScore_ = bestScore;
      lastSwitchTimeSec_ = timeSec;
      return SalvageSweepStatus{targetId_, bestIdx, bestScore, true};
    }
  }

  return SalvageSweepStatus{targetId_, curIdx, targetScore_, true};
}

} // namespace stellar::sim
