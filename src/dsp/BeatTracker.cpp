#include "stellar/dsp/BeatTracker.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace stellar::dsp {
namespace {

static float clamp01(float v) {
  return std::clamp(v, 0.0f, 1.0f);
}

// Standard 1-pole EMA alpha for a given time constant.
// y += alpha * (x - y)
static float emaAlpha(float dtSec, float tauSec) {
  if (!std::isfinite(dtSec) || dtSec <= 0.0f) return 1.0f;
  if (!std::isfinite(tauSec) || tauSec <= 1.0e-6f) return 1.0f;
  return 1.0f - std::exp(-dtSec / tauSec);
}

static float safeDiv(float a, float b, float fallback) {
  if (!std::isfinite(a) || !std::isfinite(b) || std::abs(b) <= 1.0e-12f) return fallback;
  return a / b;
}

static float median(std::vector<float>& xs) {
  if (xs.empty()) return 0.0f;
  const std::size_t n = xs.size();
  const std::size_t mid = n / 2;
  std::nth_element(xs.begin(), xs.begin() + (std::ptrdiff_t)mid, xs.end());
  float m = xs[mid];
  if ((n & 1u) == 0u) {
    // For even counts, average the two middle elements.
    const float upper = m;
    std::nth_element(xs.begin(), xs.begin() + (std::ptrdiff_t)(mid - 1), xs.end());
    const float lower = xs[mid - 1];
    m = 0.5f * (lower + upper);
  }
  return m;
}

} // namespace

void BeatTracker::setParams(const BeatTrackerParams& pIn) {
  BeatTrackerParams p = pIn;

  // Sanitize.
  if (!std::isfinite(p.meanTauSec) || p.meanTauSec < 0.0f) p.meanTauSec = 0.50f;
  if (!std::isfinite(p.devTauSec) || p.devTauSec < 0.0f) p.devTauSec = 0.50f;

  if (!std::isfinite(p.sensitivity) || p.sensitivity < 0.0f) p.sensitivity = 1.60f;

  if (!std::isfinite(p.minBeatIntervalSec) || p.minBeatIntervalSec < 0.0f) p.minBeatIntervalSec = 0.25f;

  if (!std::isfinite(p.beatPulseDecaySec) || p.beatPulseDecaySec <= 0.0f) p.beatPulseDecaySec = 0.25f;

  if (!std::isfinite(p.bpmMin) || p.bpmMin <= 0.0f) p.bpmMin = 55.0f;
  if (!std::isfinite(p.bpmMax) || p.bpmMax <= 0.0f) p.bpmMax = 200.0f;
  if (p.bpmMax < p.bpmMin) std::swap(p.bpmMin, p.bpmMax);

  p.bpmHistory = std::clamp(p.bpmHistory, 2, 64);

  params_ = p;

  // If history capacity changed, clamp.
  if ((int)intervals_.size() > params_.bpmHistory) {
    intervals_.erase(intervals_.begin(), intervals_.begin() + ((int)intervals_.size() - params_.bpmHistory));
  }

  // Re-evaluate tempo estimate.
  updateTempo_();
}

void BeatTracker::reset() {
  last_ = {};
  prevMag_.clear();

  timeSec_ = 0.0;

  onsetPrev2_ = onsetPrev1_ = 0.0f;
  thrPrev2_ = thrPrev1_ = 0.0f;
  tPrev2_ = tPrev1_ = 0.0;

  mean_ = 0.0f;
  dev_ = 0.0f;

  lastBeatSec_ = -1.0;
  intervals_.clear();
  bpm_ = 0.0f;
  bpmConf01_ = 0.0f;

  beatPulse01_ = 0.0f;
}

BeatTrackerOutput BeatTracker::processSpectrum(const float* mag, int count, float dtSec) {
  BeatTrackerOutput out{};

  if (!std::isfinite(dtSec) || dtSec <= 0.0f) {
    dtSec = 1.0f / 60.0f;
  }

  // Keep an internal clock even for empty inputs so the beatPulse can decay.
  timeSec_ += (double)dtSec;

  // Beat pulse decay (linear).
  if (params_.beatPulseDecaySec > 1.0e-6f) {
    beatPulse01_ = std::max(0.0f, beatPulse01_ - dtSec / params_.beatPulseDecaySec);
  } else {
    beatPulse01_ = 0.0f;
  }

  // Sanitize spectrum.
  if (!mag || count <= 0) {
    // Still update output with decayed beatPulse and statistics.
    out.onset = 0.0f;
    out.mean = mean_;
    out.deviation = dev_;
    out.threshold = mean_ + params_.sensitivity * dev_;
    out.beatPulse01 = beatPulse01_;
    out.bpm = bpm_;
    out.bpmConfidence01 = bpmConf01_;
    last_ = out;

    // Shift peak-pick buffers with a zero onset.
    const float onsetCur = 0.0f;
    const float thrCur = out.threshold;

    // Peak picking against an all-zero stream can't trigger.
    onsetPrev2_ = onsetPrev1_;
    onsetPrev1_ = onsetCur;
    thrPrev2_ = thrPrev1_;
    thrPrev1_ = thrCur;
    tPrev2_ = tPrev1_;
    tPrev1_ = timeSec_;

    return last_;
  }

  if ((int)prevMag_.size() != count) {
    prevMag_.assign((std::size_t)count, 0.0f);
  }

  // Spectral flux: sum of positive changes.
  float flux = 0.0f;
  for (int i = 0; i < count; ++i) {
    const float cur = std::isfinite(mag[i]) ? mag[i] : 0.0f;
    const float prev = prevMag_[(std::size_t)i];
    const float d = cur - prev;
    if (d > 0.0f) flux += d;
    prevMag_[(std::size_t)i] = cur;
  }

  if (!std::isfinite(flux) || flux < 0.0f) flux = 0.0f;

  // Adaptive statistics.
  {
    const float aMean = emaAlpha(dtSec, params_.meanTauSec);
    mean_ += aMean * (flux - mean_);

    const float aDev = emaAlpha(dtSec, params_.devTauSec);
    const float absDev = std::abs(flux - mean_);
    dev_ += aDev * (absDev - dev_);

    if (!std::isfinite(mean_) || mean_ < 0.0f) mean_ = 0.0f;
    if (!std::isfinite(dev_) || dev_ < 0.0f) dev_ = 0.0f;
  }

  const float threshold = mean_ + params_.sensitivity * dev_;

  // Peak picking: detect whether onsetPrev1_ was a local maximum above thrPrev1_.
  bool beat = false;
  {
    const float onsetCur = flux;

    const bool have3 = (timeSec_ > 2.0 * (double)dtSec);
    if (have3) {
      const bool isLocalMax = (onsetPrev1_ > onsetPrev2_) && (onsetPrev1_ >= onsetCur);
      const bool aboveThr = (onsetPrev1_ > thrPrev1_);

      if (isLocalMax && aboveThr) {
        const double beatTime = tPrev1_;

        const double sinceLast = (lastBeatSec_ >= 0.0) ? (beatTime - lastBeatSec_) : std::numeric_limits<double>::infinity();
        if (sinceLast >= (double)params_.minBeatIntervalSec) {
          // Accept.
          beat = true;

          const double prevBeat = lastBeatSec_;
          lastBeatSec_ = beatTime;

          beatPulse01_ = 1.0f;

          if (prevBeat >= 0.0) {
            const float intervalSec = (float)(beatTime - prevBeat);
            const float bpm = safeDiv(60.0f, intervalSec, 0.0f);

            if (bpm >= params_.bpmMin && bpm <= params_.bpmMax && std::isfinite(intervalSec) && intervalSec > 0.0f) {
              intervals_.push_back(intervalSec);
              if ((int)intervals_.size() > params_.bpmHistory) {
                intervals_.erase(intervals_.begin());
              }
              updateTempo_();
            }
          }
        }
      }
    }

    // Shift peak-pick history.
    onsetPrev2_ = onsetPrev1_;
    onsetPrev1_ = onsetCur;
    thrPrev2_ = thrPrev1_;
    thrPrev1_ = threshold;

    tPrev2_ = tPrev1_;
    tPrev1_ = timeSec_;
  }

  out.onset = flux;
  out.mean = mean_;
  out.deviation = dev_;
  out.threshold = threshold;
  out.beat = beat;
  out.beatPulse01 = beatPulse01_;
  out.bpm = bpm_;
  out.bpmConfidence01 = bpmConf01_;

  last_ = out;
  return last_;
}

void BeatTracker::updateTempo_() {
  bpm_ = 0.0f;
  bpmConf01_ = 0.0f;

  if (intervals_.size() < 2) return;

  // Median interval gives robustness to occasional mis-detections.
  std::vector<float> tmp = intervals_;
  const float med = median(tmp);
  if (!std::isfinite(med) || med <= 0.0f) return;

  const float bpm = safeDiv(60.0f, med, 0.0f);
  if (!std::isfinite(bpm) || bpm < params_.bpmMin || bpm > params_.bpmMax) return;

  // Compute a simple confidence from coefficient of variation + sample count.
  double sum = 0.0;
  for (float v : intervals_) sum += (double)v;
  const double mean = sum / (double)intervals_.size();

  double var = 0.0;
  for (float v : intervals_) {
    const double d = (double)v - mean;
    var += d * d;
  }
  var /= std::max<std::size_t>(1, intervals_.size() - 1);

  const double stdDev = std::sqrt(var);
  const float cv = (float)(stdDev / std::max(1.0e-9, mean));

  // cv ~= 0 -> stable, cv >= 0.25 -> low confidence
  const float stability = clamp01(1.0f - cv * 4.0f);
  const float countBoost = clamp01((float)intervals_.size() / 8.0f);

  bpm_ = bpm;
  bpmConf01_ = clamp01(stability * countBoost);
}

} // namespace stellar::dsp
