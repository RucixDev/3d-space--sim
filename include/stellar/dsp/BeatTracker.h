#pragma once

#include <cstddef>
#include <vector>

namespace stellar::dsp {

// -----------------------------------------------------------------------------
// Beat / Onset tracking (spectral flux)
// -----------------------------------------------------------------------------
//
// A tiny beat/onset tracker intended for:
//  - debug visualization tooling (Audio Analyzer window)
//  - driving audio-reactive procedural visuals
//
// The tracker consumes a *magnitude spectrum* stream (0..Nyquist), e.g. the
// output of AudioAnalyzer::spectrumMag(). It computes a spectral-flux onset
// function, applies an adaptive threshold (EMA mean + EMA deviation), then
// performs simple peak picking with a refractory period.
//
// This is dependency-free and real-time friendly:
//  - no allocations in processSpectrum() after the first call
//  - deterministic behavior for identical inputs

struct BeatTrackerParams {
  // EMA time constants (seconds).
  float meanTauSec{0.50f};
  float devTauSec{0.50f};

  // Threshold = mean + sensitivity * deviation.
  float sensitivity{1.60f};

  // Minimum time between accepted peaks (seconds).
  float minBeatIntervalSec{0.25f};

  // Beat pulse response: on beat, beatPulse01 is set to 1 and then decays.
  // A linear decay is used: beatPulse01 -= dt / beatPulseDecaySec.
  float beatPulseDecaySec{0.25f};

  // BPM estimation window (using recent beat-to-beat intervals).
  float bpmMin{55.0f};
  float bpmMax{200.0f};
  int bpmHistory{12};
};

struct BeatTrackerOutput {
  // Raw onset (spectral flux) for the current frame.
  float onset{0.0f};

  // Adaptive threshold for onset.
  float threshold{0.0f};

  // True on the frame a beat is detected.
  bool beat{false};

  // A decaying pulse in [0,1], useful for UI and shader driving.
  float beatPulse01{0.0f};

  // Estimated tempo (0 if unknown).
  float bpm{0.0f};

  // Heuristic confidence for bpm in [0,1].
  float bpmConfidence01{0.0f};

  // Debug: running mean/deviation of onset.
  float mean{0.0f};
  float deviation{0.0f};
};

class BeatTracker {
public:
  BeatTracker() = default;
  explicit BeatTracker(const BeatTrackerParams& p) { setParams(p); }

  void setParams(const BeatTrackerParams& p);
  const BeatTrackerParams& params() const { return params_; }

  void reset();

  // Process a magnitude spectrum.
  //
  // mag: pointer to spectrum magnitudes (size = count)
  // count: number of bins (typically fftSize/2)
  // dtSec: time since previous call (seconds)
  BeatTrackerOutput processSpectrum(const float* mag, int count, float dtSec);

  const BeatTrackerOutput& last() const { return last_; }

private:
  BeatTrackerParams params_{};
  BeatTrackerOutput last_{};

  // Previous spectrum for flux.
  std::vector<float> prevMag_{};

  // Internal time (seconds).
  double timeSec_{0.0};

  // Peak picking state (3-sample local maxima test).
  float onsetPrev2_{0.0f};
  float onsetPrev1_{0.0f};
  float thrPrev2_{0.0f};
  float thrPrev1_{0.0f};

  double tPrev2_{0.0};
  double tPrev1_{0.0};

  // Adaptive statistics.
  float mean_{0.0f};
  float dev_{0.0f};

  // Beat timing.
  double lastBeatSec_{-1.0};

  // Recent beat-to-beat intervals (seconds).
  std::vector<float> intervals_{};

  // Cached tempo estimate.
  float bpm_{0.0f};
  float bpmConf01_{0.0f};

  // Beat pulse.
  float beatPulse01_{0.0f};

  void updateTempo_();
};

} // namespace stellar::dsp
