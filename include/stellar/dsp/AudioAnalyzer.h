#pragma once

#include <complex>
#include <vector>

namespace stellar::dsp {

// Lightweight real-time friendly audio analyzer.
//
// Intended uses:
//  - in-game visualization (oscilloscope / spectrum)
//  - driving audio-reactive procedural visuals
//
// This analyzer operates on a provided mono sample window.
// It does not own any audio devices and is safe for headless builds.
struct AudioAnalyzerParams {
  // FFT size (must be a power of two). Larger sizes have better frequency
  // resolution but more latency.
  int fftSize{1024};

  // Number of log-frequency bands to aggregate for UI plots.
  int bands{64};

  // Log-band range in Hz (clamped to [1, Nyquist]).
  float minHz{20.0f};
  float maxHz{20000.0f};

  // Exponential smoothing for bands (0 = none, 0.95 = very smooth).
  float smoothing{0.80f};

  // Floor used when mapping band magnitudes to 0..1 dB values.
  // (bandsDb01 uses 0dB as "full" and floorDb as 0.)
  float floorDb{-80.0f};
};

class AudioAnalyzer {
public:
  AudioAnalyzer() = default;
  explicit AudioAnalyzer(const AudioAnalyzerParams& params) { setParams(params); }

  void setParams(const AudioAnalyzerParams& params);
  const AudioAnalyzerParams& params() const { return params_; }

  // Analyze the last params().fftSize samples from `samples`.
  //
  // samples: pointer to mono float samples (typical range [-1, 1])
  // count: number of samples available
  // sampleRate: sample rate in Hz
  void analyzeWindow(const float* samples, int count, float sampleRate);

  // Results (updated by analyzeWindow).
  const std::vector<float>& waveform() const { return waveform_; }        // size = fftSize
  const std::vector<float>& spectrumMag() const { return spectrumMag_; }  // size = fftSize/2
  const std::vector<float>& bandsMag() const { return bandsMag_; }        // size = bands
  const std::vector<float>& bandsDb01() const { return bandsDb01_; }      // size = bands, 0..1

  float rms() const { return rms_; }
  float peak() const { return peak_; }

private:
  AudioAnalyzerParams params_{};

  std::vector<float> window_;
  std::vector<float> waveform_;
  std::vector<std::complex<float>> fftBuf_;
  std::vector<float> spectrumMag_;

  std::vector<float> bandsMag_;
  std::vector<float> bandsDb01_;
  std::vector<float> bandsTmp_;

  float rms_{0.0f};
  float peak_{0.0f};

  void rebuild();
};

} // namespace stellar::dsp
