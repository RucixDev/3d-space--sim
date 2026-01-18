#include "test_harness.h"

#include "stellar/dsp/AudioAnalyzer.h"

#include <algorithm>
#include <cmath>
#include <vector>

int test_audio_analyzer() {
  int failures = 0;

  using stellar::dsp::AudioAnalyzer;
  using stellar::dsp::AudioAnalyzerParams;

  constexpr float sampleRate = 48000.0f;
  constexpr int fftN = 1024;

  // Pick an FFT bin exactly so the expected peak is deterministic.
  constexpr int targetBin = 23;
  const float freq = sampleRate * (float)targetBin / (float)fftN;

  std::vector<float> samples;
  samples.resize(4096);

  const float amp = 0.85f;
  for (int i = 0; i < (int)samples.size(); ++i) {
    const float t = (float)i / sampleRate;
    samples[(std::size_t)i] = amp * std::sin(2.0f * 3.14159265358979323846f * freq * t);
  }

  AudioAnalyzerParams p;
  p.fftSize = fftN;
  p.bands = 64;
  p.smoothing = 0.0f;
  p.minHz = 20.0f;
  p.maxHz = 20000.0f;
  p.floorDb = -80.0f;

  AudioAnalyzer az(p);
  az.analyzeWindow(samples.data(), (int)samples.size(), sampleRate);

  CHECK((int)az.waveform().size() == fftN);
  CHECK((int)az.spectrumMag().size() == fftN / 2);
  CHECK((int)az.bandsMag().size() == p.bands);
  CHECK((int)az.bandsDb01().size() == p.bands);

  // Peak should be close to the sine amplitude.
  CHECK(std::isfinite(az.peak()));
  CHECK(az.peak() > 0.70f);
  CHECK(az.peak() <= 1.0f);

  // Find the maximum spectrum bin.
  int maxIdx = 0;
  float maxVal = -1.0f;
  const auto& mag = az.spectrumMag();
  for (int k = 0; k < (int)mag.size(); ++k) {
    if (mag[(std::size_t)k] > maxVal) {
      maxVal = mag[(std::size_t)k];
      maxIdx = k;
    }
  }

  CHECK(std::isfinite(maxVal));
  CHECK(maxVal > 0.0f);
  CHECK(std::abs(maxIdx - targetBin) <= 1);

  // dB-normalized bands should be in [0,1].
  for (float v : az.bandsDb01()) {
    CHECK(std::isfinite(v));
    CHECK(v >= 0.0f);
    CHECK(v <= 1.0f);
  }

  return failures;
}
