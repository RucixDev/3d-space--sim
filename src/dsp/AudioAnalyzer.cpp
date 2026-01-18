#include "stellar/dsp/AudioAnalyzer.h"

#include <algorithm>
#include <cmath>

namespace stellar::dsp {
namespace {

constexpr float kPi = 3.14159265358979323846f;

static bool isPow2(int v) {
  return v > 0 && (v & (v - 1)) == 0;
}

static int nextPow2(int v) {
  int p = 1;
  while (p < v) p <<= 1;
  return p;
}

// In-place radix-2 decimation-in-time FFT.
//
// This is intentionally small and dependency-free (no FFTW/KISSFFT) because the
// project targets lightweight builds and headless unit tests.
static void fftRadix2(std::vector<std::complex<float>>& a) {
  const int n = (int)a.size();
  if (n <= 1) return;

  // Bit-reversal permutation.
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::swap(a[(std::size_t)i], a[(std::size_t)j]);
    }
  }

  for (int len = 2; len <= n; len <<= 1) {
    const float ang = -2.0f * kPi / (float)len;
    const std::complex<float> wLen(std::cos(ang), std::sin(ang));

    for (int i = 0; i < n; i += len) {
      std::complex<float> w(1.0f, 0.0f);
      const int half = len >> 1;
      for (int j = 0; j < half; ++j) {
        const std::complex<float> u = a[(std::size_t)(i + j)];
        const std::complex<float> v = a[(std::size_t)(i + j + half)] * w;
        a[(std::size_t)(i + j)] = u + v;
        a[(std::size_t)(i + j + half)] = u - v;
        w *= wLen;
      }
    }
  }
}

static float clamp01(float v) {
  return std::clamp(v, 0.0f, 1.0f);
}

static float clamp11(float v) {
  return std::clamp(v, -1.0f, 1.0f);
}

} // namespace

void AudioAnalyzer::setParams(const AudioAnalyzerParams& params) {
  AudioAnalyzerParams p = params;

  // Sanitize.
  p.fftSize = std::clamp(p.fftSize, 64, 16384);
  if (!isPow2(p.fftSize)) {
    p.fftSize = nextPow2(p.fftSize);
    p.fftSize = std::clamp(p.fftSize, 64, 16384);
  }

  p.bands = std::clamp(p.bands, 8, 512);
  p.minHz = std::max(1.0f, p.minHz);
  p.maxHz = std::max(p.minHz, p.maxHz);
  p.smoothing = std::clamp(p.smoothing, 0.0f, 0.99f);
  p.floorDb = std::min(p.floorDb, -6.0f); // keep a sensible range

  const bool sizeChanged = (p.fftSize != params_.fftSize) || (p.bands != params_.bands);

  params_ = p;
  if (sizeChanged || window_.empty()) {
    rebuild();
  }
}

void AudioAnalyzer::rebuild() {
  const int n = params_.fftSize;
  const int bands = params_.bands;

  window_.assign((std::size_t)n, 1.0f);
  waveform_.assign((std::size_t)n, 0.0f);
  fftBuf_.assign((std::size_t)n, std::complex<float>(0.0f, 0.0f));
  spectrumMag_.assign((std::size_t)(n / 2), 0.0f);

  bandsMag_.assign((std::size_t)bands, 0.0f);
  bandsDb01_.assign((std::size_t)bands, 0.0f);
  bandsTmp_.assign((std::size_t)bands, 0.0f);

  // Hann window (raised cosine).
  if (n > 1) {
    for (int i = 0; i < n; ++i) {
      const float t = (float)i / (float)(n - 1);
      window_[(std::size_t)i] = 0.5f - 0.5f * std::cos(2.0f * kPi * t);
    }
  }
}

void AudioAnalyzer::analyzeWindow(const float* samples, int count, float sampleRate) {
  if (params_.fftSize <= 0) return;
  if (sampleRate <= 0.0f) return;

  const int n = params_.fftSize;

  if ((int)waveform_.size() != n || (int)window_.size() != n || (int)fftBuf_.size() != n) {
    rebuild();
  }

  // Copy last N samples (pad leading zeros if necessary).
  std::fill(waveform_.begin(), waveform_.end(), 0.0f);
  if (samples && count > 0) {
    const int copyCount = std::min(n, count);
    const int srcStart = count - copyCount;
    const int dstStart = n - copyCount;
    for (int i = 0; i < copyCount; ++i) {
      waveform_[(std::size_t)(dstStart + i)] = clamp11(samples[srcStart + i]);
    }
  }

  // Time-domain stats.
  double sumSq = 0.0;
  float peak = 0.0f;
  for (int i = 0; i < n; ++i) {
    const float s = waveform_[(std::size_t)i];
    sumSq += (double)s * (double)s;
    peak = std::max(peak, std::abs(s));
  }
  rms_ = (n > 0) ? (float)std::sqrt(sumSq / (double)n) : 0.0f;
  peak_ = peak;

  // Window and FFT.
  for (int i = 0; i < n; ++i) {
    const float w = waveform_[(std::size_t)i] * window_[(std::size_t)i];
    fftBuf_[(std::size_t)i] = std::complex<float>(w, 0.0f);
  }

  fftRadix2(fftBuf_);

  // Magnitude spectrum (0..Nyquist).
  const int half = n / 2;
  if ((int)spectrumMag_.size() != half) spectrumMag_.assign((std::size_t)half, 0.0f);

  const float norm = 1.0f / (float)n;
  for (int k = 0; k < half; ++k) {
    spectrumMag_[(std::size_t)k] = std::abs(fftBuf_[(std::size_t)k]) * norm;
  }

  // Log-frequency bands (for plotting / audio-reactive uses).
  const float nyquist = 0.5f * sampleRate;
  float minHz = std::clamp(params_.minHz, 1.0f, nyquist);
  float maxHz = std::clamp(params_.maxHz, minHz, nyquist);

  const int bands = params_.bands;
  if ((int)bandsMag_.size() != bands) {
    bandsMag_.assign((std::size_t)bands, 0.0f);
    bandsDb01_.assign((std::size_t)bands, 0.0f);
    bandsTmp_.assign((std::size_t)bands, 0.0f);
  }

  const float logMin = std::log(minHz);
  const float logMax = std::log(maxHz);
  const float invBands = (bands > 0) ? (1.0f / (float)bands) : 0.0f;

  for (int b = 0; b < bands; ++b) {
    const float t0 = (float)b * invBands;
    const float t1 = (float)(b + 1) * invBands;

    const float f0 = std::exp(logMin + (logMax - logMin) * t0);
    const float f1 = std::exp(logMin + (logMax - logMin) * t1);

    int k0 = (int)std::floor(f0 * (float)n / sampleRate);
    int k1 = (int)std::ceil(f1 * (float)n / sampleRate);
    k0 = std::clamp(k0, 0, half - 1);
    k1 = std::clamp(k1, k0 + 1, half);

    float acc = 0.0f;
    for (int k = k0; k < k1; ++k) {
      acc += spectrumMag_[(std::size_t)k];
    }
    const float mag = acc / (float)(k1 - k0);
    bandsTmp_[(std::size_t)b] = mag;

    const float a = params_.smoothing;
    bandsMag_[(std::size_t)b] = a * bandsMag_[(std::size_t)b] + (1.0f - a) * mag;

    const float eps = 1.0e-6f;
    const float db = 20.0f * std::log10(bandsMag_[(std::size_t)b] + eps);
    const float floorDb = params_.floorDb;
    const float clampedDb = std::clamp(db, floorDb, 0.0f);
    const float db01 = (clampedDb - floorDb) / (0.0f - floorDb);
    bandsDb01_[(std::size_t)b] = clamp01(db01);
  }
}

} // namespace stellar::dsp
