#pragma once

#include "AudioEngine.h"
#include "stellar/dsp/AudioAnalyzer.h"
#include "stellar/dsp/BeatTracker.h"

#include <functional>
#include <string>
#include <vector>

namespace stellar::game {

namespace dsp = ::stellar::dsp;

// Audio Analyzer / Oscilloscope window.
//
// Provides:
//  - waveform (time domain) view
//  - log-banded spectrum (FFT + Hann window)
//  - quick SFX preview buttons
//  - export of the recent capture to a mono 16-bit WAV

struct AudioAnalyzerWindowState {
  bool open{false};

  bool pause{false};
  bool showWaveform{true};
  bool showSpectrum{true};

  // How many samples to fetch from the engine capture ring for analysis/export.
  int captureSamples{8192};

  // Analyzer parameters exposed to UI.
  dsp::AudioAnalyzerParams params{};
  dsp::AudioAnalyzer analyzer{params};

  // Beat / onset tracking (spectral flux).
  bool beatEnabled{true};
  bool beatPrevEnabled{true};
  dsp::BeatTrackerParams beatParams{};
  dsp::BeatTracker beat{beatParams};

  bool beatShowPlot{true};
  bool beatShowDebug{false};

  // Onset plot history.
  int onsetHistoryMax{420};
  std::vector<float> onsetHistory{};
  std::vector<float> thresholdHistory{};

  // Display options
  bool spectrumUseDb{true};
  bool spectrumLogBands{true};

  // Export path (relative to working directory by default).
  char exportPath[256]{"captures/audio_capture.wav"};

  // Scratch buffer for pulled audio.
  std::vector<float> capture{};
};

using ToastFn = std::function<void(const std::string& msg, double ttlSec)>;

void drawAudioAnalyzerWindow(AudioAnalyzerWindowState& st,
                             AudioEngine& audio,
                             const ToastFn& toast);

} // namespace stellar::game
