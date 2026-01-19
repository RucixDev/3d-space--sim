#include "test_harness.h"

#include "stellar/dsp/BeatTracker.h"

#include <vector>

int test_beat_tracker() {
  int failures = 0;

  using stellar::dsp::BeatTracker;
  using stellar::dsp::BeatTrackerParams;

  BeatTrackerParams p;
  p.meanTauSec = 0.25f;
  p.devTauSec = 0.25f;
  p.sensitivity = 1.20f;
  p.minBeatIntervalSec = 0.30f;
  p.beatPulseDecaySec = 0.25f;
  p.bpmMin = 90.0f;
  p.bpmMax = 150.0f;
  p.bpmHistory = 12;

  BeatTracker bt(p);

  const float dt = 1.0f / 60.0f;
  const int frames = 600; // 10 seconds
  const int bins = 128;

  std::vector<float> spec((std::size_t)bins, 0.05f);

  int beatCount = 0;
  for (int i = 0; i < frames; ++i) {
    const bool isBeatFrame = (i % 30) == 0; // 120 BPM
    const float energy = isBeatFrame ? 1.0f : 0.05f;

    for (int k = 0; k < bins; ++k) {
      // Add a slight tilt so bins are not perfectly uniform.
      spec[(std::size_t)k] = energy * (1.0f + 0.002f * (float)k);
    }

    const auto out = bt.processSpectrum(spec.data(), bins, dt);
    if (out.beat) ++beatCount;

    CHECK(out.beatPulse01 >= 0.0f);
    CHECK(out.beatPulse01 <= 1.0f);
  }

  // Expect close to 20 beats in 10 seconds. Allow a small warm-up tolerance.
  CHECK(beatCount >= 17);
  CHECK(beatCount <= 22);

  const auto& last = bt.last();
  CHECK(last.bpm > 115.0f);
  CHECK(last.bpm < 125.0f);
  CHECK(last.bpmConfidence01 >= 0.10f);
  CHECK(last.bpmConfidence01 <= 1.0f);

  return failures;
}
