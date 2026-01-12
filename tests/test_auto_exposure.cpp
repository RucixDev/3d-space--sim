#include "stellar/render/AutoExposure.h"

#include <cmath>
#include <iostream>

int test_auto_exposure() {
  int fails = 0;

  stellar::render::AutoExposureConfig cfg{};
  cfg.key = 0.18f;
  cfg.minExposure = 0.10f;
  cfg.maxExposure = 10.0f;
  cfg.speedUp = 4.0f;
  cfg.speedDown = 1.0f;

  // Target exposure sanity: when avg luminance == key, exposure should be ~1.
  {
    const float t = stellar::render::computeTargetExposure(0.18f, cfg);
    if (std::fabs(t - 1.0f) > 1e-4f) {
      std::cerr << "[test_auto_exposure] target exposure expected ~1, got " << t << "\n";
      ++fails;
    }
  }

  stellar::render::AutoExposureState st{};

  // First step should initialize state to the target (no smoothing).
  {
    const float logLum = std::log(0.18f);
    const float e = stellar::render::stepAutoExposure(st, 0.18f, logLum, 0.016f, cfg);
    if (!st.valid || std::fabs(e - 1.0f) > 1e-3f) {
      std::cerr << "[test_auto_exposure] init step failed: valid=" << st.valid << " e=" << e << "\n";
      ++fails;
    }
  }

  // Darker scene -> target exposure goes up, should move upward but not overshoot.
  {
    const float target = stellar::render::computeTargetExposure(0.09f, cfg); // ~2.0
    const float e0 = st.exposure;
    const float e1 = stellar::render::stepAutoExposure(st, 0.09f, std::log(0.09f), 0.10f, cfg);
    if (!(e1 > e0 && e1 < target + 1e-4f)) {
      std::cerr << "[test_auto_exposure] brighten step unexpected: e0=" << e0 << " e1=" << e1
                << " target=" << target << "\n";
      ++fails;
    }
  }

  // Brighter scene -> target exposure goes down, should move downward.
  {
    const float target = stellar::render::computeTargetExposure(0.36f, cfg); // ~0.5
    const float e0 = st.exposure;
    const float e1 = stellar::render::stepAutoExposure(st, 0.36f, std::log(0.36f), 0.10f, cfg);
    if (!(e1 < e0 && e1 > target - 1e-4f)) {
      std::cerr << "[test_auto_exposure] darken step unexpected: e0=" << e0 << " e1=" << e1
                << " target=" << target << "\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_auto_exposure] pass\n";
  return fails;
}
