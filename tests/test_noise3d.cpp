#include "stellar/proc/Noise.h"

#include <cmath>
#include <iostream>

int test_noise3d() {
  int fails = 0;

  const stellar::core::u64 seed = 1337;

  // Integer lattice consistency: smoothNoise3D at int coordinates should match
  // the underlying valueNoise3D exactly.
  {
    const int x = 3, y = -7, z = 11;
    const double a = stellar::proc::smoothNoise3D(seed, (double)x, (double)y, (double)z);
    const double b = stellar::proc::valueNoise3D(seed, x, y, z);
    if (a != b) {
      std::cerr << "[test_noise3d] smoothNoise3D != valueNoise3D at lattice point\n";
      ++fails;
    }
  }

  // Range checks + determinism at fractional coordinates.
  {
    const double x = 12.25, y = -3.75, z = 0.125;
    const double n0 = stellar::proc::smoothNoise3D(seed, x, y, z);
    const double n1 = stellar::proc::smoothNoise3D(seed, x, y, z);
    if (n0 != n1) {
      std::cerr << "[test_noise3d] smoothNoise3D not deterministic\n";
      ++fails;
    }
    if (!(n0 >= 0.0 && n0 <= 1.0)) {
      std::cerr << "[test_noise3d] smoothNoise3D out of range: " << n0 << "\n";
      ++fails;
    }
  }

  // fBm range: with the default gain=0.5 and octaves=5, the amplitude sum is 0.96875.
  {
    const int oct = 5;
    const double gain = 0.5;
    double amp = 0.5;
    double sum = 0.0;
    for (int i = 0; i < oct; ++i) {
      sum += amp;
      amp *= gain;
    }

    const double f = stellar::proc::fbm3D(seed, 0.33, 1.25, -2.75, oct, 2.0, gain);
    if (!(f >= 0.0 && f <= sum + 1e-12)) {
      std::cerr << "[test_noise3d] fbm3D out of expected range: " << f << " (sumAmp=" << sum << ")\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_noise3d] pass\n";
  return fails;
}
