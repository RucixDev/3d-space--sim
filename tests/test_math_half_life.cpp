#include "stellar/math/Math.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approxEq(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

int test_math_half_life() {
  int fails = 0;

  // dt==0 should be a no-op: factor=1, alpha=0.
  {
    const double f0 = math::halfLifeDecayFactor(0.0, 10.0);
    if (!approxEq(f0, 1.0)) {
      std::cerr << "[test_math_half_life] expected dt=0 decayFactor=1 got=" << f0 << "\n";
      ++fails;
    }

    const double a0 = math::halfLifeAlpha(0.0, 10.0);
    if (!approxEq(a0, 0.0)) {
      std::cerr << "[test_math_half_life] expected dt=0 alpha=0 got=" << a0 << "\n";
      ++fails;
    }
  }

  // dt==halfLife => half remaining, half moved.
  {
    const double fHalf = math::halfLifeDecayFactor(5.0, 5.0);
    if (fHalf < 0.49 || fHalf > 0.51) {
      std::cerr << "[test_math_half_life] expected dt==halfLife decayFactor~0.5 got=" << fHalf << "\n";
      ++fails;
    }

    const double aHalf = math::halfLifeAlpha(5.0, 5.0);
    if (aHalf < 0.49 || aHalf > 0.51) {
      std::cerr << "[test_math_half_life] expected dt==halfLife alpha~0.5 got=" << aHalf << "\n";
      ++fails;
    }
  }

  // halfLife<=0 => instant (for dt>0): factor=0, alpha=1.
  {
    const double fDead = math::halfLifeDecayFactor(1.0, 0.0);
    if (!approxEq(fDead, 0.0)) {
      std::cerr << "[test_math_half_life] expected halfLife<=0 decayFactor=0 got=" << fDead << "\n";
      ++fails;
    }

    const double aDead = math::halfLifeAlpha(1.0, 0.0);
    if (!approxEq(aDead, 1.0)) {
      std::cerr << "[test_math_half_life] expected halfLife<=0 alpha=1 got=" << aDead << "\n";
      ++fails;
    }
  }

  // dt<0 should behave like dt=0 (defensive).
  {
    const double fNeg = math::halfLifeDecayFactor(-1.0, 10.0);
    if (!approxEq(fNeg, 1.0)) {
      std::cerr << "[test_math_half_life] expected dt<0 decayFactor=1 got=" << fNeg << "\n";
      ++fails;
    }

    const double aNeg = math::halfLifeAlpha(-1.0, 10.0);
    if (!approxEq(aNeg, 0.0)) {
      std::cerr << "[test_math_half_life] expected dt<0 alpha=0 got=" << aNeg << "\n";
      ++fails;
    }
  }

  return fails;
}
