#include "stellar/sim/SensorModel.h"

#include <cmath>
#include <iostream>
#include <limits>

using namespace stellar;

static bool approxRel(double a, double b, double rel = 1e-9) {
  const double da = std::fabs(a - b);
  const double s = std::max({1.0, std::fabs(a), std::fabs(b)});
  return da <= rel * s;
}

int test_sensor_inversion() {
  int fails = 0;

  sim::SensorParams p{};
  p.halfRangeKm = 50000.0;
  p.minDistKm = 5.0;
  p.occlusionAtten = 0.10;

  const double sensorPower = 1.0;
  const double sig = 1.0;

  auto roundTrip = [&](double distKm, double occ01) {
    const double s01 = sim::computeSensorStrength01Occlusion01(distKm, sig, sensorPower, occ01, p);
    const double estKm = sim::estimateSensorRangeKmFromStrength01(s01, sig, sensorPower, occ01, p);

    if (!std::isfinite(estKm)) {
      std::cerr << "[test_sensor_inversion] expected finite estimate. distKm=" << distKm << " occ01=" << occ01
                << " strength=" << s01 << "\n";
      ++fails;
      return;
    }

    // The model + inversion are analytic; allow tiny numerical error.
    if (!approxRel(estKm, std::max(distKm, p.minDistKm), 1e-10)) {
      std::cerr << "[test_sensor_inversion] bad inversion. distKm=" << distKm << " occ01=" << occ01
                << " strength=" << s01 << " estKm=" << estKm << "\n";
      ++fails;
    }
  };

  // No occlusion
  roundTrip(1000.0, 0.0);
  roundTrip(25000.0, 0.0);
  roundTrip(75000.0, 0.0);

  // Partial occlusion
  roundTrip(1200.0, 0.25);
  roundTrip(25000.0, 0.60);
  roundTrip(90000.0, 0.90);

  // Edge cases
  {
    const double est = sim::estimateSensorRangeKmFromStrength01(0.0, 1.0, 1.0, 0.0, p);
    if (!std::isinf(est)) {
      std::cerr << "[test_sensor_inversion] expected +inf for strength=0. got=" << est << "\n";
      ++fails;
    }
  }

  {
    const double s01 = 0.999999999999;
    const double est = sim::estimateSensorRangeKmFromStrength01(s01, 1.0, 1.0, 0.0, p);
    if (!(std::isfinite(est) && est >= p.minDistKm - 1e-9)) {
      std::cerr << "[test_sensor_inversion] expected finite >= minDist for near-1 strength. got=" << est << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_sensor_inversion] PASS\n";
  }
  return fails;
}
