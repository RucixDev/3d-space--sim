#include "stellar/sim/SensorModel.h"

#include "test_harness.h"

#include <cmath>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_sensor_model() {
  int failures = 0;

  // --- Strength curve sanity: halfRange -> 0.5 ---
  {
    sim::SensorParams p{};
    p.halfRangeKm = 1000.0;

    const double s = sim::computeSensorStrength01(/*distKm=*/1000.0,
                                                  /*signature=*/1.0,
                                                  /*sensorPower=*/1.0,
                                                  /*occluded=*/false,
                                                  p);
    CHECK(approx(s, 0.5, 1e-12));
  }

  // --- Inverse-square behavior: double range => strength 0.2 (snr=0.25) ---
  {
    sim::SensorParams p{};
    p.halfRangeKm = 1000.0;

    const double s = sim::computeSensorStrength01(/*distKm=*/2000.0,
                                                  /*signature=*/1.0,
                                                  /*sensorPower=*/1.0,
                                                  /*occluded=*/false,
                                                  p);
    // snr=0.25 => 0.25/(1.25)=0.2
    CHECK(approx(s, 0.2, 1e-12));
  }

  // --- Occlusion attenuation reduces strength. ---
  {
    sim::SensorParams p{};
    p.halfRangeKm = 1000.0;
    p.occlusionAtten = 0.10;

    const double sVis = sim::computeSensorStrength01(/*distKm=*/1000.0, 1.0, 1.0, false, p);
    const double sOcc = sim::computeSensorStrength01(/*distKm=*/1000.0, 1.0, 1.0, true, p);
    CHECK(sOcc < sVis);
    CHECK(sOcc < 0.20);
  }

  // --- Track smoothing: dt==halfLife => move halfway. ---
  {
    sim::SensorTrack t{};
    sim::SensorTrackParams tp{};
    tp.riseHalfLifeSec = 0.30;

    const auto out = sim::updateSensorTrack(t, /*dtSec=*/0.30, /*measured=*/1.0, tp);
    CHECK(approx(out.strength01, 0.5, 1e-9));
    CHECK(out.visible); // should be above ghost threshold
  }

  // --- Identification hysteresis: identify, keep, drop. ---
  {
    sim::SensorTrack t{};
    sim::SensorTrackParams tp{};
    tp.identifyThreshold = 0.65;
    tp.maintainIdentifyThreshold = 0.45;

    t.strength01 = 0.66;
    t.identified = false;
    auto a = sim::updateSensorTrack(t, /*dtSec=*/0.0, /*measured=*/0.66, tp);
    CHECK(a.identified);

    t.strength01 = 0.46;
    t.identified = true;
    auto b = sim::updateSensorTrack(t, /*dtSec=*/0.0, /*measured=*/0.46, tp);
    CHECK(b.identified);

    t.strength01 = 0.44;
    t.identified = true;
    auto c = sim::updateSensorTrack(t, /*dtSec=*/0.0, /*measured=*/0.44, tp);
    CHECK(!c.identified);
  }

  return failures;
}
