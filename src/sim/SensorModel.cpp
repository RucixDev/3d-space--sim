#include "stellar/sim/SensorModel.h"

#include <cmath>

namespace stellar::sim {

static inline double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

double computeSensorStrength01(double distKm,
                               double signature,
                               double sensorPower,
                               bool occluded,
                               const SensorParams& params) {
  distKm = std::max(distKm, params.minDistKm);
  signature = std::max(0.0, signature);
  sensorPower = std::max(0.0, sensorPower);

  const double r = std::max(params.halfRangeKm, 1e-6);
  double snr = sensorPower * signature * (r * r) / (distKm * distKm);
  if (occluded) {
    snr *= clamp01(params.occlusionAtten);
  }

  if (!std::isfinite(snr) || snr <= 0.0) return 0.0;

  const double strength = snr / (1.0 + snr);
  return clamp01(strength);
}

static inline double expSmoothingAlpha(double dtSec, double halfLifeSec) {
  dtSec = std::max(0.0, dtSec);
  halfLifeSec = std::max(1e-6, halfLifeSec);

  // For y += a*(x-y) with constant x, choosing:
  //   a = 1 - exp(-dt/tau), tau = halfLife/ln(2)
  // ensures y reaches half the gap after 'halfLife' seconds.
  const double tau = halfLifeSec / std::log(2.0);
  return 1.0 - std::exp(-dtSec / tau);
}

SensorTrackResult updateSensorTrack(SensorTrack& track,
                                    double dtSec,
                                    double measuredStrength01,
                                    const SensorTrackParams& params) {
  measuredStrength01 = clamp01(measuredStrength01);

  const double halfLife = (measuredStrength01 > track.strength01)
                            ? params.riseHalfLifeSec
                            : params.fallHalfLifeSec;

  const double a = expSmoothingAlpha(dtSec, halfLife);
  track.strength01 += a * (measuredStrength01 - track.strength01);
  track.strength01 = clamp01(track.strength01);

  const double idTh = clamp01(params.identifyThreshold);
  const double keepTh = clamp01(params.maintainIdentifyThreshold);

  if (!track.identified && track.strength01 >= idTh) {
    track.identified = true;
  } else if (track.identified && track.strength01 < keepTh) {
    track.identified = false;
  }

  SensorTrackResult out{};
  out.strength01 = track.strength01;
  out.visible = track.strength01 >= clamp01(params.ghostThreshold);
  out.identified = track.identified;
  return out;
}

} // namespace stellar::sim
