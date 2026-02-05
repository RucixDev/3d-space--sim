#include "stellar/sim/SensorModel.h"

#include "stellar/math/Math.h"

#include <cmath>
#include <limits>

namespace stellar::sim {

static inline double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

double computeSensorStrength01(double distKm,
                               double signature,
                               double sensorPower,
                               bool occluded,
                               const SensorParams& params) {
  return computeSensorStrength01Occlusion01(distKm,
                                            signature,
                                            sensorPower,
                                            occluded ? 1.0 : 0.0,
                                            params);
}

double computeSensorStrength01Occlusion01(double distKm,
                                          double signature,
                                          double sensorPower,
                                          double occlusion01,
                                          const SensorParams& params) {
  distKm = std::max(distKm, params.minDistKm);
  signature = std::max(0.0, signature);
  sensorPower = std::max(0.0, sensorPower);

  const double r = std::max(params.halfRangeKm, 1e-6);
  double snr = sensorPower * signature * (r * r) / (distKm * distKm);

  occlusion01 = clamp01(occlusion01);
  const double occAtten = clamp01(params.occlusionAtten);
  const double atten = (1.0 - occlusion01) + occlusion01 * occAtten;
  snr *= atten;

  if (!std::isfinite(snr) || snr <= 0.0) return 0.0;

  const double strength = snr / (1.0 + snr);
  return clamp01(strength);
}


double estimateSensorRangeKmFromStrength01(double strength01,
                                          double signature,
                                          double sensorPower,
                                          double occlusion01,
                                          const SensorParams& params) {
  strength01 = clamp01(strength01);
  signature = std::max(0.0, signature);
  sensorPower = std::max(0.0, sensorPower);

  if (!std::isfinite(strength01) || strength01 <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }

  // Avoid division by zero as strength approaches 1.
  if (strength01 >= 1.0) strength01 = 1.0 - 1.0e-12;

  const double snr = strength01 / std::max(1.0e-12, 1.0 - strength01);

  occlusion01 = clamp01(occlusion01);
  const double occAtten = clamp01(params.occlusionAtten);
  const double atten = (1.0 - occlusion01) + occlusion01 * occAtten;

  const double gain = sensorPower * signature * atten;
  if (!std::isfinite(gain) || gain <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }

  const double r = std::max(1.0e-6, params.halfRangeKm);
  double distKm = r * std::sqrt(std::max(0.0, gain) / std::max(1.0e-18, snr));
  distKm = std::max(distKm, params.minDistKm);

  if (!std::isfinite(distKm)) {
    return std::numeric_limits<double>::infinity();
  }
  return distKm;
}

SensorTrackResult updateSensorTrack(SensorTrack& track,
                                    double dtSec,
                                    double measuredStrength01,
                                    const SensorTrackParams& params) {
  measuredStrength01 = clamp01(measuredStrength01);

  const double halfLife = (measuredStrength01 > track.strength01)
                            ? params.riseHalfLifeSec
                            : params.fallHalfLifeSec;

  const double a = math::halfLifeAlpha(dtSec, halfLife);
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
