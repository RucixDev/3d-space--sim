#pragma once

#include "stellar/sim/KinematicTrack.h"

#include <algorithm>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// Fire control helpers
// -----------------------------------------------------------------------------
//
// These helpers let gameplay/UI tie weapons assistance (lead indicators, soft aim assist)
// to sensor quality in a deterministic way.
//
// The radar HUD already maintains a kinematic track (alpha-beta filter) with a simple
// scalar uncertainty model (sigmaKm + ageSinceMeasSec). Fire control can consume that
// track to decide whether guidance is trustworthy and how strongly to apply it.

inline double fireControlQuality01(const KinematicTrack3d& tr,
                                  double maxAgeSec,
                                  double maxSigmaKm) {
  if (!tr.initialized) return 0.0;
  const double qAge = (maxAgeSec > 1e-9) ? std::clamp(1.0 - tr.ageSinceMeasSec / maxAgeSec, 0.0, 1.0) : 0.0;
  const double qSig = (maxSigmaKm > 1e-9) ? std::clamp(1.0 - tr.sigmaKm / maxSigmaKm, 0.0, 1.0) : 0.0;
  return std::min(qAge, qSig);
}

struct FireControlEstimate {
  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};
  double quality01{0.0};
  double sigmaKm{0.0};
  double ageSinceMeasSec{0.0};
  bool valid{false};
};

inline FireControlEstimate estimateFireControlFromTrack(const KinematicTrack3d& tr,
                                                       double maxAgeSec,
                                                       double maxSigmaKm) {
  FireControlEstimate out{};
  if (!tr.initialized) return out;

  out.posKm = tr.posKm;
  out.velKmS = tr.velKmS;
  out.sigmaKm = tr.sigmaKm;
  out.ageSinceMeasSec = tr.ageSinceMeasSec;
  out.quality01 = fireControlQuality01(tr, maxAgeSec, maxSigmaKm);
  out.valid = out.quality01 > 1e-4;
  return out;
}

} // namespace stellar::sim
