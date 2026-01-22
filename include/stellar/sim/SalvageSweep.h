#pragma once

#include "stellar/core/Types.h"
#include "stellar/econ/Commodity.h"
#include "stellar/math/Vec3.h"

#include <vector>

namespace stellar::sim {

// Candidate salvage pod for the SalvageSweepComputer.
//
// This structure is intentionally lightweight and headless; the caller is
// responsible for constructing a list of candidates from their world state.
struct SalvagePodCandidate {
  core::u64 id{0};
  econ::CommodityId commodity{econ::CommodityId::Food};
  double units{0.0};

  math::Vec3d posKm{0, 0, 0};
  math::Vec3d velKmS{0, 0, 0};

  // Optional: non-zero if the pod was spawned as part of a mission objective.
  core::u64 missionId{0};

  // Optional: estimated value. If <= 0, it is treated as unknown/zero.
  double estimatedValueCr{0.0};
};

struct SalvageSweepParams {
  // Normalization scales used to convert raw value/distance into 0..1-ish terms.
  double valueScaleCr{4500.0};
  double distScaleKm{65000.0};

  // Score weights:
  //   score = missionBonus + valueWeight*valueTerm - distWeight*distTerm
  double valueWeight{1.0};
  double distWeight{1.0};
  double missionBonus{0.65};

  // Switching behavior.
  // Prevents thrashing between similar candidates when values are close.
  double switchHysteresis{0.18};

  // Minimum time between target switches (seconds).
  double minHoldSeconds{1.25};

  // Skip duration when the user requests to skip a target (seconds).
  double skipSeconds{45.0};
};

struct SalvageSweepStatus {
  core::u64 targetId{0};
  int targetCandidateIndex{-1}; // index into the candidate list passed to update()
  double score{0.0};
  bool hasTarget{false};
};

class SalvageSweepComputer {
public:
  SalvageSweepComputer();

  const SalvageSweepParams& params() const { return params_; }
  void setParams(const SalvageSweepParams& p) { params_ = p; }

  void reset();

  // Skip a target id for params.skipSeconds. If id==0, skips the current target.
  void skip(core::u64 id, double timeSec);

  // Pick/update a target.
  // `timeSec` should be monotonically increasing (e.g. real time seconds).
  SalvageSweepStatus update(double timeSec,
                            const math::Vec3d& shipPosKm,
                            const std::vector<SalvagePodCandidate>& candidates);

private:
  SalvageSweepParams params_{};

  core::u64 targetId_{0};
  double targetScore_{0.0};
  double lastSwitchTimeSec_{-1.0e9};

  struct SkipEntry {
    core::u64 id{0};
    double untilTimeSec{0.0};
  };
  std::vector<SkipEntry> skip_;
};

} // namespace stellar::sim
