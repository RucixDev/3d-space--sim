#pragma once

#include "stellar/core/Types.h"
#include "stellar/sim/ManeuverComputer.h"

#include <vector>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// ManeuverSequenceComputer — execute multiple maneuver nodes in order
// -----------------------------------------------------------------------------
//
// This is a thin wrapper around ManeuverComputer that allows chaining multiple
// maneuver plans (e.g. a Lambert rendezvous with an arrival "match velocity"
// burn).
//
// Design goals:
//  - Deterministic + headless (usable in unit tests and tools).
//  - Minimal policy: nodes are executed in time order; if a node aborts, the
//    whole sequence aborts.
//  - Simple integration: output is a ShipInput override like ManeuverComputer.

enum class ManeuverSequencePhase : core::u8 {
  Off = 0,
  Executing = 1,
  Complete = 2,
  Aborted = 3,
};

struct ManeuverSequenceOutput {
  ShipInput input{};

  ManeuverSequencePhase phase{ManeuverSequencePhase::Off};

  // 0-based index of the currently armed node. -1 when Off.
  int nodeIndex{-1};
  int nodeCount{0};

  // Output from the underlying ManeuverComputer for the current node.
  ManeuverComputerOutput node{};

  // True if the *sequence* completed/aborted during this update.
  bool finished{false};
  bool aborted{false};

  // True if the sequence advanced to the next node during this update.
  bool advanced{false};
};

class ManeuverSequenceComputer {
public:
  // Reset to Off and clear any stored plan.
  void disengage();

  // Engage with a set of maneuver plans.
  //
  // If sortByTime is true, plans are executed in ascending nodeTimeDays.
  // If the vector is empty, the computer remains Off.
  void engage(const Ship& ship,
              const std::vector<ManeuverPlan>& plans,
              bool sortByTime = true);

  bool active() const { return phase_ == ManeuverSequencePhase::Executing; }
  ManeuverSequencePhase phase() const { return phase_; }

  int nodeIndex() const { return nodeIndex_; }
  int nodeCount() const { return (int)plans_.size(); }

  const std::vector<ManeuverPlan>& plans() const { return plans_; }
  const ManeuverPlan* currentPlan() const;

  // Step the sequence and return recommended ship inputs.
  ManeuverSequenceOutput update(const Ship& ship,
                                double nowTimeDays,
                                double dtSimSec,
                                const ManeuverComputerParams& nodeParams = {});

private:
  void startNode_(const Ship& ship);

  std::vector<ManeuverPlan> plans_{};
  int nodeIndex_{-1};

  ManeuverSequencePhase phase_{ManeuverSequencePhase::Off};

  ManeuverComputer nodeComputer_{};
};

} // namespace stellar::sim
