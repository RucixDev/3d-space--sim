#include "stellar/sim/ManeuverSequenceComputer.h"

#include <algorithm>

namespace stellar::sim {

const ManeuverPlan* ManeuverSequenceComputer::currentPlan() const {
  if (nodeIndex_ < 0 || nodeIndex_ >= (int)plans_.size()) return nullptr;
  return &plans_[(std::size_t)nodeIndex_];
}

void ManeuverSequenceComputer::disengage() {
  phase_ = ManeuverSequencePhase::Off;
  plans_.clear();
  nodeIndex_ = -1;
  nodeComputer_.disengage();
}

void ManeuverSequenceComputer::startNode_(const Ship& ship) {
  nodeComputer_.disengage();

  if (nodeIndex_ < 0 || nodeIndex_ >= (int)plans_.size()) {
    return;
  }

  nodeComputer_.engage(ship, plans_[(std::size_t)nodeIndex_]);
}

void ManeuverSequenceComputer::engage(const Ship& ship,
                                     const std::vector<ManeuverPlan>& plans,
                                     bool sortByTime) {
  plans_ = plans;

  if (sortByTime) {
    std::stable_sort(plans_.begin(), plans_.end(),
                     [](const ManeuverPlan& a, const ManeuverPlan& b) {
                       return a.nodeTimeDays < b.nodeTimeDays;
                     });
  }

  if (plans_.empty()) {
    disengage();
    return;
  }

  phase_ = ManeuverSequencePhase::Executing;
  nodeIndex_ = 0;
  startNode_(ship);

  // If the first node completes immediately (zero Δv), allow update() to advance
  // across any remaining nodes.
  if (plans_.size() == 1 && nodeComputer_.phase() == ManeuverComputerPhase::Complete) {
    phase_ = ManeuverSequencePhase::Complete;
  }
}

ManeuverSequenceOutput ManeuverSequenceComputer::update(const Ship& ship,
                                                       double nowTimeDays,
                                                       double dtSimSec,
                                                       const ManeuverComputerParams& nodeParams) {
  ManeuverSequenceOutput out{};
  out.phase = phase_;
  out.nodeIndex = nodeIndex_;
  out.nodeCount = (int)plans_.size();

  if (phase_ == ManeuverSequencePhase::Off) {
    return out;
  }

  if (phase_ == ManeuverSequencePhase::Complete) {
    out.finished = true;
    return out;
  }

  if (phase_ == ManeuverSequencePhase::Aborted) {
    out.finished = true;
    out.aborted = true;
    return out;
  }

  if (nodeIndex_ < 0 || nodeIndex_ >= (int)plans_.size()) {
    phase_ = ManeuverSequencePhase::Aborted;
    out.phase = phase_;
    out.finished = true;
    out.aborted = true;
    return out;
  }

  const auto nodeOut = nodeComputer_.update(ship, nowTimeDays, dtSimSec, nodeParams);
  out.node = nodeOut;
  out.input = nodeOut.input;

  if (nodeOut.phase == ManeuverComputerPhase::Aborted) {
    phase_ = ManeuverSequencePhase::Aborted;
    out.phase = phase_;
    out.finished = true;
    out.aborted = true;
    return out;
  }

  if (nodeOut.phase == ManeuverComputerPhase::Complete) {
    const int next = nodeIndex_ + 1;
    if (next < (int)plans_.size()) {
      nodeIndex_ = next;
      out.advanced = true;

      startNode_(ship);

      // The output of this frame corresponds to the completed node (typically no
      // thrust). The newly armed node will drive the next frame.
      out.phase = phase_;
      out.nodeIndex = nodeIndex_;
      out.nodeCount = (int)plans_.size();
      return out;
    }

    phase_ = ManeuverSequencePhase::Complete;
    out.phase = phase_;
    out.finished = true;
    return out;
  }

  out.phase = phase_;
  out.nodeIndex = nodeIndex_;
  out.nodeCount = (int)plans_.size();
  return out;
}

} // namespace stellar::sim
