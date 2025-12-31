#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"
#include "stellar/sim/Docking.h"
#include "stellar/sim/FlightController.h"
#include "stellar/sim/Ship.h"
#include "stellar/sim/System.h"

namespace stellar::sim {

enum class DockingComputerPhase : core::u8 {
  Staging   = 0,
  AlignEntry = 1,
  SlotRun   = 2,
};

// Parameters for the docking computer approach logic.
// Values are tuned to match the feel of the SDL prototype while keeping the
// implementation headless and testable.
struct DockingComputerParams {
  // --- Transition thresholds (fractions of station-defined geometry) ---
  double corridorEnterFrac{0.65};
  double corridorEnterZFrac{1.05};

  double alignSlotFrac{0.75};
  double alignZExtraFrac{0.25};

  double driftReAlignFrac{0.95};
  double driftReAlignZExtraFrac{0.50};

  double corridorResetFrac{1.25};

  // Stage point offset beyond the corridor end.
  double stageOffsetRadiusFrac{0.80};

  // --- Speed shaping ---
  double alignLateralSpeedPenalty{0.55};
  double slotLateralSpeedPenalty{0.70};

  double stageSpeedMult{2.80};
  double alignSpeedMult{1.35};
  double slotSpeedMult{0.75};

  double stageSpeedMin{0.45};
  double alignSpeedMin{0.28};
  double slotSpeedMin{0.12};

  // Controller gains used by sim::approachTarget.
  // Higher than the initial prototype value to ensure the approach converges
  // in a bounded time (the headless docking test runs a finite horizon).
  double stageSpeedGain{0.008};
  double slotSpeedGain{0.006};
  double stageVelGain{1.35};
  double slotVelGain{1.65};

  double faceGain{1.80};
  double rollGain{1.60};
};

struct DockingComputerResult {
  ShipInput input{};
  DockingComputerPhase phase{DockingComputerPhase::Staging};

  bool holdingForClearance{false};
  bool docked{false};

  // Optional debug value for tools/UI.
  math::Vec3d desiredPointKm{0,0,0};
};

class DockingComputer {
public:
  void reset() { phase_ = DockingComputerPhase::Staging; }
  DockingComputerPhase phase() const { return phase_; }

  DockingComputerResult update(const Ship& ship,
                               const Station& st,
                               const math::Vec3d& stPosKm,
                               const math::Quatd& stOrient,
                               const math::Vec3d& stVelKmS,
                               bool clearanceGranted,
                               const DockingComputerParams& params = {});

private:
  DockingComputerPhase phase_{DockingComputerPhase::Staging};
};

} // namespace stellar::sim
