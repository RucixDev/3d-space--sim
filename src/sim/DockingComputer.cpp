#include "stellar/sim/DockingComputer.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

DockingComputerResult DockingComputer::update(const Ship& ship,
                                             const Station& st,
                                             const math::Vec3d& stPosKm,
                                             const math::Quatd& stOrient,
                                             const math::Vec3d& stVelKmS,
                                             bool clearanceGranted,
                                             const DockingComputerParams& params) {
  DockingComputerResult out{};

  // Geometry in station-local space (+Z points out of the slot).
  const double zEntrance = st.radiusKm * 1.10;
  const double zMin = zEntrance - st.slotDepthKm;
  const double zCorridorEnd = zEntrance + st.approachLengthKm;
  const double zStage = zCorridorEnd + st.radiusKm * params.stageOffsetRadiusFrac;

  const math::Vec3d relWorldKm = ship.positionKm() - stPosKm;
  const math::Vec3d relLocal = stOrient.conjugate().rotate(relWorldKm);
  const double lateralDist = math::Vec3d{relLocal.x, relLocal.y, 0.0}.length();

  const double slotHalfW = st.slotWidthKm * 0.5;
  const double slotHalfH = st.slotHeightKm * 0.5;

  // Desired points for each phase in station-local space.
  const double zEntryHold = zEntrance + std::max(st.radiusKm * 0.35, st.slotWidthKm * 0.35);
  const double zAlignOk = zEntryHold + st.radiusKm * params.alignZExtraFrac;

  // If clearance isn't available yet, hold at staging point.
  if (!clearanceGranted) {
    phase_ = DockingComputerPhase::Staging;
    out.holdingForClearance = true;
  }

  // Phase switching:
  //  - Staging: hold at corridor end until centered, then transition to AlignEntry.
  //  - AlignEntry: align to a point just outside the slot. When centered + close, transition to SlotRun.
  //  - SlotRun: commit through the slot tunnel. If drifting close to the frame, fall back to AlignEntry.
  if (phase_ == DockingComputerPhase::Staging) {
    if (clearanceGranted &&
        // NOTE: we stage *beyond* the corridor end. Use the stage point (not the
        // corridor end) as the longitudinal threshold, otherwise some station
        // geometries can get stuck in Staging forever.
        relLocal.z < zStage * params.corridorEnterZFrac &&
        lateralDist < st.approachRadiusKm * params.corridorEnterFrac) {
      phase_ = DockingComputerPhase::AlignEntry;
    }
  } else if (phase_ == DockingComputerPhase::AlignEntry) {
    if (!clearanceGranted) {
      phase_ = DockingComputerPhase::Staging;
    } else if (std::abs(relLocal.x) < slotHalfW * params.alignSlotFrac &&
               std::abs(relLocal.y) < slotHalfH * params.alignSlotFrac &&
               relLocal.z < zAlignOk) {
      phase_ = DockingComputerPhase::SlotRun;
    } else if (relLocal.z > zStage ||
               lateralDist > st.approachRadiusKm * params.corridorResetFrac) {
      phase_ = DockingComputerPhase::Staging;
    }
  } else {
    // SlotRun
    if (!clearanceGranted) {
      phase_ = DockingComputerPhase::Staging;
    } else if (std::abs(relLocal.x) > slotHalfW * params.driftReAlignFrac ||
               std::abs(relLocal.y) > slotHalfH * params.driftReAlignFrac) {
      phase_ = DockingComputerPhase::AlignEntry;
    } else if (relLocal.z > zEntryHold + st.radiusKm * params.driftReAlignZExtraFrac) {
      phase_ = DockingComputerPhase::AlignEntry;
    }
  }

  double zDockAim = zMin + st.radiusKm * 0.25;
  // Clamp into the docking tunnel band used by dockingSlotConditions.
  zDockAim = std::clamp(zDockAim,
                        zMin + st.radiusKm * 0.15,
                        zEntrance - st.radiusKm * 0.15);

  math::Vec3d desiredLocal{0, 0, zStage};
  if (phase_ == DockingComputerPhase::AlignEntry) desiredLocal = {0, 0, zEntryHold};
  else if (phase_ == DockingComputerPhase::SlotRun) desiredLocal = {0, 0, zDockAim};

  const math::Vec3d desiredPointKm = stPosKm + stOrient.rotate(desiredLocal);
  out.desiredPointKm = desiredPointKm;

  // Speed scaling: be more cautious when not centered.
  double speedScale = 1.0;
  if (phase_ == DockingComputerPhase::AlignEntry) {
    const double frac = std::clamp(lateralDist / std::max(1.0, st.approachRadiusKm), 0.0, 1.0);
    speedScale = (1.0 - params.alignLateralSpeedPenalty * frac);
  } else if (phase_ == DockingComputerPhase::SlotRun) {
    const double fracX = std::clamp(std::abs(relLocal.x) / std::max(1.0, slotHalfW), 0.0, 1.0);
    const double fracY = std::clamp(std::abs(relLocal.y) / std::max(1.0, slotHalfH), 0.0, 1.0);
    const double frac = std::max(fracX, fracY);
    speedScale = (1.0 - params.slotLateralSpeedPenalty * frac);
  }

  // Phase-specific max speeds (relative to station).
  double baseMaxV = std::max(st.maxApproachSpeedKmS * params.stageSpeedMult, params.stageSpeedMin);
  if (phase_ == DockingComputerPhase::AlignEntry) {
    baseMaxV = std::max(st.maxApproachSpeedKmS * params.alignSpeedMult, params.alignSpeedMin);
  }
  if (phase_ == DockingComputerPhase::SlotRun) {
    baseMaxV = std::max(st.maxApproachSpeedKmS * params.slotSpeedMult, params.slotSpeedMin);
  }

  FlightControlParams fp{};
  fp.maxSpeedKmS = baseMaxV * speedScale;
  fp.speedGain = (phase_ == DockingComputerPhase::SlotRun ? params.slotSpeedGain : params.stageSpeedGain) * speedScale;
  fp.velGain = (phase_ == DockingComputerPhase::SlotRun ? params.slotVelGain : params.stageVelGain);
  fp.desiredDistKm = 0.0;
  fp.allowBoost = false;
  fp.dampers = true;

  AttitudeControlParams ap{};
  ap.faceGain = params.faceGain;
  ap.rollGain = params.rollGain;
  ap.alignUp = true;

  // Align forward into station (towards -axisOut) and keep ship "up" aligned to station +Y.
  const math::Vec3d axisOut = stOrient.rotate({0, 0, 1});
  const math::Vec3d desiredFwdWorld = -axisOut;
  const math::Vec3d desiredUpWorld = stOrient.rotate({0, 1, 0});

  const auto fc = approachTarget(ship,
                                 desiredPointKm,
                                 stVelKmS,
                                 fp,
                                 ap,
                                 desiredFwdWorld,
                                 &desiredUpWorld);

  out.input = fc.input;
  out.phase = phase_;

  // Auto-dock: once we are inside the slot tunnel under speed limit + aligned.
  if (phase_ == DockingComputerPhase::SlotRun && clearanceGranted) {
    const math::Vec3d relVelWorld = ship.velocityKmS() - stVelKmS;
    const math::Vec3d relVelLocal = stOrient.conjugate().rotate(relVelWorld);
    const math::Vec3d fwdLocal = stOrient.conjugate().rotate(ship.forward());
    const math::Vec3d upLocal  = stOrient.conjugate().rotate(ship.up());

    if (dockingSlotConditions(st, relLocal, relVelLocal, fwdLocal, upLocal, true)) {
      out.docked = true;
    }
  }

  return out;
}

} // namespace stellar::sim
