#pragma once

#include "stellar/core/Types.h"

#include <algorithm>

namespace stellar::sim {

// -----------------------------------------------------------------------------
// FieldRepair
// -----------------------------------------------------------------------------
// Deterministic helper for in-flight "jury-rig" hull repairs using basic
// materials (Metals + Machinery).
//
// Design intent:
//  - Keep gameplay logic (when allowed, UI, consequences) in the game layer.
//  - Keep the math deterministic and testable in the core sim library.
//  - Provide a soft cap so full repairs still require station services.
//
// Units:
//  - hull is an abstract "hull points" scalar (same scale as stellar_game).
//  - inventory is in commodity units.
//  - dtRealSec is in real seconds (not sim-time).

struct FieldRepairParams {
  // Maximum hull fraction reachable by field repairs (0..1).
  // Example: 0.75 means you can patch up to 75% hull.
  double maxRepairFrac{0.75};

  // Repair rate (hull points per real second).
  double hullPerRealSec{2.5};

  // Material recipe (commodity units per hull point).
  // Defaults are intentionally *less efficient* than station services.
  double metalsPerHull{0.30};
  double machineryPerHull{0.22};

  // Heat added per hull point repaired (ThermalSystem "heat units").
  // Caller typically adds `heatAdded` to ThermalInputs::heatImpulse.
  double heatPerHull{4.0};

  // Numerical epsilon.
  double eps{1e-9};
};

struct FieldRepairInventory {
  double metals{0.0};
  double machinery{0.0};
};

struct FieldRepairQuote {
  bool ok{false};
  const char* reason{nullptr};

  double hullCurrent{0.0};
  double hullMax{0.0};

  double capHull{0.0};       // hull limit under field repairs
  double hullMissingToCap{0.0};
  double hullToRepair{0.0};

  // Material plan.
  double metalsNeeded{0.0};
  double machineryNeeded{0.0};

  // Time to complete hullToRepair at hullPerRealSec.
  double timeSec{0.0};

  bool limitedByCap{false};
  bool limitedByStock{false};
};

// Compute an actionable plan to repair toward the "field repair cap".
// Returns ok=false with a short reason if repair is impossible.
FieldRepairQuote quoteFieldRepairToCap(const FieldRepairInventory& inv,
                                       double hullCurrent,
                                       double hullMax,
                                       const FieldRepairParams& params = {});

struct FieldRepairStepResult {
  bool progressed{false};
  bool done{false};
  bool reachedCap{false};
  bool outOfSupplies{false};

  double capHull{0.0};
  double hullRepaired{0.0};

  double metalsTaken{0.0};
  double machineryTaken{0.0};

  double heatAdded{0.0};
};

// Apply a single repair step for dtRealSec.
//
// - Consumes from ioInv.
// - Increases ioHullCurrent.
// - Clamps repairs to `params.maxRepairFrac` of hullMax.
// - Returns heatAdded for this step (caller usually adds to ThermalInputs::heatImpulse).
FieldRepairStepResult stepFieldRepair(double dtRealSec,
                                      double hullMax,
                                      FieldRepairInventory& ioInv,
                                      double& ioHullCurrent,
                                      const FieldRepairParams& params = {});

} // namespace stellar::sim
