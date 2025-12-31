#include "stellar/sim/Docking.h"
#include "stellar/sim/DockingComputer.h"

#include "stellar/sim/Ship.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool inRange(double v, double lo, double hi, double eps = 1e-9) {
  return v >= (lo - eps) && v <= (hi + eps);
}

int test_docking() {
  int fails = 0;

  // --- Station hull helper (box hull with tunnel cutout) ---
  {
    sim::Station st;
    st.radiusKm = 100.0;
    st.slotWidthKm = 60.0;
    st.slotHeightKm = 25.0;
    st.slotDepthKm = 120.0;
    st.approachLengthKm = 500.0;
    st.approachRadiusKm = 80.0;
    st.maxApproachSpeedKmS = 4.0;

    // Outside hull box.
    const math::Vec3d outside{st.radiusKm * 0.80, 0.0, 0.0}; // wx=0.70*R -> 70
    if (sim::insideStationHullExceptSlot(st, outside)) {
      std::cerr << "[test_docking] insideStationHullExceptSlot: expected outside box.\n";
      ++fails;
    }

    // Inside box and outside tunnel -> should collide.
    const double wx = st.radiusKm * 0.70;
    const double slotHalfW = st.slotWidthKm * 0.5;
    const math::Vec3d insideSolid{std::min(wx * 0.95, slotHalfW * 1.20), 0.0, 0.0};
    if (!sim::insideStationHullExceptSlot(st, insideSolid)) {
      std::cerr << "[test_docking] insideStationHullExceptSlot: expected inside hull solid region.\n";
      ++fails;
    }

    // Inside tunnel volume -> should NOT collide (cutout).
    const double zEntrance = st.radiusKm * 1.10;
    const double zMin = zEntrance - st.slotDepthKm;
    const math::Vec3d insideTunnel{0.0, 0.0, (zMin + zEntrance) * 0.5};
    if (sim::insideStationHullExceptSlot(st, insideTunnel)) {
      std::cerr << "[test_docking] insideStationHullExceptSlot: expected tunnel cutout.\n";
      ++fails;
    }
  }

  // --- Docking slot conditions helper ---
  {
    sim::Station st;
    st.radiusKm = 100.0;
    st.slotWidthKm = 60.0;
    st.slotHeightKm = 25.0;
    st.slotDepthKm = 120.0;
    st.maxApproachSpeedKmS = 4.0;

    const double zEntrance = st.radiusKm * 1.10;
    const double zMin = zEntrance - st.slotDepthKm;

    // Well inside the tunnel band.
    math::Vec3d relLocal{0.0, 0.0, zEntrance - 0.10 * st.radiusKm};
    if (!inRange(relLocal.z, zMin + 0.10 * st.radiusKm, zEntrance - 0.05 * st.radiusKm)) {
      std::cerr << "[test_docking] internal sanity: relLocal.z not in tunnel band.\n";
      ++fails;
    }

    const math::Vec3d velLocal{0.0, 0.0, -1.0};
    const math::Vec3d fwdLocal{0.0, 0.0, -1.0};
    const math::Vec3d upLocal{0.0, 1.0, 0.0};

    if (!sim::dockingSlotConditions(st, relLocal, velLocal, fwdLocal, upLocal, /*clearanceGranted=*/true)) {
      std::cerr << "[test_docking] dockingSlotConditions: expected PASS in nominal case.\n";
      ++fails;
    }

    if (sim::dockingSlotConditions(st, relLocal, velLocal, fwdLocal, upLocal, /*clearanceGranted=*/false)) {
      std::cerr << "[test_docking] dockingSlotConditions: expected FAIL without clearance.\n";
      ++fails;
    }

    const math::Vec3d tooFast{0.0, 0.0, -10.0};
    if (sim::dockingSlotConditions(st, relLocal, tooFast, fwdLocal, upLocal, true)) {
      std::cerr << "[test_docking] dockingSlotConditions: expected FAIL when too fast.\n";
      ++fails;
    }

    const math::Vec3d wrongFacing{0.0, 0.0, 1.0};
    if (sim::dockingSlotConditions(st, relLocal, velLocal, wrongFacing, upLocal, true)) {
      std::cerr << "[test_docking] dockingSlotConditions: expected FAIL when facing wrong way.\n";
      ++fails;
    }
  }

  // --- Docking computer: should be able to drive the ship into the slot ---
  {
    sim::Station st;
    st.radiusKm = 100.0;
    st.slotWidthKm = 60.0;
    st.slotHeightKm = 25.0;
    st.slotDepthKm = 120.0;
    st.commsRangeKm = 2000.0;
    st.approachLengthKm = 500.0;
    st.approachRadiusKm = 80.0;
    st.maxApproachSpeedKmS = 4.0;

    sim::Ship ship;
    ship.setPositionKm({0.0, 0.0, 800.0});
    ship.setVelocityKmS({0.0, 0.0, 0.0});
    ship.setOrientation({1, 0, 0, 0});
    ship.setAngularVelocityRadS({0.0, 0.0, 0.0});
    ship.setMaxLinearAccelKmS2(2.0);
    ship.setMaxLinearAccelBoostKmS2(4.0);
    ship.setMaxAngularAccelRadS2(10.0);

    sim::DockingComputer dc;
    dc.reset();

    const math::Vec3d stPos{0.0, 0.0, 0.0};
    const math::Quatd stQ{1, 0, 0, 0};
    const math::Vec3d stV{0.0, 0.0, 0.0};

    bool docked = false;
    const double dt = 0.25;
    for (int i = 0; i < 4000; ++i) {
      const auto out = dc.update(ship, st, stPos, stQ, stV, /*clearanceGranted=*/true);
      ship.step(dt, out.input);

      if (out.docked) {
        docked = true;
        break;
      }
    }

    if (!docked) {
      const math::Vec3d relLocal = stQ.conjugate().rotate(ship.positionKm() - stPos);
      std::cerr << "[test_docking] docking computer did not dock. "
                << "finalPos=(" << ship.positionKm().x << "," << ship.positionKm().y << "," << ship.positionKm().z << ") "
                << "relLocal=(" << relLocal.x << "," << relLocal.y << "," << relLocal.z << ") "
                << "phase=" << (int)dc.phase() << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_docking] PASS\n";
  }
  return fails;
}
