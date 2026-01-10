#include "stellar/sim/SupercruiseComputer.h"

#include "stellar/sim/Ship.h"

#include <iostream>

using namespace stellar;

int test_supercruise() {
  int fails = 0;

  // A very lightweight "approach" simulation:
  // Start far away at max speed and ensure nav-assist can decelerate early enough
  // to reach a safe-drop window (without requiring an emergency drop).
  {
    sim::Ship ship;
    ship.setPositionKm({0.0, 0.0, -50'000'000.0});
    ship.setVelocityKmS({0.0, 0.0, 18'000.0});
    ship.setOrientation(math::Quatd::identity());
    ship.setAngularVelocityRadS({0.0, 0.0, 0.0});

    // Keep dampers in the inertial frame for this test.
    ship.setDampingFrameVelocityKmS({0.0, 0.0, 0.0});
    ship.setMaxLinearAccelKmS2(6.0);
    ship.setMaxAngularAccelRadS2(1.2);

    const math::Vec3d destPosKm{0.0, 0.0, 0.0};
    const math::Vec3d destVelKmS{0.0, 0.0, 0.0};
    const double dropKm = 15'000.0;

    sim::SupercruiseParams params;
    params.safeTtaSec = 7.0;
    params.safeWindowSlackSec = 2.0;
    params.maxSpeedKmS = 18'000.0;
    params.accelCapKmS2 = 6.0;
    params.angularCapRadS2 = 1.2;
    params.useBrakingDistanceLimit = true;

    bool sawSafeWindow = false;
    bool droppedNormally = false;
    bool emergency = false;

    const double dt = 0.5;
    for (int i = 0; i < 40'000; ++i) {
      const auto sc = sim::guideSupercruise(ship, destPosKm, destVelKmS, dropKm,
                                            /*navAssistEnabled=*/true,
                                            /*dropRequested=*/false,
                                            /*interdicted=*/false,
                                            params);

      if (sc.hud.safeDropReady) {
        sawSafeWindow = true;
      }

      if (sc.dropNow) {
        droppedNormally = !sc.emergencyDrop;
        emergency = sc.emergencyDrop;
        break;
      }

      ship.setMaxLinearAccelKmS2(sc.recommendedMaxLinearAccelKmS2);
      ship.setMaxAngularAccelRadS2(sc.recommendedMaxAngularAccelRadS2);
      ship.step(dt, sc.input);
    }

    if (!sawSafeWindow) {
      std::cerr << "[test_supercruise] expected to see safe-drop window at least once.\n";
      ++fails;
    }

    if (!droppedNormally) {
      const double dist = (destPosKm - ship.positionKm()).length();
      const double spd = ship.velocityKmS().length();
      std::cerr << "[test_supercruise] nav-assist did not request a normal drop. "
                << "distKm=" << dist << " speedKmS=" << spd << " emergency=" << emergency << "\n";
      ++fails;
    }
  }


  // Moving destination test: a target with significant lateral velocity.
  //
  // This is representative of chasing a traffic convoy / moving station in supercruise.
  // We expect the controller to use intercept-course lead (when enabled) and still
  // reach a normal (non-emergency) drop.
  {
    sim::Ship ship;
    ship.setPositionKm({0.0, 0.0, -35'000'000.0});
    ship.setVelocityKmS({0.0, 0.0, 18'000.0});
    ship.setOrientation(math::Quatd::identity());
    ship.setAngularVelocityRadS({0.0, 0.0, 0.0});

    ship.setDampingFrameVelocityKmS({0.0, 0.0, 0.0});
    ship.setMaxLinearAccelKmS2(6.0);
    ship.setMaxAngularAccelRadS2(1.2);

    math::Vec3d destPosKm{0.0, 0.0, 0.0};
    const math::Vec3d destVelKmS{130.0, 0.0, 0.0}; // lateral motion
    const double dropKm = 18'000.0;

    sim::SupercruiseParams params;
    params.safeTtaSec = 7.0;
    params.safeWindowSlackSec = 2.0;
    params.maxSpeedKmS = 18'000.0;
    params.accelCapKmS2 = 6.0;
    params.angularCapRadS2 = 1.2;
    params.useBrakingDistanceLimit = true;
    params.interceptEnabled = true;

    bool sawSafeWindow = false;
    bool droppedNormally = false;
    bool emergency = false;
    bool usedLead = false;

    const double dt = 0.5;
    for (int i = 0; i < 50'000; ++i) {
      // Move the destination forward in time.
      destPosKm += destVelKmS * dt;

      const auto sc = sim::guideSupercruise(ship, destPosKm, destVelKmS, dropKm,
                                            /*navAssistEnabled=*/true,
                                            /*dropRequested=*/false,
                                            /*interdicted=*/false,
                                            params);

      usedLead = usedLead || sc.hud.leadUsed;
      sawSafeWindow = sawSafeWindow || sc.hud.safeDropReady;

      if (sc.dropNow) {
        droppedNormally = !sc.emergencyDrop;
        emergency = sc.emergencyDrop;
        break;
      }

      ship.setMaxLinearAccelKmS2(sc.recommendedMaxLinearAccelKmS2);
      ship.setMaxAngularAccelRadS2(sc.recommendedMaxAngularAccelRadS2);
      ship.step(dt, sc.input);
    }

    if (!usedLead) {
      std::cerr << "[test_supercruise] expected moving-target case to use intercept lead.\n";
      ++fails;
    }

    if (!sawSafeWindow) {
      std::cerr << "[test_supercruise] moving-target case: expected to see safe-drop window at least once.\n";
      ++fails;
    }

    if (!droppedNormally) {
      const double dist = (destPosKm - ship.positionKm()).length();
      const double spd = ship.velocityKmS().length();
      std::cerr << "[test_supercruise] moving-target case did not request a normal drop. "
                << "distKm=" << dist << " speedKmS=" << spd << " emergency=" << emergency << "\n";
      ++fails;
    }
  }


  // Manual drop exactly at drop radius should be treated as a *normal* (non-emergency) drop.
  // This guards against strict '<' comparisons causing unnecessary emergency drops at the boundary.
  {
    sim::Ship ship;
    ship.setPositionKm({0.0, 0.0, -15'000.0});
    ship.setVelocityKmS({0.0, 0.0, 15'000.0 / 7.0});
    ship.setOrientation(math::Quatd::identity());
    ship.setAngularVelocityRadS({0.0, 0.0, 0.0});

    ship.setDampingFrameVelocityKmS({0.0, 0.0, 0.0});
    ship.setMaxLinearAccelKmS2(6.0);
    ship.setMaxAngularAccelRadS2(1.2);

    const math::Vec3d destPosKm{0.0, 0.0, 0.0};
    const math::Vec3d destVelKmS{0.0, 0.0, 0.0};
    const double dropKm = 15'000.0;

    sim::SupercruiseParams params;
    params.safeTtaSec = 7.0;
    params.safeWindowSlackSec = 2.0;
    params.maxSpeedKmS = 18'000.0;
    params.accelCapKmS2 = 6.0;
    params.angularCapRadS2 = 1.2;
    params.useBrakingDistanceLimit = true;

    const auto sc = sim::guideSupercruise(ship, destPosKm, destVelKmS, dropKm,
                                          /*navAssistEnabled=*/false,
                                          /*dropRequested=*/true,
                                          /*interdicted=*/false,
                                          params);

    if (!sc.dropNow || sc.emergencyDrop) {
      std::cerr << "[test_supercruise] expected manual boundary drop to be normal (non-emergency).\n";
      ++fails;
    }
  }

  // If interdicted, guideSupercruise should not force a drop even in degenerate cases.
  {
    sim::Ship ship;
    ship.setPositionKm({0.0, 0.0, 0.0});
    ship.setVelocityKmS({0.0, 0.0, 0.0});
    ship.setOrientation(math::Quatd::identity());
    ship.setAngularVelocityRadS({0.0, 0.0, 0.0});

    ship.setDampingFrameVelocityKmS({0.0, 0.0, 0.0});
    ship.setMaxLinearAccelKmS2(6.0);
    ship.setMaxAngularAccelRadS2(1.2);

    const math::Vec3d destPosKm{0.0, 0.0, 0.0};
    const math::Vec3d destVelKmS{0.0, 0.0, 0.0};
    const double dropKm = 15'000.0;

    sim::SupercruiseParams params;
    params.safeTtaSec = 7.0;
    params.safeWindowSlackSec = 2.0;
    params.maxSpeedKmS = 18'000.0;
    params.accelCapKmS2 = 6.0;
    params.angularCapRadS2 = 1.2;
    params.useBrakingDistanceLimit = true;

    const auto sc = sim::guideSupercruise(ship, destPosKm, destVelKmS, dropKm,
                                          /*navAssistEnabled=*/true,
                                          /*dropRequested=*/false,
                                          /*interdicted=*/true,
                                          params);

    if (sc.dropNow) {
      std::cerr << "[test_supercruise] expected interdiction to suppress forced drops in degenerate case.\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_supercruise] PASS\n";
  }
  return fails;
}
