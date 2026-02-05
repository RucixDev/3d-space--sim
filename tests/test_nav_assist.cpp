#include "stellar/sim/NavAssistComputer.h"

#include <cmath>
#include <iostream>

using namespace stellar;

int test_nav_assist() {
  int fails = 0;

  // --- Approach: converge to a standoff distance while facing the target.
  {
    sim::Ship s;
    s.setPositionKm({0, 0, 0});
    s.setVelocityKmS({0, 0, 0});
    s.setOrientation({1, 0, 0, 0});
    s.setAngularVelocityRadS({0, 0, 0});

    sim::NavAssistComputer nav;
    sim::NavAssistParams p = nav.params();
    p.approachMaxSpeedKmS = 0.60;
    p.approachSlowDownRangeKm = 15.0; // start slowing down within ~15km
    p.approachVelGain = 2.2;
    p.approachAllowBoost = false;
    p.arriveDistEpsKm = 0.6;
    p.arriveRelSpeedEpsKmS = 0.03;
    nav.setParams(p);

    const math::Vec3d targetPos{0, 0, 20.0};
    const math::Vec3d targetVel{0, 0, 0};
    const double desiredDist = 2.0;

    nav.engageApproach(desiredDist);

    const double dt = 0.1;
    for (int i = 0; i < 1200; ++i) {
      const auto out = nav.update(s, targetPos, targetVel, dt);
      s.step(dt, out.input);
    }

    const double dist = (targetPos - s.positionKm()).length();
    const double relSpeed = (s.velocityKmS() - targetVel).length();

    if (!(dist > 1.0 && dist < 3.5)) {
      std::cerr << "[test_nav_assist] Approach: distance did not converge. dist=" << dist << "\n";
      ++fails;
    }
    if (relSpeed > 0.08) {
      std::cerr << "[test_nav_assist] Approach: relative speed too high. relSpeed=" << relSpeed << "\n";
      ++fails;
    }
  }

  // --- Match velocity: hold separation while matching the target's motion.
  {
    sim::Ship s;
    s.setPositionKm({0, 0, 0});
    s.setVelocityKmS({0, 0, 0});
    s.setOrientation({1, 0, 0, 0});
    s.setAngularVelocityRadS({0, 0, 0});

    sim::NavAssistComputer nav;
    sim::NavAssistParams p = nav.params();
    p.matchMaxSpeedKmS = 0.60;
    p.matchSlowDownRangeKm = 18.0;
    p.matchVelGain = 2.8;
    p.matchAllowBoost = false;
    p.arriveDistEpsKm = 2.0;
    p.arriveRelSpeedEpsKmS = 0.03;
    nav.setParams(p);

    math::Vec3d targetPos{0, 0, 20.0};
    const math::Vec3d targetVel{0.20, 0.0, 0.0};

    nav.engageMatchVelocity(s, targetPos);

    const double desiredDist = (targetPos - s.positionKm()).length();
    const double dt = 0.1;
    for (int i = 0; i < 1200; ++i) {
      targetPos = targetPos + targetVel * dt;

      const auto out = nav.update(s, targetPos, targetVel, dt);
      s.step(dt, out.input);
    }

    const double dist = (targetPos - s.positionKm()).length();
    const double relSpeed = (s.velocityKmS() - targetVel).length();

    if (std::abs(dist - desiredDist) > 5.0) {
      std::cerr << "[test_nav_assist] MatchVelocity: distance drifted too far. dist=" << dist
                << " desired=" << desiredDist << "\n";
      ++fails;
    }
    if (relSpeed > 0.08) {
      std::cerr << "[test_nav_assist] MatchVelocity: relative speed too high. relSpeed=" << relSpeed << "\n";
      ++fails;
    }
  }



  // --- Orbit: hold a standoff distance while maintaining tangential motion.
  {
    sim::Ship s;
    s.setPositionKm({10.0, 0, 0});
    s.setVelocityKmS({0, 0, 0});
    s.setOrientation({1, 0, 0, 0});
    s.setAngularVelocityRadS({0, 0, 0});

    sim::NavAssistComputer nav;
    sim::NavAssistParams p = nav.params();
    p.orbitMaxSpeedKmS = 0.65;
    p.orbitSlowDownRangeKm = 12.0;
    p.orbitVelGain = 2.5;
    p.orbitAllowBoost = false;
    p.orbitTangentialSpeedKmS = 0.25;
    p.orbitLeadTimeSec = 3.0;
    p.orbitLeadMaxFrac = 0.60;
    p.orbitFaceTarget = true;

    // Allow the controller to declare a stable orbit if needed.
    p.arriveDistEpsKm = 2.0;
    p.arriveRelSpeedEpsKmS = 0.06;
    nav.setParams(p);

    const math::Vec3d targetPos{0, 0, 0};
    const math::Vec3d targetVel{0, 0, 0};
    const double desiredDist = 10.0;

    nav.engageOrbit(s, targetPos, targetVel, desiredDist);

    const double dt = 0.1;
    for (int i = 0; i < 2500; ++i) {
      const auto out = nav.update(s, targetPos, targetVel, dt);
      s.step(dt, out.input);
    }

    const math::Vec3d relPos = s.positionKm() - targetPos;
    const math::Vec3d relVel = s.velocityKmS() - targetVel;

    const double dist = relPos.length();
    const double speed = relVel.length();
    const math::Vec3d radial = relPos.normalized();
    const double radialSpeed = math::dot(relVel, radial);
    const math::Vec3d h = math::cross(relPos, relVel);

    if (std::abs(dist - desiredDist) > 3.0) {
      std::cerr << "[test_nav_assist] Orbit: standoff distance not maintained. dist=" << dist
                << " desired=" << desiredDist << "\n";
      ++fails;
    }

    if (speed < 0.05) {
      std::cerr << "[test_nav_assist] Orbit: did not build tangential velocity. speed=" << speed << "\n";
      ++fails;
    }

    // Tangential motion should dominate radial drift.
    if (std::abs(radialSpeed) > speed * 0.75) {
      std::cerr << "[test_nav_assist] Orbit: motion too radial. radialSpeed=" << radialSpeed
                << " speed=" << speed << "\n";
      ++fails;
    }

    // With identity orientation, orbit normal is expected to be +Y.
    if (h.y <= 0.0) {
      std::cerr << "[test_nav_assist] Orbit: wrong orbit direction / plane. h.y=" << h.y << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_nav_assist] PASS\n";
  }
  return fails;
}
