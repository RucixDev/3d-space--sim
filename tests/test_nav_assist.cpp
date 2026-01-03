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

  if (fails == 0) {
    std::cout << "[test_nav_assist] PASS\n";
  }
  return fails;
}
