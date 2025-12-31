#include "stellar/sim/Combat.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

static bool approx(double a, double b, double eps = 1e-9) {
  return std::fabs(a - b) <= eps;
}

int test_combat() {
  int fails = 0;

  // --- Damage splits shield->hull ---
  {
    double shield = 10.0;
    double hull = 25.0;
    sim::applyDamage(6.0, shield, hull);
    if (!approx(shield, 4.0) || !approx(hull, 25.0)) {
      std::cerr << "[test_combat] damage should drain shield first. got shield="
                << shield << " hull=" << hull << "\n";
      ++fails;
    }
    sim::applyDamage(9.0, shield, hull);
    if (!approx(shield, 0.0) || !approx(hull, 20.0)) {
      std::cerr << "[test_combat] overflow should spill into hull. got shield="
                << shield << " hull=" << hull << "\n";
      ++fails;
    }
  }

  // --- Ray-sphere intersection (entry distance) ---
  {
    const math::Vec3d o{0, 0, 0};
    const math::Vec3d d{0, 0, 1};
    const math::Vec3d c{0, 0, 10};
    double t = 0.0;
    const bool hit = sim::raySphereIntersectKm(o, d, c, 1.0, t);
    if (!hit || !approx(t, 9.0, 1e-9)) {
      std::cerr << "[test_combat] raySphereIntersectKm wrong. hit=" << hit << " t=" << t << " expected t=9\n";
      ++fails;
    }
  }

  // --- Nearest target selection ---
  {
    sim::SphereTarget targets[2]{};
    targets[0].kind = sim::CombatTargetKind::Ship;
    targets[0].index = 0;
    targets[0].id = 111;
    targets[0].centerKm = {0, 0, 10};
    targets[0].radiusKm = 1.0;

    targets[1].kind = sim::CombatTargetKind::Asteroid;
    targets[1].index = 1;
    targets[1].id = 222;
    targets[1].centerKm = {0, 0, 6};
    targets[1].radiusKm = 1.0;

    const auto hit = sim::raycastNearestSphereKm({0, 0, 0}, {0, 0, 1}, 50.0, targets, 2);
    if (!hit.hit || hit.id != 222) {
      std::cerr << "[test_combat] raycastNearestSphereKm should hit nearer target first. got hit="
                << hit.hit << " id=" << hit.id << "\n";
      ++fails;
    }
  }

  // --- Projectile step + hit emission ---
  {
    std::vector<sim::Projectile> ps;
    sim::Projectile p{};
    p.posKm = {0, 0, 0};
    p.prevKm = p.posKm;
    p.velKmS = {0, 0, 10};
    p.ttlSimSec = 2.0;
    p.radiusKm = 0.1;
    p.dmg = 5.0;
    p.fromPlayer = true;
    p.shooterId = 0;
    ps.push_back(p);

    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 7;
    tgt.id = 999;
    tgt.centerKm = {0, 0, 5};
    tgt.radiusKm = 0.5;

    std::vector<sim::ProjectileHit> hits;
    sim::stepProjectiles(ps, 1.0, &tgt, 1, hits);

    if (hits.size() != 1 || hits[0].targetId != 999 || !approx(hits[0].dmg, 5.0)) {
      std::cerr << "[test_combat] stepProjectiles should emit one hit. hits=" << hits.size()
                << " id=" << (hits.empty() ? 0 : hits[0].targetId) << "\n";
      ++fails;
    }

    if (ps.empty() || ps[0].ttlSimSec > 0.0) {
      std::cerr << "[test_combat] projectile should be consumed on hit. ttl="
                << (ps.empty() ? -1.0 : ps[0].ttlSimSec) << "\n";
      ++fails;
    }
  }

  // --- Missile step + detonation + splash hit ---
  {
    std::vector<sim::Missile> ms;
    sim::Missile m{};
    m.posKm = {0, 0, 0};
    m.prevKm = m.posKm;
    m.velKmS = {0, 0, 10};
    m.ttlSimSec = 2.0;
    m.radiusKm = 0.1;
    m.blastRadiusKm = 2.0;
    m.dmg = 5.0;
    m.fromPlayer = true;
    m.shooterId = 0;
    m.hasTarget = true;
    m.targetKind = sim::CombatTargetKind::Ship;
    m.targetId = 999;
    ms.push_back(m);

    sim::SphereTarget tgt{};
    tgt.kind = sim::CombatTargetKind::Ship;
    tgt.index = 7;
    tgt.id = 999;
    tgt.centerKm = {0, 0, 5};
    tgt.velKmS = {0, 0, 0};
    tgt.radiusKm = 0.5;

    std::vector<sim::MissileDetonation> dets;
    std::vector<sim::MissileHit> hits;
    sim::stepMissiles(ms, 1.0, &tgt, 1, dets, hits);

    if (dets.size() != 1) {
      std::cerr << "[test_combat] stepMissiles should emit one detonation. got=" << dets.size() << "\n";
      ++fails;
    }
    if (hits.size() != 1 || hits[0].targetId != 999 || !approx(hits[0].dmg, 5.0)) {
      std::cerr << "[test_combat] stepMissiles should splash-hit the target. hits=" << hits.size()
                << " id=" << (hits.empty() ? 0 : hits[0].targetId)
                << " dmg=" << (hits.empty() ? 0.0 : hits[0].dmg) << "\n";
      ++fails;
    }
    if (ms.empty() || ms[0].ttlSimSec > 0.0) {
      std::cerr << "[test_combat] missile should be consumed on detonation. ttl="
                << (ms.empty() ? -1.0 : ms[0].ttlSimSec) << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_combat] PASS\n";
  }
  return fails;
}
