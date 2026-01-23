#include "stellar/sim/StationSecurity.h"

#include "test_harness.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static sim::Station makeTestStation() {
  sim::Station st{};
  st.id = 1;
  st.factionId = 42;
  st.name = "Test Station";
  st.radiusKm = 6000.0;
  st.slotWidthKm = 500.0;
  st.slotHeightKm = 250.0;
  st.slotDepthKm = 2000.0;
  st.maxApproachSpeedKmS = 0.20;
  return st;
}

static sim::LawProfile makeTestLaw() {
  sim::LawProfile law{};
  law.scanStrictness = 1.25;
  law.corruption = 0.20;
  law.fineBaseCr = 300.0;
  law.fineRate = 0.75;
  law.bribeBaseCr = 260.0;
  law.bribeRate = 0.55;
  law.repBase = 3.2;
  law.repDiv = 3500.0;
  law.repMin = 2.0;
  law.repMax = 10.0;
  law.evadeRepMult = 1.35;
  return law;
}

int test_station_security() {
  int failures = 0;

  const sim::Station st = makeTestStation();
  const sim::LawProfile law = makeTestLaw();

  sim::StationSecurityParams params{};
  params.weaponEventCooldownSec = 0.50;
  params.speedingTriggerSec = 0.60;
  params.trespassTriggerSec = 0.40;
  params.genericCooldownSec = 0.0; // tests want fast re-trigger

  const math::Vec3d stPos{0, 0, 0};
  const math::Vec3d stVel{0, 0, 0};
  const math::Quatd stQ = math::Quatd::identity();

  // --- Weapon discharge escalation ---
  {
    sim::StationSecurityState sec{};
    const double dist = sim::stationNoFireZoneRadiusKm(st, params) * 0.5;
    const math::Vec3d shipPos{0, 0, dist};
    const math::Vec3d shipVel{0, 0, 0};

    double tDays = 100.0;

    auto e1 = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                        shipPos, shipVel,
                                        /*hasClearance=*/false,
                                        /*shipDocked=*/false,
                                        /*shipFiredWeaponThisTick=*/true,
                                        tDays, /*dtSimSec=*/1.0);
    CHECK(e1.has_value());
    CHECK(e1->kind == sim::StationOffenseKind::WeaponDischarge);
    CHECK(e1->action == sim::StationOffenseAction::Warning);
    CHECK(e1->strikes == 1);

    // Within cooldown: should not emit again.
    tDays += 0.25 / 86400.0;
    auto eSpam = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                           shipPos, shipVel,
                                           false, false, true, tDays, 0.25);
    CHECK(!eSpam.has_value());

    // Beyond cooldown: second strike should produce a fine.
    tDays += 1.00 / 86400.0;
    auto e2 = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                        shipPos, shipVel,
                                        false, false, true, tDays, 1.0);
    CHECK(e2.has_value());
    CHECK(e2->action == sim::StationOffenseAction::Fine);
    CHECK(e2->fineCr > 1.0);
    CHECK(e2->repPenalty < -1e-6);
    CHECK(e2->strikes == 2);

    // Third strike: bounty.
    tDays += 1.00 / 86400.0;
    auto e3 = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                        shipPos, shipVel,
                                        false, false, true, tDays, 1.0);
    CHECK(e3.has_value());
    CHECK(e3->action == sim::StationOffenseAction::Bounty);
    CHECK(e3->bountyCr > 1.0);
    CHECK(e3->repPenalty < -1e-6);
    CHECK(e3->strikes == 3);
  }

  // --- Speeding accumulation ---
  {
    sim::StationSecurityState sec{};
    const double dist = sim::stationSpeedZoneRadiusKm(st, params) * 0.55;
    const math::Vec3d shipPos{0, 0, dist};

    const double limit = sim::stationSpeedLimitKmS(st, params, dist);
    const double relSpeed = limit * (1.0 + params.speedToleranceFrac + 0.60);
    const math::Vec3d shipVel{relSpeed, 0, 0};

    double tDays = 200.0;

    // Accumulate just under trigger -> no event.
    auto e0 = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                        shipPos, shipVel,
                                        false, false, false, tDays, 0.30);
    CHECK(!e0.has_value());

    // Cross trigger: warning.
    tDays += 0.30 / 86400.0;
    auto e1 = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                        shipPos, shipVel,
                                        false, false, false, tDays, 0.30);
    CHECK(e1.has_value());
    CHECK(e1->kind == sim::StationOffenseKind::Speeding);
    CHECK(e1->action == sim::StationOffenseAction::Warning);
    CHECK(e1->strikes == 1);
    CHECK(e1->measuredSpeedKmS > e1->speedLimitKmS);

    // Next trigger: fine.
    tDays += 0.60 / 86400.0;
    auto e2 = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                        shipPos, shipVel,
                                        false, false, false, tDays, 0.60);
    CHECK(e2.has_value());
    CHECK(e2->action == sim::StationOffenseAction::Fine);
    CHECK(e2->fineCr > 1.0);
    CHECK(e2->repPenalty < -1e-6);
    CHECK(e2->strikes == 2);
  }

  // --- Trespass (slot tunnel entry without clearance) ---
  {
    sim::StationSecurityState sec{};

    const double wz = st.radiusKm * 1.10;
    const math::Vec3d relLocal{0, 0, wz - 100.0};
    CHECK(sim::insideStationSlotTunnel(st, relLocal));

    const math::Vec3d shipPos = relLocal; // station at origin, identity orient
    const math::Vec3d shipVel{0, 0, 0};

    double tDays = 300.0;

    auto e0 = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                        shipPos, shipVel,
                                        /*hasClearance=*/false,
                                        /*shipDocked=*/false,
                                        /*shipFiredWeaponThisTick=*/false,
                                        tDays, 0.20);
    CHECK(!e0.has_value());

    // Cross trigger -> warning.
    tDays += 0.20 / 86400.0;
    auto e1 = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                        shipPos, shipVel,
                                        false, false, false, tDays, 0.20);
    CHECK(e1.has_value());
    CHECK(e1->kind == sim::StationOffenseKind::Trespass);
    CHECK(e1->action == sim::StationOffenseAction::Warning);
    CHECK(e1->strikes == 1);

    // Next trigger -> fine.
    tDays += 0.40 / 86400.0;
    auto e2 = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                        shipPos, shipVel,
                                        false, false, false, tDays, 0.40);
    CHECK(e2.has_value());
    CHECK(e2->action == sim::StationOffenseAction::Fine);
    CHECK(e2->fineCr > 1.0);
    CHECK(e2->repPenalty < -1e-6);
    CHECK(e2->strikes == 2);

    // Clearance suppresses trespass.
    tDays += 0.40 / 86400.0;
    auto eSuppressed = sim::updateStationSecurity(sec, params, st, law, stPos, stVel, stQ,
                                                 shipPos, shipVel,
                                                 /*hasClearance=*/true,
                                                 /*shipDocked=*/false,
                                                 /*shipFiredWeaponThisTick=*/false,
                                                 tDays, 0.80);
    CHECK(!eSuppressed.has_value());
  }

  if (failures == 0) {
    std::cout << "[test_station_security] PASS\n";
  }
  return failures;
}
