#include "test_harness.h"

#include "stellar/sim/ExplorationData.h"

#include <cmath>

using namespace stellar;

int test_exploration_data() {
  int failures = 0;

  // ---------------------------------------------------------------------------
  // Scan key stability + basic collision resistance
  // ---------------------------------------------------------------------------
  {
    const sim::SystemId sysA = 111;
    const sim::SystemId sysB = 222;

    const core::u64 starA = sim::scanKeyStar(sysA);
    const core::u64 starA2 = sim::scanKeyStar(sysA);
    const core::u64 starB = sim::scanKeyStar(sysB);

    CHECK(starA == starA2);
    CHECK(starA != starB);

    const core::u64 p0 = sim::scanKeyPlanet(sysA, 0);
    const core::u64 p1 = sim::scanKeyPlanet(sysA, 1);
    const core::u64 p0b = sim::scanKeyPlanet(sysB, 0);

    CHECK(p0 != p1);
    CHECK(p0 != p0b);
    CHECK(p0 != starA);

    const core::u64 st = sim::scanKeyStation((sim::StationId)9001);
    const core::u64 sg = sim::scanKeySignal((core::u64)5555);
    const core::u64 as = sim::scanKeyAsteroid((core::u64)7777);
    const core::u64 comp = sim::scanKeySystemComplete(sysA);

    CHECK(st != starA);
    CHECK(st != p0);
    CHECK(sg != st);
    CHECK(as != sg);
    CHECK(comp != starA);
    CHECK(comp != p0);
  }

  // ---------------------------------------------------------------------------
  // Value ladders: finite + positive + rough ordering
  // ---------------------------------------------------------------------------
  {
    sim::Star m{};
    m.cls = sim::StarClass::M;
    m.massSol = 0.3;
    m.luminositySol = 0.05;
    m.temperatureK = 3200.0;

    sim::Star o{};
    o.cls = sim::StarClass::O;
    o.massSol = 35.0;
    o.luminositySol = 120.0;
    o.temperatureK = 35000.0;

    const double vM = sim::scanValueStar(m);
    const double vO = sim::scanValueStar(o);

    CHECK(std::isfinite(vM) && vM > 0.0);
    CHECK(std::isfinite(vO) && vO > 0.0);
    CHECK(vO > vM);

    sim::Planet rocky{};
    rocky.type = sim::PlanetType::Rocky;
    rocky.radiusEarth = 1.0;
    rocky.massEarth = 1.0;

    sim::Planet ocean{};
    ocean.type = sim::PlanetType::Ocean;
    ocean.radiusEarth = 1.0;
    ocean.massEarth = 1.0;

    const double vR = sim::scanValuePlanet(rocky);
    const double vOc = sim::scanValuePlanet(ocean);

    CHECK(std::isfinite(vR) && vR > 0.0);
    CHECK(std::isfinite(vOc) && vOc > 0.0);
    CHECK(vOc > vR);

    sim::Station out{};
    out.type = econ::StationType::Outpost;
    out.radiusKm = 4000.0;

    sim::Station res{};
    res.type = econ::StationType::Research;
    res.radiusKm = 4000.0;

    const double vOut = sim::scanValueStation(out);
    const double vRes = sim::scanValueStation(res);

    CHECK(std::isfinite(vOut) && vOut > 0.0);
    CHECK(std::isfinite(vRes) && vRes > 0.0);
    CHECK(vRes > vOut);

    const double vSig = sim::scanValueSignal(sim::SignalKind::Derelict);
    CHECK(std::isfinite(vSig) && vSig > 0.0);

    const double vAst = sim::scanValueAsteroidProspect();
    CHECK(std::isfinite(vAst) && vAst > 0.0);

    const double vBonus = sim::scanValueSystemSurveyBonus(6, 2);
    CHECK(std::isfinite(vBonus) && vBonus > 0.0);
  }

  // ---------------------------------------------------------------------------
  // Broker multipliers: distance premium, faction bonus, station demand
  // ---------------------------------------------------------------------------
  {
    sim::SystemStub sale{};
    sale.id = 1;
    sale.posLy = {0, 0, 0};
    sale.factionId = 5;

    sim::SystemStub scan{};
    scan.id = 2;
    scan.posLy = {100, 0, 0};
    scan.factionId = 5;

    sim::Station st{};
    st.factionId = 5;
    st.type = econ::StationType::Research;

    sim::ExplorationDataBrokerParams p{};
    p.distanceScaleLy = 100.0;
    p.maxDistancePremium = 0.30;
    p.sameFactionBonus = 0.10;
    p.otherFactionPenalty = 0.05;
    p.enableStationDemand = true;
    p.demandStrength = 1.0;

    const double mStar = sim::explorationDataBrokerMultiplier(p, sale, st, scan,
                                                             sim::LogbookEntryKind::StarScan, 0);
    CHECK(std::isfinite(mStar));
    CHECK(mStar > 1.0);

    // Different faction should reduce vs same faction.
    scan.factionId = 9;
    const double mOther = sim::explorationDataBrokerMultiplier(p, sale, st, scan,
                                                              sim::LogbookEntryKind::StarScan, 0);
    CHECK(std::isfinite(mOther));
    CHECK(mOther < mStar);
  }

  return failures;
}
