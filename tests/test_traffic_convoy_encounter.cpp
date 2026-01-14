#include "test_harness.h"

#include "stellar/sim/TrafficConvoyEncounter.h"

#include <cmath>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_traffic_convoy_encounter() {
  int failures = 0;

  using namespace stellar;
  using namespace stellar::sim;

  // --- Determinism: same inputs => identical plan --------------------------
  {
    SystemSecurityProfile sec{};
    sec.security01 = 0.78;
    sec.piracy01 = 0.25;
    sec.traffic01 = 0.60;

    TrafficConvoyEncounterParams p{};
    p.npcShieldRegenPerSec = 0.035;

    const auto a = planTrafficConvoyEncounter(1337, 42, 12000.0, sec, /*lawful=*/true, /*alivePiratesNow=*/0, p);
    const auto b = planTrafficConvoyEncounter(1337, 42, 12000.0, sec, /*lawful=*/true, /*alivePiratesNow=*/0, p);

    CHECK(a.ok);
    CHECK(b.ok);
    CHECK(a.seed == b.seed);
    CHECK(nearly(a.value01, b.value01));
    CHECK(nearly(a.risk01, b.risk01));
    CHECK(a.escortCount == b.escortCount);
    CHECK(a.pirateCount == b.pirateCount);
    CHECK(a.ambush == b.ambush);
    CHECK(a.npcs.size() == b.npcs.size());

    for (std::size_t i = 0; i < a.npcs.size(); ++i) {
      const auto& x = a.npcs[i];
      const auto& y = b.npcs[i];

      CHECK(x.role == y.role);
      CHECK(x.leader == y.leader);
      CHECK(x.hullClass == y.hullClass);
      CHECK(x.weapon == y.weapon);
      CHECK(x.thrusterMk == y.thrusterMk);
      CHECK(x.shieldMk == y.shieldMk);
      CHECK(x.distributorMk == y.distributorMk);

      CHECK(nearly(x.posOffsetKm.x, y.posOffsetKm.x));
      CHECK(nearly(x.posOffsetKm.y, y.posOffsetKm.y));
      CHECK(nearly(x.posOffsetKm.z, y.posOffsetKm.z));
      CHECK(nearly(x.velOffsetKmS.x, y.velOffsetKmS.x));
      CHECK(nearly(x.velOffsetKmS.y, y.velOffsetKmS.y));
      CHECK(nearly(x.velOffsetKmS.z, y.velOffsetKmS.z));
    }
  }

  // --- Lawful vs anarchy: escorts suppressed when lawfulSystem=false -------
  {
    SystemSecurityProfile sec{};
    sec.security01 = 0.95;
    sec.piracy01 = 0.10;
    sec.traffic01 = 0.60;

    TrafficConvoyEncounterParams p{};
    const auto lawful = planTrafficConvoyEncounter(9, 777, 20000.0, sec, /*lawful=*/true, /*alivePiratesNow=*/0, p);
    const auto anarchy = planTrafficConvoyEncounter(9, 777, 20000.0, sec, /*lawful=*/false, /*alivePiratesNow=*/0, p);

    CHECK(lawful.ok);
    CHECK(anarchy.ok);
    CHECK(lawful.escortCount >= 1);
    CHECK(anarchy.escortCount == 0);
  }

  // --- Spawn bands: escorts/pirates spawn relative to trader within ranges -
  {
    SystemSecurityProfile sec{};
    sec.security01 = 0.65;
    sec.piracy01 = 0.90;
    sec.traffic01 = 0.70;

    TrafficConvoyEncounterParams p{};
    p.pirateAliveCap = 100; // ensure ambush can happen

    const auto plan = planTrafficConvoyEncounter(555, 999, 50000.0, sec, /*lawful=*/true, /*alivePiratesNow=*/0, p);
    CHECK(plan.ok);
    CHECK(!plan.npcs.empty());
    CHECK(plan.npcs[0].role == TrafficEncounterNpcRole::Trader);

    const math::Vec3d traderOff = plan.npcs[0].posOffsetKm;
    const double traderDist = traderOff.length();
    CHECK(traderDist >= p.convoyOffsetMinKm - 1e-6);
    CHECK(traderDist <= p.convoyOffsetMaxKm + 1e-6);

    // Verify escort offsets are near the expected band around the trader.
    for (std::size_t i = 1; i < plan.npcs.size(); ++i) {
      const auto& n = plan.npcs[i];
      if (n.role == TrafficEncounterNpcRole::Police) {
        const double d = (n.posOffsetKm - traderOff).length();
        CHECK(d >= p.escortOffsetMinKm - 1e-6);
        CHECK(d <= p.escortOffsetMaxKm + 1e-6);
        CHECK(n.velOffsetKmS.length() <= p.escortVelJitterMaxKmS + 1e-6);
      }
      if (n.role == TrafficEncounterNpcRole::Pirate) {
        const double d = (n.posOffsetKm - traderOff).length();
        CHECK(d >= p.pirateOffsetMinKm - 1e-6);
        CHECK(d <= p.pirateOffsetMaxKm + 1e-6);
        CHECK(n.velOffsetKmS.length() <= p.pirateVelJitterMaxKmS + 1e-6);
      }
    }
  }

  // --- Pirate cap suppresses additional raiders ----------------------------
  {
    SystemSecurityProfile sec{};
    sec.security01 = 0.10;
    sec.piracy01 = 1.00;
    sec.traffic01 = 0.60;

    TrafficConvoyEncounterParams p{};
    p.pirateAliveCap = 5;

    const auto plan = planTrafficConvoyEncounter(777, 1234, 80000.0, sec, /*lawful=*/true, /*alivePiratesNow=*/5, p);
    CHECK(plan.ok);
    CHECK(plan.pirateCount == 0);
    CHECK(plan.ambush == false);
  }

  return failures;
}
