#include "stellar/sim/FactionProfile.h"
#include "stellar/sim/Universe.h"

#include <cmath>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

static void check01(const char* label, double v, int& fails) {
  if (!(v >= 0.0 - 1e-12 && v <= 1.0 + 1e-12) || !std::isfinite(v)) {
    std::cerr << "[test_faction_profile] expected " << label << " in [0,1], got=" << v << "\n";
    ++fails;
  }
}

int test_faction_profile() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  Universe u(1337ull);
  const core::u64 seed = u.seed();

  const auto& factions = u.factions();
  if (factions.size() < 3) {
    std::cerr << "[test_faction_profile] expected at least 3 factions (including independent)\n";
    return 1;
  }

  const core::u32 idA = factions[1].id;
  const core::u32 idB = factions[2].id;

  // --- Determinism -----------------------------------------------------------
  {
    const auto a1 = factionProfile(seed, idA);
    const auto a2 = factionProfile(seed, idA);

    if (!nearly(a1.authority, a2.authority) ||
        !nearly(a1.corruption, a2.corruption) ||
        !nearly(a1.wealth, a2.wealth) ||
        !nearly(a1.stability, a2.stability) ||
        !nearly(a1.tech, a2.tech) ||
        !nearly(a1.militarism, a2.militarism)) {
      std::cerr << "[test_faction_profile] factionProfile not deterministic\n";
      ++fails;
    }
  }

  // --- Bounds ---------------------------------------------------------------
  {
    const auto p = factionProfile(seed, idA);
    check01("authority", p.authority, fails);
    check01("corruption", p.corruption, fails);
    check01("wealth", p.wealth, fails);
    check01("stability", p.stability, fails);
    check01("tech", p.tech, fails);
    check01("militarism", p.militarism, fails);
  }

  // --- Different factions differ --------------------------------------------
  {
    const auto a = factionProfile(seed, idA);
    const auto b = factionProfile(seed, idB);

    const bool allSame =
        nearly(a.authority, b.authority) &&
        nearly(a.corruption, b.corruption) &&
        nearly(a.wealth, b.wealth) &&
        nearly(a.stability, b.stability) &&
        nearly(a.tech, b.tech) &&
        nearly(a.militarism, b.militarism);

    if (allSame) {
      std::cerr << "[test_faction_profile] expected faction profiles to differ\n";
      ++fails;
    }
  }

  // --- Relation symmetry + bounds -------------------------------------------
  {
    const double rAB = factionRelation(seed, factions[1], factions[2]);
    const double rBA = factionRelation(seed, factions[2], factions[1]);

    if (!nearly(rAB, rBA, 1e-12)) {
      std::cerr << "[test_faction_profile] factionRelation not symmetric\n";
      ++fails;
    }

    if (!(rAB >= -1.0 - 1e-12 && rAB <= 1.0 + 1e-12) || !std::isfinite(rAB)) {
      std::cerr << "[test_faction_profile] factionRelation out of bounds: " << rAB << "\n";
      ++fails;
    }

    const double rAA = factionRelation(seed, factions[1], factions[1]);
    if (!nearly(rAA, 1.0, 1e-12)) {
      std::cerr << "[test_faction_profile] expected factionRelation(self)=1, got=" << rAA << "\n";
      ++fails;
    }

    const auto k = classifyFactionRelation(rAB);
    const char* name = factionRelationKindName(k);
    if (name == nullptr || name[0] == '\0') {
      std::cerr << "[test_faction_profile] expected relation kind name\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_faction_profile] PASS\n";
  }
  return fails;
}
