#include "stellar/sim/Extortion.h"

#include <cmath>
#include <iostream>

int test_extortion() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::sim;

  const core::u64 seed = 42ull;
  const core::u64 targetId = 1234567ull;

  // No cargo => no offer.
  {
    const auto p = planExtortionDemand(seed, targetId, /*targetCargoValueCr=*/50.0,
                                       /*attackerStrength=*/100.0,
                                       /*defenderStrength=*/100.0,
                                       /*systemSecurity01=*/0.5,
                                       /*policeNearby=*/false);
    if (p.offer) {
      std::cerr << "[test_extortion] expected offer=false for low cargo value\n";
      ++fails;
    }
  }

  // Baseline sanity.
  {
    const auto p = planExtortionDemand(seed, targetId, /*targetCargoValueCr=*/10000.0,
                                       /*attackerStrength=*/300.0,
                                       /*defenderStrength=*/300.0,
                                       /*systemSecurity01=*/0.2,
                                       /*policeNearby=*/false);
    if (!p.offer) {
      std::cerr << "[test_extortion] expected offer=true for meaningful cargo value\n";
      ++fails;
    }
    if (!(p.demandedValueCr >= 80.0 && p.demandedValueCr <= 12000.0)) {
      std::cerr << "[test_extortion] demandedValueCr out of bounds: " << p.demandedValueCr << "\n";
      ++fails;
    }
    if (!(p.complyChance01 >= 0.0 && p.complyChance01 <= 1.0)) {
      std::cerr << "[test_extortion] complyChance01 out of bounds: " << p.complyChance01 << "\n";
      ++fails;
    }
    if (!(p.suggestedCooldownSec >= 30.0 && p.suggestedCooldownSec <= 140.0)) {
      std::cerr << "[test_extortion] suggestedCooldownSec out of bounds: " << p.suggestedCooldownSec << "\n";
      ++fails;
    }
    if (!(p.suggestedFleeSec >= 60.0 && p.suggestedFleeSec <= 260.0)) {
      std::cerr << "[test_extortion] suggestedFleeSec out of bounds: " << p.suggestedFleeSec << "\n";
      ++fails;
    }
  }

  // Stronger attacker => higher compliance and (usually) higher demanded value.
  {
    const auto weak = planExtortionDemand(seed, targetId, 10000.0, 220.0, 800.0, 0.2, false);
    const auto strong = planExtortionDemand(seed, targetId, 10000.0, 800.0, 220.0, 0.2, false);

    if (!(strong.complyChance01 > weak.complyChance01 + 1e-6)) {
      std::cerr << "[test_extortion] expected strong attacker to have higher complyChance01\n";
      std::cerr << "  weak=" << weak.complyChance01 << " strong=" << strong.complyChance01 << "\n";
      ++fails;
    }

    if (!(strong.demandedValueCr > weak.demandedValueCr - 1e-6)) {
      std::cerr << "[test_extortion] expected strong attacker to demand >= weak attacker\n";
      std::cerr << "  weak=" << weak.demandedValueCr << " strong=" << strong.demandedValueCr << "\n";
      ++fails;
    }
  }

  // Higher security => lower compliance and lower demand.
  {
    const auto lowSec = planExtortionDemand(seed, targetId, 10000.0, 500.0, 250.0, 0.0, false);
    const auto highSec = planExtortionDemand(seed, targetId, 10000.0, 500.0, 250.0, 1.0, false);

    if (!(lowSec.complyChance01 > highSec.complyChance01 + 1e-6)) {
      std::cerr << "[test_extortion] expected low security to increase compliance\n";
      std::cerr << "  lowSec=" << lowSec.complyChance01 << " highSec=" << highSec.complyChance01 << "\n";
      ++fails;
    }
    if (!(lowSec.demandedValueCr > highSec.demandedValueCr - 1e-6)) {
      std::cerr << "[test_extortion] expected low security to increase demanded value\n";
      std::cerr << "  lowSec=" << lowSec.demandedValueCr << " highSec=" << highSec.demandedValueCr << "\n";
      ++fails;
    }
  }

  // Nearby authorities => lower compliance.
  {
    const auto clear = planExtortionDemand(seed, targetId, 10000.0, 600.0, 260.0, 0.2, false);
    const auto police = planExtortionDemand(seed, targetId, 10000.0, 600.0, 260.0, 0.2, true);
    if (!(clear.complyChance01 > police.complyChance01 + 1e-6)) {
      std::cerr << "[test_extortion] expected policeNearby to reduce compliance\n";
      std::cerr << "  clear=" << clear.complyChance01 << " police=" << police.complyChance01 << "\n";
      ++fails;
    }
  }

  return fails;
}
