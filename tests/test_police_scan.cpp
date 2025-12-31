#include "stellar/sim/PoliceScan.h"

#include "stellar/sim/Contraband.h"

#include <algorithm>
#include <cmath>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_police_scan() {
  int fails = 0;

  using namespace stellar;
  using namespace stellar::econ;
  using namespace stellar::sim;

  const core::u64 seed = 1337ull;
  const core::u32 factionId = 7;

  // --- scanIllegalCargo ------------------------------------------------------
  {
    std::array<double, kCommodityCount> cargo{};
    cargo.fill(0.0);

    const core::u32 mask = illegalCommodityMask(seed, factionId);
    if (mask == 0u) {
      std::cerr << "[test_police_scan] expected some illegal commodities for deterministic test\n";
      ++fails;
    } else {
      // Pick up to two illegal commodities deterministically.
      std::size_t i1 = kCommodityCount;
      std::size_t i2 = kCommodityCount;
      for (std::size_t i = 0; i < kCommodityCount; ++i) {
        if ((mask & ((core::u32)1u << (core::u32)i)) != 0u) {
          if (i1 == kCommodityCount) i1 = i;
          else { i2 = i; break; }
        }
      }

      if (i1 == kCommodityCount) {
        std::cerr << "[test_police_scan] failed to pick an illegal commodity\n";
        ++fails;
      } else {
        cargo[i1] = 10.0;
        if (i2 != kCommodityCount) cargo[i2] = 5.0;

        // Also set a legal commodity for noise.
        for (std::size_t i = 0; i < kCommodityCount; ++i) {
          if ((mask & ((core::u32)1u << (core::u32)i)) == 0u) {
            cargo[i] = 3.0;
            break;
          }
        }

        std::array<double, kCommodityCount> mid{};
        for (std::size_t i = 0; i < kCommodityCount; ++i) {
          mid[i] = commodityDef((CommodityId)i).basePrice;
        }
        mid[i1] = 100.0;
        if (i2 != kCommodityCount) mid[i2] = 250.0;

        const auto res = scanIllegalCargo(seed, factionId, cargo, &mid);

        const double expected = 10.0 * 100.0 + ((i2 != kCommodityCount) ? 5.0 * 250.0 : 0.0);
        if (!nearly(res.illegalValueCr, expected, 1e-6)) {
          std::cerr << "[test_police_scan] illegalValueCr mismatch: got=" << res.illegalValueCr
                    << " expected=" << expected << "\n";
          ++fails;
        }

        if (!nearly(res.scannedIllegalUnits[i1], 10.0)) {
          std::cerr << "[test_police_scan] scannedIllegalUnits[i1] mismatch\n";
          ++fails;
        }
        if (i2 != kCommodityCount && !nearly(res.scannedIllegalUnits[i2], 5.0)) {
          std::cerr << "[test_police_scan] scannedIllegalUnits[i2] mismatch\n";
          ++fails;
        }

        // Legal commodities should not be flagged.
        for (std::size_t i = 0; i < kCommodityCount; ++i) {
          if ((mask & ((core::u32)1u << (core::u32)i)) == 0u) {
            if (res.scannedIllegalUnits[i] > 1e-9) {
              std::cerr << "[test_police_scan] legal commodity was marked illegal\n";
              ++fails;
              break;
            }
          }
        }
      }
    }

    // factionId=0 should always scan clean.
    {
      std::array<double, kCommodityCount> cargo2{};
      cargo2.fill(5.0);
      const auto res2 = scanIllegalCargo(seed, 0, cargo2, nullptr);
      if (res2.illegalValueCr > 1e-9) {
        std::cerr << "[test_police_scan] expected illegalValueCr=0 for factionId=0\n";
        ++fails;
      }
    }
  }

  // --- Scan rates / durations ------------------------------------------------
  {
    LawProfile law{};
    law.scanStrictness = 1.2;

    const double r0 = cargoScanStartRatePerSec(false, law, 0);
    const double r1 = cargoScanStartRatePerSec(false, law, 1);
    const double r2 = cargoScanStartRatePerSec(false, law, 2);
    const double r3 = cargoScanStartRatePerSec(false, law, 3);

    if (!(r0 > 0.0 && r1 > 0.0 && r2 > 0.0 && r3 > 0.0)) {
      std::cerr << "[test_police_scan] expected positive scan rates\n";
      ++fails;
    }
    if (!(r0 > r1 && r1 > r2 && r2 > r3)) {
      std::cerr << "[test_police_scan] expected scan rate to decrease with better holds\n";
      ++fails;
    }

    const double d0 = cargoScanDurationSecStation(true, 0);
    const double d3 = cargoScanDurationSecStation(true, 3);
    if (!(d3 > d0)) {
      std::cerr << "[test_police_scan] expected scan duration to increase with better holds\n";
      ++fails;
    }
  }

  // --- Bribe offer chance ----------------------------------------------------
  {
    LawProfile law{};
    law.corruption = 0.8;
    law.scanStrictness = 1.1;

    const double rep = 10.0;
    const double heat = 20.0;
    const double illegalValue = 5000.0;

    const double chance = bribeOfferChance(law, rep, heat, illegalValue);
    if (chance < 0.0 || chance > 1.0) {
      std::cerr << "[test_police_scan] bribeOfferChance out of bounds\n";
      ++fails;
    }

    // Compare against the reference formula in the prototype.
    double expected = 0.33;
    expected *= std::clamp(law.corruption, 0.0, 1.0);
    expected *= std::clamp(1.10 - (heat / 110.0), 0.25, 1.10);
    expected *= std::clamp(1.05 - (illegalValue / 20000.0), 0.35, 1.05);
    expected = std::clamp(expected, 0.0, 1.0);

    if (!nearly(chance, expected, 1e-12)) {
      std::cerr << "[test_police_scan] bribeOfferChance mismatch\n";
      ++fails;
    }

    // Bad rep should strongly reduce the chance.
    const double bad = bribeOfferChance(law, -50.0, heat, illegalValue);
    if (!(bad < chance)) {
      std::cerr << "[test_police_scan] expected bad rep to reduce bribe chance\n";
      ++fails;
    }
  }

  // --- Enforcement math ------------------------------------------------------
  {
    LawProfile law{};
    law.scanStrictness = 1.4;
    law.fineBaseCr = 100.0;
    law.fineRate = 0.50;
    law.repBase = 3.0;
    law.repDiv = 4000.0;
    law.repMin = 2.0;
    law.repMax = 9.0;
    law.evadeRepMult = 1.35;

    std::array<double, kCommodityCount> cargo{};
    cargo.fill(0.0);

    // Mark commodity 0 as "scanned illegal" for the purposes of this pure function.
    std::array<double, kCommodityCount> scanned{};
    scanned.fill(0.0);
    cargo[0] = 12.0;
    scanned[0] = 10.0; // scanner saw 10 units

    const double credits = 120.0;
    const double illegalValue = 1000.0;

    const auto res = enforceContraband(law, credits, cargo, scanned, illegalValue);

    const double expectedFine = law.contrabandFineCr(illegalValue);
    const double expectedPaid = std::min(credits, expectedFine);
    const double expectedUnpaid = expectedFine - expectedPaid;

    if (!nearly(res.fineCr, expectedFine, 1e-9)) {
      std::cerr << "[test_police_scan] enforceContraband fine mismatch\n";
      ++fails;
    }
    if (!nearly(res.paidCr, expectedPaid, 1e-9) || !nearly(res.unpaidCr, expectedUnpaid, 1e-9)) {
      std::cerr << "[test_police_scan] enforceContraband paid/unpaid mismatch\n";
      ++fails;
    }

    if (!nearly(res.creditsAfter, credits - expectedPaid, 1e-9)) {
      std::cerr << "[test_police_scan] enforceContraband creditsAfter mismatch\n";
      ++fails;
    }

    // Confiscate min(have, scanned)
    if (!nearly(res.cargoAfter[0], 2.0, 1e-9)) {
      std::cerr << "[test_police_scan] enforceContraband cargoAfter mismatch\n";
      ++fails;
    }

    if (!(res.policeHeatDelta >= 0.0)) {
      std::cerr << "[test_police_scan] expected non-negative policeHeatDelta\n";
      ++fails;
    }

    if (!nearly(res.bountyAddedCr, std::max(0.0, expectedUnpaid), 1e-9)) {
      std::cerr << "[test_police_scan] expected bountyAddedCr==unpaid\n";
      ++fails;
    }
  }

  // --- Evade math ------------------------------------------------------------
  {
    LawProfile law{};
    law.scanStrictness = 1.2;
    law.evadeRepMult = 1.4;
    law.repBase = 3.0;
    law.repDiv = 3500.0;
    law.repMin = 2.0;
    law.repMax = 9.0;

    const double fine = 900.0;
    const double illegalValue = 4000.0;

    const auto resA = evadeContraband(law, fine, illegalValue, true);
    if (!nearly(resA.bountyAddedCr, fine, 1e-9)) {
      std::cerr << "[test_police_scan] evadeContraband bounty mismatch\n";
      ++fails;
    }
    if (!(resA.shipHeatDelta > 0.0)) {
      std::cerr << "[test_police_scan] expected shipHeatDelta > 0\n";
      ++fails;
    }

    const auto resB = evadeContraband(law, fine, illegalValue, false);
    if (!(nearly(resB.policeHeatDelta, 0.0) && nearly(resB.policeAlertSeconds, 0.0))) {
      std::cerr << "[test_police_scan] expected no local escalation when escalateLocal=false\n";
      ++fails;
    }
  }

  return fails;
}
