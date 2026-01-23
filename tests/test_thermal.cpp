#include "stellar/sim/ThermalSystem.h"

#include <cmath>
#include <iostream>

static bool nearly(double a, double b, double eps = 1e-9) {
  return std::abs(a - b) <= eps;
}

int test_thermal() {
  int fails = 0;

  using namespace stellar::sim;

  // Baseline: no heating, undocked, nominal cooling.
  {
    ThermalInputs in{};
    in.dtReal = 2.0;
    in.docked = false;
    in.heatCoolRate = 10.0;

    const auto r = stepThermal(/*heat*/50.0, in);
    const double expected = 50.0 + (0.0 - 10.0) * 2.0;
    if (!nearly(r.heat, expected)) {
      std::cerr << "[test_thermal] baseline cooling mismatch. got=" << r.heat
                << " expected=" << expected << "\n";
      ++fails;
    }
    if (r.hullDamage != 0.0) {
      std::cerr << "[test_thermal] baseline cooling should not cause hull damage\n";
      ++fails;
    }
  }

  // Docked cooling should be 2x the undocked base (by default params).
  {
    ThermalInputs in{};
    in.dtReal = 1.0;
    in.docked = true;
    in.heatCoolRate = 10.0;

    const auto r = stepThermal(/*heat*/40.0, in);
    const double expected = 40.0 + (0.0 - 20.0) * 1.0;
    if (!nearly(r.heat, expected)) {
      std::cerr << "[test_thermal] docked cooling mismatch. got=" << r.heat
                << " expected=" << expected << "\n";
      ++fails;
    }
  }

  // Impulse heat should apply before the integration step.
  {
    ThermalInputs in{};
    in.dtReal = 1.0;
    in.heatImpulse = 15.0;
    in.docked = false;
    in.heatCoolRate = 10.0;

    const auto r = stepThermal(/*heat*/10.0, in);
    const double expected = (10.0 + 15.0) + (0.0 - 10.0) * 1.0;
    if (!nearly(r.heat, expected)) {
      std::cerr << "[test_thermal] impulse mismatch. got=" << r.heat
                << " expected=" << expected << "\n";
      ++fails;
    }
  }

  // Heating sources: boost + supercruise + FSD.
  {
    ThermalInputs in{};
    in.dtReal = 1.0;
    in.docked = false;
    in.heatCoolRate = 10.0;
    in.boostAppliedFrac = 0.5;
    in.supercruiseActive = true;
    in.fsd = ThermalFsdState::Charging;

    // heatIn = 18*0.5 + 6 + 12 = 27
    const auto r = stepThermal(/*heat*/0.0, in);
    const double expected = 0.0 + (27.0 - 10.0) * 1.0;
    if (!nearly(r.heat, expected)) {
      std::cerr << "[test_thermal] heating sources mismatch. got=" << r.heat
                << " expected=" << expected << "\n";
      ++fails;
    }
  }

  // Overheat damage should scale with hullMax and dt.
  {
    ThermalInputs in{};
    in.dtReal = 2.0;
    in.docked = false;
    in.heatCoolRate = 0.0; // keep heat stable
    in.hullMax = 100.0;

    const auto r = stepThermal(/*heat*/110.0, in);
    // over=10, rate=0.0008, hullMax=100, dt=2 => 1.6
    const double expectedDamage = 10.0 * 0.0008 * 100.0 * 2.0;
    if (!nearly(r.hullDamage, expectedDamage, 1e-6)) {
      std::cerr << "[test_thermal] overheat damage mismatch. got=" << r.hullDamage
                << " expected=" << expectedDamage << "\n";
      ++fails;
    }
    if (!r.overheated) {
      std::cerr << "[test_thermal] expected overheated=true\n";
      ++fails;
    }
  }

  // Silent running: adds baseline heat and reduces cooling.
  {
    ThermalInputs in{};
    in.dtReal = 1.0;
    in.docked = false;
    in.silentRunning = true;
    in.heatCoolRate = 10.0;

    ThermalParams p{};
    p.baseCoolUndocked = 10.0;
    p.referenceCoolRate = 10.0;
    p.silentCoolMult = 0.25;
    p.heatPerSilentSec = 5.0;

    // Disable other continuous heating sources for this test.
    p.heatPerBoostSec = 0.0;
    p.heatPerSupercruiseSec = 0.0;
    p.heatPerFsdChargeSec = 0.0;
    p.heatPerFsdJumpSec = 0.0;

    const auto r = stepThermal(/*heat*/50.0, in, p);
    const double expected = 50.0 + (5.0 - (10.0 * 0.25)) * 1.0;
    if (!nearly(r.heat, expected)) {
      std::cerr << "[test_thermal] silent running mismatch. got=" << r.heat
                << " expected=" << expected << "\n";
      ++fails;
    }
    if (!nearly(r.heatInRate, 5.0)) {
      std::cerr << "[test_thermal] silent running heatInRate mismatch. got=" << r.heatInRate
                << " expected=5.0\n";
      ++fails;
    }
    if (!nearly(r.coolRate, 2.5)) {
      std::cerr << "[test_thermal] silent running coolRate mismatch. got=" << r.coolRate
                << " expected=2.5\n";
      ++fails;
    }
  }

  return fails;
}
