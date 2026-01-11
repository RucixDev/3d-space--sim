#include "stellar/sim/CargoJettisonPlanner.h"

#include <cmath>
#include <iostream>

static bool near(double a, double b, double eps) {
  return std::abs(a - b) <= eps;
}

static int findLineUnits(const stellar::sim::CargoJettisonPlan& plan, stellar::econ::CommodityId cid) {
  for (const auto& ln : plan.lines) {
    if (ln.commodity == cid) return (int)std::lround(ln.units);
  }
  return 0;
}

int test_cargo_jettison_planner() {
  int fails = 0;

  using stellar::econ::CommodityId;

  // --- Avoid reserved mission cargo when possible ---
  {
    std::array<double, stellar::econ::kCommodityCount> cargo{};
    std::array<double, stellar::econ::kCommodityCount> reserved{};
    cargo.fill(0.0);
    reserved.fill(0.0);

    cargo[(std::size_t)CommodityId::Food] = 100.0;
    reserved[(std::size_t)CommodityId::Food] = 50.0;

    const auto plan = stellar::sim::planCargoJettisonForValue(cargo, reserved, /*requiredValueCr=*/600.0,
                                                             /*allowUsingReserved=*/false);

    if (!plan.success) {
      std::cerr << "[test_cargo_jettison_planner] expected success without reserved cargo\n";
      ++fails;
    }
    if (plan.usedReserved) {
      std::cerr << "[test_cargo_jettison_planner] expected usedReserved=false\n";
      ++fails;
    }
    if (!near(plan.plannedValueCr, 600.0, 1e-6)) {
      std::cerr << "[test_cargo_jettison_planner] plannedValueCr mismatch: got=" << plan.plannedValueCr << "\n";
      ++fails;
    }
    const int foodUnits = findLineUnits(plan, CommodityId::Food);
    if (foodUnits != 50) {
      std::cerr << "[test_cargo_jettison_planner] expected 50 Food units, got " << foodUnits << "\n";
      ++fails;
    }
  }

  // --- Insufficient free cargo: returns dump-everything plan (success=false) ---
  {
    std::array<double, stellar::econ::kCommodityCount> cargo{};
    std::array<double, stellar::econ::kCommodityCount> reserved{};
    cargo.fill(0.0);
    reserved.fill(0.0);

    cargo[(std::size_t)CommodityId::Food] = 100.0;
    reserved[(std::size_t)CommodityId::Food] = 50.0;

    const auto plan = stellar::sim::planCargoJettisonForValue(cargo, reserved, /*requiredValueCr=*/900.0,
                                                             /*allowUsingReserved=*/false);

    if (plan.success) {
      std::cerr << "[test_cargo_jettison_planner] expected failure without reserved cargo\n";
      ++fails;
    }
    if (!near(plan.plannedValueCr, 600.0, 1e-6)) {
      std::cerr << "[test_cargo_jettison_planner] expected plannedValueCr==600 (dump all free), got "
                << plan.plannedValueCr << "\n";
      ++fails;
    }
    const int foodUnits = findLineUnits(plan, CommodityId::Food);
    if (foodUnits != 50) {
      std::cerr << "[test_cargo_jettison_planner] expected dump 50 Food units, got " << foodUnits << "\n";
      ++fails;
    }
  }

  // --- Allow reserved cargo: succeeds and marks usedReserved=true ---
  {
    std::array<double, stellar::econ::kCommodityCount> cargo{};
    std::array<double, stellar::econ::kCommodityCount> reserved{};
    cargo.fill(0.0);
    reserved.fill(0.0);

    cargo[(std::size_t)CommodityId::Food] = 100.0;
    reserved[(std::size_t)CommodityId::Food] = 50.0;

    const auto plan = stellar::sim::planCargoJettisonForValue(cargo, reserved, /*requiredValueCr=*/900.0,
                                                             /*allowUsingReserved=*/true);

    if (!plan.success) {
      std::cerr << "[test_cargo_jettison_planner] expected success when reserved cargo allowed\n";
      ++fails;
    }
    if (!plan.usedReserved) {
      std::cerr << "[test_cargo_jettison_planner] expected usedReserved=true\n";
      ++fails;
    }
    if (!near(plan.plannedValueCr, 900.0, 1e-6)) {
      std::cerr << "[test_cargo_jettison_planner] plannedValueCr mismatch: got=" << plan.plannedValueCr << "\n";
      ++fails;
    }
    const int foodUnits = findLineUnits(plan, CommodityId::Food);
    if (foodUnits != 75) {
      std::cerr << "[test_cargo_jettison_planner] expected 75 Food units, got " << foodUnits << "\n";
      ++fails;
    }
  }

  // --- Minimizes overpay first (even if that increases mass) ---
  {
    std::array<double, stellar::econ::kCommodityCount> cargo{};
    std::array<double, stellar::econ::kCommodityCount> reserved{};
    cargo.fill(0.0);
    reserved.fill(0.0);

    cargo[(std::size_t)CommodityId::Fuel] = 10.0;       // 35 cr each
    cargo[(std::size_t)CommodityId::Stimulants] = 2.0;  // 165 cr each

    const auto plan = stellar::sim::planCargoJettisonForValue(cargo, reserved, /*requiredValueCr=*/205.0,
                                                             /*allowUsingReserved=*/false);

    if (!plan.success) {
      std::cerr << "[test_cargo_jettison_planner] expected success\n";
      ++fails;
    }
    // Best achievable is 210 via 6 Fuel.
    if (!near(plan.plannedValueCr, 210.0, 1e-6)) {
      std::cerr << "[test_cargo_jettison_planner] expected plannedValueCr==210, got " << plan.plannedValueCr << "\n";
      ++fails;
    }
    const int fuelUnits = findLineUnits(plan, CommodityId::Fuel);
    if (fuelUnits != 6) {
      std::cerr << "[test_cargo_jettison_planner] expected 6 Fuel units, got " << fuelUnits << "\n";
      ++fails;
    }
    const int stimUnits = findLineUnits(plan, CommodityId::Stimulants);
    if (stimUnits != 0) {
      std::cerr << "[test_cargo_jettison_planner] expected 0 Stimulants units, got " << stimUnits << "\n";
      ++fails;
    }
  }

  // --- For the same minimal overpay, chooses lower-mass combination ---
  {
    std::array<double, stellar::econ::kCommodityCount> cargo{};
    std::array<double, stellar::econ::kCommodityCount> reserved{};
    cargo.fill(0.0);
    reserved.fill(0.0);

    cargo[(std::size_t)CommodityId::Food] = 100.0;  // 12 cr
    cargo[(std::size_t)CommodityId::Water] = 100.0; // 6 cr

    const auto plan = stellar::sim::planCargoJettisonForValue(cargo, reserved, /*requiredValueCr=*/200.0,
                                                             /*allowUsingReserved=*/false);

    if (!plan.success) {
      std::cerr << "[test_cargo_jettison_planner] expected success\n";
      ++fails;
    }
    // Smallest reachable >=200 is 204.
    if (!near(plan.plannedValueCr, 204.0, 1e-6)) {
      std::cerr << "[test_cargo_jettison_planner] expected plannedValueCr==204, got " << plan.plannedValueCr << "\n";
      ++fails;
    }
    // 17 Food == 204 (mass 17) is lighter than 34 Water == 204 (mass 34).
    const int foodUnits = findLineUnits(plan, CommodityId::Food);
    if (foodUnits != 17) {
      std::cerr << "[test_cargo_jettison_planner] expected 17 Food units, got " << foodUnits << "\n";
      ++fails;
    }
  }

  return fails;
}
