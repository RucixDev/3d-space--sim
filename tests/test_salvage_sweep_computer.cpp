#include "stellar/sim/SalvageSweep.h"

#include <cmath>
#include <iostream>
#include <vector>

using namespace stellar;

namespace {
math::Vec3d V(double x, double y, double z) { return math::Vec3d{x, y, z}; }
}

int test_salvage_sweep_computer() {
  int fails = 0;

  // 1) Selection: prefers mission pods when missionBonus outweighs a slightly better value term.
  {
    sim::SalvageSweepParams p{};
    p.missionBonus = 0.65;
    p.valueScaleCr = 100.0;
    p.valueWeight = 1.0;
    p.distScaleKm = 100.0;
    p.distWeight = 1.0;
    p.switchHysteresis = 0.25;
    p.skipSeconds = 10.0;

    sim::SalvageSweepComputer comp;
    comp.setParams(p);

    const math::Vec3d shipPos = V(0.0, 0.0, 0.0);
    std::vector<sim::SalvagePodCandidate> cands;
    cands.push_back(sim::SalvagePodCandidate{.id = 1,
                                             .commodity = econ::CommodityId::Water,
                                             .units = 1.0,
                                             .posKm = shipPos,
                                             .velKmS = V(0.0, 0.0, 0.0),
                                             .missionId = 0,
                                             .estimatedValueCr = 100.0});
    cands.push_back(sim::SalvagePodCandidate{.id = 2,
                                             .commodity = econ::CommodityId::Water,
                                             .units = 1.0,
                                             .posKm = shipPos,
                                             .velKmS = V(0.0, 0.0, 0.0),
                                             .missionId = 123,
                                             .estimatedValueCr = 0.0});

    const auto st = comp.update(/*timeSec=*/0.0, shipPos, cands);
    if (!st.hasTarget || st.targetId != 2 || st.targetCandidateIndex != 1) {
      std::cerr << "[test_salvage_sweep_computer] expected mission candidate id=2 to win; got hasTarget="
                << st.hasTarget << " id=" << st.targetId << " idx=" << st.targetCandidateIndex << "\n";
      ++fails;
    }
  }

  // 2) Skip: temporarily excludes the current target.
  {
    sim::SalvageSweepParams p{};
    p.missionBonus = 0.0;
    p.valueScaleCr = 100.0;
    p.valueWeight = 1.0;
    p.distScaleKm = 100.0;
    p.distWeight = 0.0;
    p.switchHysteresis = 0.0;
    p.skipSeconds = 10.0;

    sim::SalvageSweepComputer comp;
    comp.setParams(p);
    const math::Vec3d shipPos = V(0.0, 0.0, 0.0);

    std::vector<sim::SalvagePodCandidate> cands;
    cands.push_back(sim::SalvagePodCandidate{.id = 10,
                                             .commodity = econ::CommodityId::Water,
                                             .units = 1.0,
                                             .posKm = shipPos,
                                             .velKmS = V(0.0, 0.0, 0.0),
                                             .missionId = 0,
                                             .estimatedValueCr = 200.0});
    cands.push_back(sim::SalvagePodCandidate{.id = 11,
                                             .commodity = econ::CommodityId::Water,
                                             .units = 1.0,
                                             .posKm = shipPos,
                                             .velKmS = V(0.0, 0.0, 0.0),
                                             .missionId = 0,
                                             .estimatedValueCr = 100.0});

    auto st = comp.update(0.0, shipPos, cands);
    if (!st.hasTarget || st.targetId != 10) {
      std::cerr << "[test_salvage_sweep_computer] expected id=10 as initial best, got hasTarget=" << st.hasTarget
                << " id=" << st.targetId << "\n";
      ++fails;
    }

    comp.skip(/*id=*/10, /*timeSec=*/0.0);
    st = comp.update(1.0, shipPos, cands);
    if (!st.hasTarget || st.targetId != 11) {
      std::cerr << "[test_salvage_sweep_computer] expected skip to force id=11, got hasTarget=" << st.hasTarget
                << " id=" << st.targetId << "\n";
      ++fails;
    }

    st = comp.update(11.0, shipPos, cands);
    if (!st.hasTarget || st.targetId != 10) {
      std::cerr << "[test_salvage_sweep_computer] expected id=10 to be eligible after skip, got hasTarget="
                << st.hasTarget << " id=" << st.targetId << "\n";
      ++fails;
    }
  }

  // 3) Hysteresis: discourages target thrash.
  {
    sim::SalvageSweepParams p{};
    p.missionBonus = 0.0;
    p.valueScaleCr = 100.0;
    p.valueWeight = 1.0;
    p.distScaleKm = 100.0;
    p.distWeight = 0.0;
    p.switchHysteresis = 0.25;
    p.skipSeconds = 10.0;

    sim::SalvageSweepComputer comp;
    comp.setParams(p);
    const math::Vec3d shipPos = V(0.0, 0.0, 0.0);

    std::vector<sim::SalvagePodCandidate> cands;
    cands.push_back(sim::SalvagePodCandidate{.id = 20,
                                             .commodity = econ::CommodityId::Water,
                                             .units = 1.0,
                                             .posKm = shipPos,
                                             .velKmS = V(0.0, 0.0, 0.0),
                                             .missionId = 0,
                                             .estimatedValueCr = 100.0});
    cands.push_back(sim::SalvagePodCandidate{.id = 21,
                                             .commodity = econ::CommodityId::Water,
                                             .units = 1.0,
                                             .posKm = shipPos,
                                             .velKmS = V(0.0, 0.0, 0.0),
                                             .missionId = 0,
                                             .estimatedValueCr = 105.0});

    auto st = comp.update(0.0, shipPos, cands);
    if (!st.hasTarget || st.targetId != 21) {
      std::cerr << "[test_salvage_sweep_computer] expected id=21 as initial best, got hasTarget=" << st.hasTarget
                << " id=" << st.targetId << "\n";
      ++fails;
    }

    // Slightly improve the runner-up; should still keep current target due to hysteresis.
    cands[0].estimatedValueCr = 110.0;
    st = comp.update(1.0, shipPos, cands);
    if (!st.hasTarget || st.targetId != 21) {
      std::cerr << "[test_salvage_sweep_computer] expected hysteresis to hold id=21, got hasTarget=" << st.hasTarget
                << " id=" << st.targetId << "\n";
      ++fails;
    }

    // Big improvement should overcome hysteresis and flip.
    cands[0].estimatedValueCr = 200.0;
    st = comp.update(2.0, shipPos, cands);
    if (!st.hasTarget || st.targetId != 20) {
      std::cerr << "[test_salvage_sweep_computer] expected hysteresis to allow switch to id=20, got hasTarget="
                << st.hasTarget << " id=" << st.targetId << "\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_salvage_sweep_computer] pass\n";
  return fails;
}
