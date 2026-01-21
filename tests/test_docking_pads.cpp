#include "stellar/sim/DockingPads.h"

#include <iostream>

using namespace stellar;

int test_docking_pads() {
  int fails = 0;

  sim::Station st{};
  st.id = 42;
  st.type = econ::StationType::TradeHub;
  st.radiusKm = 6500.0;
  st.slotWidthKm = 5400.0;
  st.slotHeightKm = 4800.0;
  st.slotDepthKm = 9000.0;

  const core::u64 seed = 0xBADC0FFEEu;

  const int n = sim::stationDockingPadCount(st);
  if (n < 4) {
    std::cerr << "[test_docking_pads] expected at least 4 pads. got=" << n << "\n";
    ++fails;
  }

  // Determinism: same seed yields identical pads.
  for (int i = 1; i <= std::min(n, 12); ++i) {
    const auto a = sim::stationDockingPad(seed, st, (core::u16)i);
    const auto b = sim::stationDockingPad(seed, st, (core::u16)i);
    const math::Vec3d d = a.localPosKm - b.localPosKm;
    if (d.lengthSq() > 1e-12) {
      std::cerr << "[test_docking_pads] nondeterministic pad " << i << "\n";
      ++fails;
      break;
    }
  }

  // Pads should lie within the station hull box used by collision code.
  const double wx = st.radiusKm * 0.70;
  const double wy = st.radiusKm * 0.70;
  const double wz = st.radiusKm * 1.10;

  for (int i = 1; i <= n; ++i) {
    const auto p = sim::stationDockingPad(seed, st, (core::u16)i);
    if (std::abs(p.localPosKm.x) > wx + 1e-6 ||
        std::abs(p.localPosKm.y) > wy + 1e-6 ||
        std::abs(p.localPosKm.z) > wz + 1e-6) {
      std::cerr << "[test_docking_pads] pad out of bounds: " << i << " pos=("
                << p.localPosKm.x << "," << p.localPosKm.y << "," << p.localPosKm.z << ")\n";
      ++fails;
      break;
    }
  }

  // Nearest-pad query should recover the original pad id (within a small offset).
  if (n >= 3) {
    const auto p3 = sim::stationDockingPad(seed, st, 3);
    const math::Vec3d probe = p3.localPosKm + math::Vec3d{15.0, -9.0, 8.0};
    const auto idx = sim::nearestDockingPadIndex(seed, st, probe);
    if (idx != 3) {
      std::cerr << "[test_docking_pads] nearestDockingPadIndex mismatch. expected 3 got="
                << idx << "\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_docking_pads] PASS\n";
  }

  return fails;
}
