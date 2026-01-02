#include "stellar/sim/SecurityModel.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approxEq(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

int test_security_model() {
  int fails = 0;

  // Controlling faction and contestedness are deterministic.
  {
    sim::StarSystem sys{};
    {
      sim::Station a{};
      a.id = 1;
      a.name = "A";
      a.type = econ::StationType::Outpost;
      a.factionId = 2;
      sys.stations.push_back(a);
    }
    {
      sim::Station b{};
      b.id = 2;
      b.name = "B";
      b.type = econ::StationType::Outpost;
      b.factionId = 2;
      sys.stations.push_back(b);
    }
    {
      sim::Station c{};
      c.id = 3;
      c.name = "C";
      c.type = econ::StationType::Outpost;
      c.factionId = 5;
      sys.stations.push_back(c);
    }

    const auto ctrlA = sim::computeSystemControl(sys);
    const auto ctrlB = sim::computeSystemControl(sys);
    if (ctrlA.controllingFactionId != 2 || ctrlB.controllingFactionId != 2) {
      std::cerr << "[test_security_model] expected controlling faction 2. got="
                << ctrlA.controllingFactionId << "/" << ctrlB.controllingFactionId << "\n";
      ++fails;
    }
    if (!approxEq(ctrlA.controlFrac, ctrlB.controlFrac) || !approxEq(ctrlA.contest01, ctrlB.contest01)) {
      std::cerr << "[test_security_model] expected deterministic control fractions.\n";
      ++fails;
    }
    if (ctrlA.controlFrac <= 0.6 || ctrlA.controlFrac >= 0.8) {
      std::cerr << "[test_security_model] expected controlFrac ~ 0.666. got=" << ctrlA.controlFrac << "\n";
      ++fails;
    }
    if (ctrlA.contest01 <= 0.2 || ctrlA.contest01 >= 0.5) {
      std::cerr << "[test_security_model] expected contest01 ~ 0.333. got=" << ctrlA.contest01 << "\n";
      ++fails;
    }
  }

  // Tie-breaks are stable: if weights are equal, lowest faction id wins.
  {
    sim::StarSystem sys{};
    sim::Station a{};
    a.id = 1;
    a.name = "A";
    a.type = econ::StationType::Outpost;
    a.factionId = 7;
    sys.stations.push_back(a);

    sim::Station b{};
    b.id = 2;
    b.name = "B";
    b.type = econ::StationType::Outpost;
    b.factionId = 3;
    sys.stations.push_back(b);

    const auto ctrl = sim::computeSystemControl(sys);
    if (ctrl.controllingFactionId != 3) {
      std::cerr << "[test_security_model] tie-break mismatch. expected=3 got="
                << ctrl.controllingFactionId << "\n";
      ++fails;
    }
  }

  // Baseline for fully independent systems is centered.
  {
    sim::StarSystem sys{};
    sim::Station a{};
    a.id = 1;
    a.name = "Ind";
    a.type = econ::StationType::Outpost;
    a.factionId = 0;
    sys.stations.push_back(a);

    const auto sp = sim::systemSecurityProfile(1337u, sys);
    if (sp.controllingFactionId != 0) {
      std::cerr << "[test_security_model] expected controllingFactionId=0 got="
                << sp.controllingFactionId << "\n";
      ++fails;
    }
    if (sp.security01 < 0.30 || sp.security01 > 0.70 ||
        sp.piracy01 < 0.30 || sp.piracy01 > 0.70 ||
        sp.traffic01 < 0.30 || sp.traffic01 > 0.70) {
      std::cerr << "[test_security_model] expected independent baseline around 0.5. got sec="
                << sp.security01 << " piracy=" << sp.piracy01 << " traffic=" << sp.traffic01 << "\n";
      ++fails;
    }
  }

  // Contestedness reduces security (all else equal).
  {
    sim::StarSystem sysA{};
    sim::StarSystem sysB{};

    // SysA: fully controlled by faction 2.
    {
      sim::Station st{};
      st.id = 1;
      st.type = econ::StationType::Outpost;
      st.factionId = 2;
      sysA.stations.push_back(st);
    }
    {
      sim::Station st{};
      st.id = 2;
      st.type = econ::StationType::Outpost;
      st.factionId = 2;
      sysA.stations.push_back(st);
    }

    // SysB: same controlling faction but contested.
    {
      sim::Station st{};
      st.id = 1;
      st.type = econ::StationType::Outpost;
      st.factionId = 2;
      sysB.stations.push_back(st);
    }
    {
      sim::Station st{};
      st.id = 2;
      st.type = econ::StationType::Outpost;
      st.factionId = 2;
      sysB.stations.push_back(st);
    }
    {
      sim::Station st{};
      st.id = 3;
      st.type = econ::StationType::Outpost;
      st.factionId = 9;
      sysB.stations.push_back(st);
    }

    const auto a = sim::systemSecurityProfile(9001u, sysA);
    const auto b = sim::systemSecurityProfile(9001u, sysB);
    if (a.controllingFactionId != b.controllingFactionId) {
      std::cerr << "[test_security_model] expected same controlling faction.\n";
      ++fails;
    }
    if (!(b.contest01 > a.contest01 + 1e-9)) {
      std::cerr << "[test_security_model] expected contested system to have higher contest01.\n";
      ++fails;
    }
    if (!(b.security01 < a.security01 + 1e-9)) {
      std::cerr << "[test_security_model] expected contested system to have lower security.\n";
      ++fails;
    }
  }

  if (fails == 0) {
    std::cout << "[test_security_model] PASS\n";
  }
  return fails;
}
