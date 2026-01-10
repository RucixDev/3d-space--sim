#include "stellar/sim/SystemSecurityDynamics.h"

#include "stellar/sim/SecurityModel.h"

#include <cmath>
#include <iostream>

using namespace stellar;

static bool approxEq(double a, double b, double eps = 1e-12) {
  return std::abs(a - b) <= eps;
}

int test_system_security_dynamics() {
  int fails = 0;

  // decayFactorDays uses half-life semantics.
  {
    const double f0 = sim::decayFactorDays(0.0, 10.0);
    if (!approxEq(f0, 1.0)) {
      std::cerr << "[test_system_security_dynamics] expected dt=0 factor=1 got=" << f0 << "\n";
      ++fails;
    }

    const double fHalf = sim::decayFactorDays(5.0, 5.0);
    if (fHalf < 0.49 || fHalf > 0.51) {
      std::cerr << "[test_system_security_dynamics] expected half-life factor~0.5 got=" << fHalf << "\n";
      ++fails;
    }

    const double fDead = sim::decayFactorDays(1.0, 0.0);
    if (!approxEq(fDead, 0.0)) {
      std::cerr << "[test_system_security_dynamics] expected halfLife<=0 => factor=0 got=" << fDead << "\n";
      ++fails;
    }
  }

  // Impulses decay-to-now then clamp.
  {
    sim::SystemSecurityDynamicsParams p{};
    p.securityHalfLifeDays = 10.0;
    p.piracyHalfLifeDays = 10.0;
    p.trafficHalfLifeDays = 10.0;
    p.maxAbsDelta = 0.25;

    sim::SystemSecurityDeltaState st{};
    st.systemId = 42;
    st.securityDelta = 0.10;
    st.piracyDelta = -0.10;
    st.trafficDelta = 0.10;
    st.lastUpdateDay = 0.0;

    // After 10 days with 10-day half-life, delta should halve.
    sim::applySystemSecurityImpulse(st, 10.0, 0.0, 0.0, 0.0, p);
    if (st.securityDelta < 0.049 || st.securityDelta > 0.051) {
      std::cerr << "[test_system_security_dynamics] expected decayed securityDelta~0.05 got=" << st.securityDelta << "\n";
      ++fails;
    }

    // Large impulse should clamp.
    sim::applySystemSecurityImpulse(st, 10.0, 10.0, -10.0, 10.0, p);
    if (st.securityDelta > 0.251 || st.securityDelta < 0.249) {
      std::cerr << "[test_system_security_dynamics] expected clamped securityDelta~0.25 got=" << st.securityDelta << "\n";
      ++fails;
    }
    if (st.piracyDelta < -0.251 || st.piracyDelta > -0.249) {
      std::cerr << "[test_system_security_dynamics] expected clamped piracyDelta~-0.25 got=" << st.piracyDelta << "\n";
      ++fails;
    }
  }

  // Applying deltas clamps the effective profile to [0,1].
  {
    sim::SystemSecurityDynamicsParams p{};
    p.securityHalfLifeDays = 1000.0;
    p.piracyHalfLifeDays = 1000.0;
    p.trafficHalfLifeDays = 1000.0;
    p.maxAbsDelta = 1.0;

    sim::SystemSecurityProfile base{};
    base.security01 = 0.95;
    base.piracy01 = 0.05;
    base.traffic01 = 0.50;

    sim::SystemSecurityDeltaState st{};
    st.systemId = 1;
    st.securityDelta = 0.20;
    st.piracyDelta = -0.20;
    st.trafficDelta = 0.30;
    st.lastUpdateDay = 0.0;

    const auto out = sim::applySystemSecurityDelta(base, st, 0.0, p);
    if (!approxEq(out.security01, 1.0)) {
      std::cerr << "[test_system_security_dynamics] expected security01 clamped to 1.0 got=" << out.security01 << "\n";
      ++fails;
    }
    if (!approxEq(out.piracy01, 0.0)) {
      std::cerr << "[test_system_security_dynamics] expected piracy01 clamped to 0.0 got=" << out.piracy01 << "\n";
      ++fails;
    }
    if (out.traffic01 < 0.79 || out.traffic01 > 0.81) {
      std::cerr << "[test_system_security_dynamics] expected traffic01~0.8 got=" << out.traffic01 << "\n";
      ++fails;
    }
  }

  // Negligible deltas can be pruned after enough decay.
  {
    sim::SystemSecurityDynamicsParams p{};
    p.securityHalfLifeDays = 1.0;
    p.piracyHalfLifeDays = 1.0;
    p.trafficHalfLifeDays = 1.0;
    p.negligibleAbs = 1e-4;

    sim::SystemSecurityDeltaState st{};
    st.systemId = 123;
    st.securityDelta = 0.01;
    st.piracyDelta = 0.0;
    st.trafficDelta = 0.0;
    st.lastUpdateDay = 0.0;

    // After 20 half-lives, 0.01 -> 0.01/2^20 ~ 9.5e-9 (below negligibleAbs).
    if (!sim::isSystemSecurityDeltaNegligible(st, 20.0, p)) {
      std::cerr << "[test_system_security_dynamics] expected delta to be negligible after decay\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_system_security_dynamics] pass\n";
  return fails;
}
