#include "stellar/sim/SystemSecurityDynamics.h"

#include "stellar/sim/SecurityModel.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {
namespace {

constexpr double kLn2 = 0.693147180559945309417232121458176568;

static double clampAbs(double v, double maxAbs) {
  const double m = std::max(0.0, maxAbs);
  return std::clamp(v, -m, m);
}

static double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

static double safeDt(double nowDays, double lastUpdateDay) {
  // Defensive: allow callers to pass slightly out-of-order timestamps.
  return std::max(0.0, nowDays - lastUpdateDay);
}

} // namespace

double decayFactorDays(double dtDays, double halfLifeDays) {
  if (dtDays <= 0.0) return 1.0;
  if (halfLifeDays <= 0.0) return 0.0;
  return std::exp(-kLn2 * dtDays / halfLifeDays);
}

SystemSecurityDeltaState decayedSystemSecurityDelta(const SystemSecurityDeltaState& st,
                                                    double nowDays,
                                                    const SystemSecurityDynamicsParams& params) {
  SystemSecurityDeltaState out = st;
  const double dt = safeDt(nowDays, st.lastUpdateDay);

  const double fSec = decayFactorDays(dt, params.securityHalfLifeDays);
  const double fPir = decayFactorDays(dt, params.piracyHalfLifeDays);
  const double fTrf = decayFactorDays(dt, params.trafficHalfLifeDays);

  out.securityDelta *= fSec;
  out.piracyDelta *= fPir;
  out.trafficDelta *= fTrf;
  out.lastUpdateDay = nowDays;
  return out;
}

void decaySystemSecurityDeltaInPlace(SystemSecurityDeltaState& st,
                                     double nowDays,
                                     const SystemSecurityDynamicsParams& params) {
  const double dt = safeDt(nowDays, st.lastUpdateDay);
  if (dt <= 0.0) {
    st.lastUpdateDay = nowDays;
    return;
  }

  st.securityDelta *= decayFactorDays(dt, params.securityHalfLifeDays);
  st.piracyDelta *= decayFactorDays(dt, params.piracyHalfLifeDays);
  st.trafficDelta *= decayFactorDays(dt, params.trafficHalfLifeDays);
  st.lastUpdateDay = nowDays;
}

void applySystemSecurityImpulse(SystemSecurityDeltaState& st,
                                double nowDays,
                                double dSecurity,
                                double dPiracy,
                                double dTraffic,
                                const SystemSecurityDynamicsParams& params) {
  decaySystemSecurityDeltaInPlace(st, nowDays, params);

  st.securityDelta = clampAbs(st.securityDelta + dSecurity, params.maxAbsDelta);
  st.piracyDelta = clampAbs(st.piracyDelta + dPiracy, params.maxAbsDelta);
  st.trafficDelta = clampAbs(st.trafficDelta + dTraffic, params.maxAbsDelta);
}

bool isSystemSecurityDeltaNegligible(const SystemSecurityDeltaState& st,
                                     double nowDays,
                                     const SystemSecurityDynamicsParams& params) {
  const auto d = decayedSystemSecurityDelta(st, nowDays, params);
  const double eps = std::max(0.0, params.negligibleAbs);
  return std::abs(d.securityDelta) <= eps && std::abs(d.piracyDelta) <= eps && std::abs(d.trafficDelta) <= eps;
}

SystemSecurityProfile applySystemSecurityDelta(const SystemSecurityProfile& base,
                                               const SystemSecurityDeltaState& st,
                                               double nowDays,
                                               const SystemSecurityDynamicsParams& params) {
  const auto d = decayedSystemSecurityDelta(st, nowDays, params);

  SystemSecurityProfile out = base;
  out.security01 = clamp01(out.security01 + d.securityDelta);
  out.piracy01 = clamp01(out.piracy01 + d.piracyDelta);
  out.traffic01 = clamp01(out.traffic01 + d.trafficDelta);
  return out;
}

} // namespace stellar::sim
