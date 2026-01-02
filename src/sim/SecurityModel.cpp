#include "stellar/sim/SecurityModel.h"

#include "stellar/econ/Economy.h"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace stellar::sim {

static double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

static double stationWeight(econ::StationType t) {
  using econ::StationType;
  // Very lightweight weighting: stations that plausibly anchor more trade/traffic
  // influence control slightly more.
  switch (t) {
    case StationType::TradeHub: return 1.35;
    case StationType::Shipyard: return 1.30;
    case StationType::Research: return 1.20;
    case StationType::Industrial: return 1.15;
    case StationType::Refinery: return 1.10;
    case StationType::Mining: return 1.05;
    case StationType::Agricultural: return 1.05;
    case StationType::Outpost: return 1.0;
    default: return 1.0;
  }
}

SystemControl computeSystemControl(const StarSystem& sys) {
  SystemControl out{};

  if (sys.stations.empty()) {
    // No hard ownership signals.
    out.controllingFactionId = 0;
    out.controlFrac = 0.0;
    out.contest01 = 0.0;
    return out;
  }

  std::unordered_map<core::u32, double> weights;
  weights.reserve(sys.stations.size());
  double total = 0.0;

  for (const auto& st : sys.stations) {
    const double w = stationWeight(st.type);
    weights[st.factionId] += w;
    total += w;
  }

  if (total <= 1e-12) {
    out.controllingFactionId = 0;
    out.controlFrac = 0.0;
    out.contest01 = 0.0;
    return out;
  }

  // Find controlling faction (max weight; deterministic tie-break by id).
  core::u32 bestId = 0;
  double bestW = -1.0;
  for (const auto& kv : weights) {
    const core::u32 fid = kv.first;
    const double w = kv.second;
    if (w > bestW + 1e-12 || (std::abs(w - bestW) <= 1e-12 && fid < bestId)) {
      bestW = w;
      bestId = fid;
    }
  }

  out.controllingFactionId = bestId;
  out.controlFrac = std::clamp(bestW / total, 0.0, 1.0);
  out.contest01 = clamp01(1.0 - out.controlFrac);

  // Sorted breakdown for tools/UI.
  out.factionWeights.reserve(weights.size());
  for (const auto& kv : weights) {
    out.factionWeights.emplace_back(kv.first, kv.second);
  }
  std::sort(out.factionWeights.begin(), out.factionWeights.end(), [](const auto& a, const auto& b) {
    if (a.second != b.second) return a.second > b.second;
    return a.first < b.first;
  });

  return out;
}

SystemSecurityProfile systemSecurityProfile(core::u64 universeSeed, const StarSystem& sys) {
  SystemSecurityProfile out{};

  const SystemControl ctrl = computeSystemControl(sys);
  out.controllingFactionId = ctrl.controllingFactionId;
  out.controlFrac = ctrl.controlFrac;
  out.contest01 = ctrl.contest01;
  out.traits = factionProfile(universeSeed, out.controllingFactionId);

  const double auth = clamp01(out.traits.authority);
  const double corr = clamp01(out.traits.corruption);
  const double wealth = clamp01(out.traits.wealth);
  const double stab = clamp01(out.traits.stability);
  const double mil = clamp01(out.traits.militarism);
  (void)corr;

  // Security effectiveness: authority/stability/militarism help, contestedness hurts.
  out.security01 = clamp01(0.35 * auth + 0.40 * stab + 0.25 * mil - 0.30 * out.contest01);

  // Piracy pressure: instability is the primary driver. Wealth and contestedness
  // attract opportunists, militarism suppresses.
  out.piracy01 = clamp01(0.25
                         + 0.55 * (1.0 - stab)
                         + 0.20 * wealth
                         + 0.20 * out.contest01
                         - 0.25 * mil);

  // Civilian/trader traffic: wealth + stability increase throughput, piracy and
  // contestedness reduce it.
  out.traffic01 = clamp01(0.05
                          + 0.45 * wealth
                          + 0.25 * stab
                          + 0.20 * (1.0 - out.piracy01)
                          - 0.25 * out.contest01);

  return out;
}

} // namespace stellar::sim
