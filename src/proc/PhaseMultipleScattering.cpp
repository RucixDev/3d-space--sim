#include "stellar/proc/PhaseMultipleScattering.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace stellar::proc {

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;

// Trapezoidal integration over an evenly-sampled μ grid.
template <class F>
double trapMu(int n, double dMu, F&& fn) {
  if (n <= 0) return 0.0;
  if (n == 1) return fn(0);
  double sum = 0.5 * (fn(0) + fn(n - 1));
  for (int i = 1; i < n - 1; ++i) sum += fn(i);
  return sum * dMu;
}

// Renormalize a phase function so ∫ p dΩ = 1.
// p is sampled over μ∈[-1,1].
void renormalizePhase(std::vector<double>& p, const std::vector<double>& mu) {
  const int n = (int)p.size();
  if (n <= 1 || (int)mu.size() != n) return;
  const double dMu = (mu.back() - mu.front()) / (double)(n - 1);
  const double integral = 2.0 * kPi * trapMu(n, dMu, [&](int i) { return p[(std::size_t)i]; });
  if (!(integral > 0.0) || !std::isfinite(integral)) return;
  const double inv = 1.0 / integral;
  for (double& v : p) v *= inv;
}

} // namespace

SpectralMiePhaseResult generateMultipleScatteringPhase(
    const SpectralMiePhaseResult& single,
    const MultipleScatteringPhaseSettings& s) {
  SpectralMiePhaseResult out;
  const int N = single.angleSamples;
  if (N <= 0 || (int)single.mu.size() != N) return out;

  const int L = std::clamp(s.legendreOrder, 0, 256);
  const double omega = std::clamp(s.scatteringAlbedo, 0.0, 0.999);

  out.angleSamples = N;
  out.mu = single.mu;
  out.qSca = single.qSca;

  // Convert μ to double.
  std::vector<double> muD((std::size_t)N);
  for (int i = 0; i < N; ++i) muD[(std::size_t)i] = (double)single.mu[(std::size_t)i];

  const double dMu = (N > 1) ? (muD.back() - muD.front()) / (double)(N - 1) : 0.0;

  auto getChannel = [&](int c) -> const std::vector<float>& {
    if (c == 0) return single.phaseR;
    if (c == 1) return single.phaseG;
    return single.phaseB;
  };

  auto setChannel = [&](int c, std::vector<float>& dst, const std::vector<double>& src, double pMaxPhys) {
    dst.resize((std::size_t)N);
    const double invMax = (s.peakNormalize && pMaxPhys > 0.0) ? (1.0 / pMaxPhys) : 1.0;
    for (int i = 0; i < N; ++i) {
      double v = src[(std::size_t)i];
      if (s.peakNormalize) v *= invMax;
      dst[(std::size_t)i] = (float)std::clamp(v, 0.0, 1e9);
    }
  };

  // Compute per-channel.
  for (int c = 0; c < 3; ++c) {
    const auto& in = getChannel(c);
    if ((int)in.size() != N) {
      // Missing channel => isotropic.
      const double pIso = 1.0 / (4.0 * kPi);
      std::vector<double> iso((std::size_t)N, pIso);
      if (c == 0) setChannel(c, out.phaseR, iso, pIso);
      else if (c == 1) setChannel(c, out.phaseG, iso, pIso);
      else setChannel(c, out.phaseB, iso, pIso);
      out.pMax[(std::size_t)c] = (float)pIso;
      out.asymmetryG[(std::size_t)c] = 0.0f;
      continue;
    }

    // Recover physically normalized p(μ).
    const double scale = (s.inputWasPeakNormalized) ? (double)single.pMax[(std::size_t)c] : 1.0;
    std::vector<double> pPhys((std::size_t)N);
    for (int i = 0; i < N; ++i) pPhys[(std::size_t)i] = (double)in[(std::size_t)i] * scale;
    renormalizePhase(pPhys, muD);

    // Compute Legendre moments χ_l = 2π ∫ p(μ) P_l(μ) dμ.
    std::vector<double> chi((std::size_t)L + 1, 0.0);
    for (int l = 0; l <= L; ++l) {
      const double moment = 2.0 * kPi * trapMu(N, dMu, [&](int i) {
        const double mu = muD[(std::size_t)i];
        // Evaluate P_l(μ) on-the-fly with a small recurrence for stability.
        if (l == 0) return pPhys[(std::size_t)i];
        if (l == 1) return pPhys[(std::size_t)i] * mu;
        double Pm2 = 1.0;
        double Pm1 = mu;
        double Pl = 0.0;
        for (int k = 1; k < l; ++k) {
          const double kp1 = (double)k + 1.0;
          const double a = (2.0 * (double)k + 1.0) * mu * Pm1;
          const double b = (double)k * Pm2;
          Pl = (a - b) / kp1;
          Pm2 = Pm1;
          Pm1 = Pl;
        }
        return pPhys[(std::size_t)i] * Pm1;
      });
      chi[(std::size_t)l] = std::isfinite(moment) ? moment : 0.0;
    }

    // Enforce χ_0 ≈ 1.
    chi[0] = 1.0;

    // Analytic multiple-scattering moment sum.
    std::vector<double> chiMs((std::size_t)L + 1, 0.0);
    for (int l = 0; l <= L; ++l) {
      const double x = chi[(std::size_t)l];
      const double denom = std::max(1e-6, 1.0 - omega * x);
      double v = 0.0;
      if (s.mode == MultipleScatteringMode::TotalOrders) {
        // (1-ω) * χ / (1-ωχ)
        v = (1.0 - omega) * x / denom;
      } else {
        // Orders n>=2: (1-ω) * χ^2 / (1-ωχ)
        v = (1.0 - omega) * x * x / denom;
      }
      chiMs[(std::size_t)l] = std::isfinite(v) ? v : 0.0;
    }
    chiMs[0] = 1.0;

    // Reconstruct p_ms(μ) = (1/4π) Σ (2l+1)χ_ms P_l(μ).
    std::vector<double> pMs((std::size_t)N, 0.0);
    for (int i = 0; i < N; ++i) {
      const double mu = muD[(std::size_t)i];
      double sum = 0.0;

      // P_0 and P_1.
      double Pm2 = 1.0;
      sum += (1.0) * chiMs[0] * Pm2;
      if (L >= 1) {
        double Pm1 = mu;
        sum += (3.0) * chiMs[1] * Pm1;
        // Higher orders.
        for (int l = 1; l < L; ++l) {
          const double lp1 = (double)l + 1.0;
          const double a = (2.0 * (double)l + 1.0) * mu * Pm1;
          const double b = (double)l * Pm2;
          const double Pl = (a - b) / lp1;
          const double w = (2.0 * lp1 + 1.0);
          sum += w * chiMs[(std::size_t)l + 1] * Pl;
          Pm2 = Pm1;
          Pm1 = Pl;
        }
      }

      pMs[(std::size_t)i] = sum / (4.0 * kPi);
    }
    renormalizePhase(pMs, muD);

    // Diagnostics: pMax and g.
    double pMaxPhys = 0.0;
    for (double v : pMs) pMaxPhys = std::max(pMaxPhys, v);
    const double g = 2.0 * kPi * trapMu(N, dMu, [&](int i) {
      return pMs[(std::size_t)i] * muD[(std::size_t)i];
    });

    out.pMax[(std::size_t)c] = (float)pMaxPhys;
    out.asymmetryG[(std::size_t)c] = (float)(std::isfinite(g) ? g : 0.0);

    if (c == 0) setChannel(c, out.phaseR, pMs, pMaxPhys);
    else if (c == 1) setChannel(c, out.phaseG, pMs, pMaxPhys);
    else setChannel(c, out.phaseB, pMs, pMaxPhys);
  }

  return out;
}

} // namespace stellar::proc
