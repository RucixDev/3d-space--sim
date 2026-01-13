#include "stellar/proc/MiePhase.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>

namespace stellar::proc {

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;

int mieNmax(double x) {
  if (x <= 1e-12) return 1;
  // Wiscombe-style truncation used in many practical Mie implementations.
  // See e.g. C. Mätzler (2002).
  const double nmaxD = std::round(2.0 + x + 4.0 * std::cbrt(x));
  int nmax = (int)std::max(1.0, nmaxD);
  // Guard rails for real-time tool usage.
  nmax = std::clamp(nmax, 1, 8192);
  return nmax;
}

void sphericalBesselJY(std::complex<double> z,
                      int nmax,
                      std::vector<std::complex<double>>& j,
                      std::vector<std::complex<double>>& y) {
  j.assign((std::size_t)nmax + 1, std::complex<double>(0.0, 0.0));
  y.assign((std::size_t)nmax + 1, std::complex<double>(0.0, 0.0));

  const double az = std::abs(z);
  if (az < 1e-12) {
    // z ~ 0: j0 ~ 1, j1 ~ z/3; y is singular (won't be used in this regime).
    j[0] = std::complex<double>(1.0, 0.0);
    if (nmax >= 1) j[1] = z / 3.0;
    for (int n = 2; n <= nmax; ++n) j[n] = std::complex<double>(0.0, 0.0);
    y[0] = std::complex<double>(-1e30, 0.0);
    for (int n = 1; n <= nmax; ++n) y[n] = std::complex<double>(-1e30, 0.0);
    return;
  }

  const std::complex<double> z2 = z * z;
  j[0] = std::sin(z) / z;
  if (nmax >= 1) {
    j[1] = std::sin(z) / z2 - std::cos(z) / z;
  }

  y[0] = -std::cos(z) / z;
  if (nmax >= 1) {
    y[1] = -std::cos(z) / z2 - std::sin(z) / z;
  }

  for (int n = 1; n < nmax; ++n) {
    const std::complex<double> coef = (2.0 * (double)n + 1.0) / z;
    j[n + 1] = coef * j[n] - j[n - 1];
    y[n + 1] = coef * y[n] - y[n - 1];
  }
}

struct MieCoeffs {
  int nmax{0};
  std::vector<std::complex<double>> a;
  std::vector<std::complex<double>> b;
  double qSca{0.0};
};

MieCoeffs mieCoefficients(std::complex<double> m, double x) {
  MieCoeffs out;
  if (x <= 1e-10) {
    out.nmax = 0;
    return out;
  }

  out.nmax = mieNmax(x);
  out.a.assign((std::size_t)out.nmax + 1, std::complex<double>(0.0, 0.0));
  out.b.assign((std::size_t)out.nmax + 1, std::complex<double>(0.0, 0.0));

  const std::complex<double> i(0.0, 1.0);
  const std::complex<double> cx(x, 0.0);
  const std::complex<double> z = m * cx;
  const std::complex<double> m2 = m * m;

  // Compute spherical Bessel functions at x (real) and z=m*x (complex).
  std::vector<std::complex<double>> jx, yx, jz, yz;
  sphericalBesselJY(cx, out.nmax, jx, yx);
  sphericalBesselJY(z, out.nmax, jz, yz);

  // Spherical Hankel (first kind) at x.
  std::vector<std::complex<double>> hx((std::size_t)out.nmax + 1);
  for (int n = 0; n <= out.nmax; ++n) {
    hx[n] = jx[n] + i * yx[n];
  }

  // Derivative helpers:
  //  d/dx [x * j_n(x)] = x * j_{n-1}(x) - n * j_n(x)
  //  d/dx [x * h_n(x)] = x * h_{n-1}(x) - n * h_n(x)
  std::vector<std::complex<double>> ax((std::size_t)out.nmax + 1);
  std::vector<std::complex<double>> ahx((std::size_t)out.nmax + 1);
  std::vector<std::complex<double>> az((std::size_t)out.nmax + 1);

  ax[0] = std::complex<double>(0.0, 0.0);
  ahx[0] = std::complex<double>(0.0, 0.0);
  az[0] = std::complex<double>(0.0, 0.0);

  for (int n = 1; n <= out.nmax; ++n) {
    ax[n] = cx * jx[n - 1] - (double)n * jx[n];
    ahx[n] = cx * hx[n - 1] - (double)n * hx[n];
    az[n] = z * jz[n - 1] - (double)n * jz[n];
  }

  // Mie coefficients a_n and b_n (Bohren & Huffman formalism).
  for (int n = 1; n <= out.nmax; ++n) {
    const std::complex<double> numA = m2 * jz[n] * ax[n] - jx[n] * az[n];
    const std::complex<double> denA = m2 * jz[n] * ahx[n] - hx[n] * az[n];

    const std::complex<double> numB = jz[n] * ax[n] - jx[n] * az[n];
    const std::complex<double> denB = jz[n] * ahx[n] - hx[n] * az[n];

    out.a[n] = numA / denA;
    out.b[n] = numB / denB;
  }

  // Scattering efficiency Q_sca (dimensionless), useful as a spectral strength hint.
  // Q_sca = (2/x^2) * Σ_{n=1..nmax} (2n+1)(|a_n|^2 + |b_n|^2)
  double sum = 0.0;
  for (int n = 1; n <= out.nmax; ++n) {
    sum += (2.0 * (double)n + 1.0) * (std::norm(out.a[n]) + std::norm(out.b[n]));
  }
  out.qSca = (2.0 * sum) / (x * x);

  return out;
}

void miePiTau(double u, int nmax, std::vector<double>& pi, std::vector<double>& tau) {
  pi.assign((std::size_t)nmax + 1, 0.0);
  tau.assign((std::size_t)nmax + 1, 0.0);

  if (nmax < 1) return;

  pi[1] = 1.0;
  tau[1] = u;

  if (nmax >= 2) {
    pi[2] = 3.0 * u;
    tau[2] = 3.0 * (2.0 * u * u - 1.0); // 3*cos(2*acos(u))
  }

  for (int n = 3; n <= nmax; ++n) {
    const double p1 = (2.0 * (double)n - 1.0) / ((double)n - 1.0) * pi[n - 1] * u;
    const double p2 = (double)n / ((double)n - 1.0) * pi[n - 2];
    pi[n] = p1 - p2;

    const double t1 = (double)n * u * pi[n];
    const double t2 = ((double)n + 1.0) * pi[n - 1];
    tau[n] = t1 - t2;
  }
}

double mieUnpolarizedIntensity(const MieCoeffs& c, double u) {
  if (c.nmax <= 0) {
    return 1.0; // treat as isotropic in the x->0 limit
  }

  std::vector<double> pi, tau;
  miePiTau(u, c.nmax, pi, tau);

  std::complex<double> S1(0.0, 0.0);
  std::complex<double> S2(0.0, 0.0);

  for (int n = 1; n <= c.nmax; ++n) {
    const double nn = (double)n;
    const double n2 = (2.0 * nn + 1.0) / (nn * (nn + 1.0));
    const double pin = n2 * pi[n];
    const double tin = n2 * tau[n];

    S1 += c.a[n] * pin + c.b[n] * tin;
    S2 += c.a[n] * tin + c.b[n] * pin;
  }

  // Unpolarized intensity (Mueller S11 up to a constant factor).
  return 0.5 * (std::norm(S1) + std::norm(S2));
}

struct RadiusSample {
  double rUm{0.0};
  double w{1.0};
};

std::vector<RadiusSample> buildRadiusSamples(const SpectralMiePhaseSettings& s) {
  std::vector<RadiusSample> out;
  if (s.radiusDist == MieRadiusDistribution::Mono || s.lognormalSigma <= 1e-6 || s.radiusSampleCount <= 1) {
    out.push_back(RadiusSample{s.radiusUm, 1.0});
    return out;
  }

  const int k = std::clamp(s.radiusSampleCount, 2, 64);
  const double sigma = std::max(1e-6, s.lognormalSigma);
  const double lnMean = std::log(std::max(1e-9, s.radiusUm));

  // Sample evenly in log-space over roughly ±3σ.
  const double lnMin = lnMean - 3.0 * sigma;
  const double lnMax = lnMean + 3.0 * sigma;

  out.reserve((std::size_t)k);
  double wSum = 0.0;
  for (int i = 0; i < k; ++i) {
    const double t = (k == 1) ? 0.0 : (double)i / (double)(k - 1);
    const double lnR = lnMin + (lnMax - lnMin) * t;
    const double r = std::exp(lnR);

    // Log-normal PDF in log-space is just a Gaussian; ΔlnR is constant, so we can
    // omit it and renormalize weights.
    const double z = (lnR - lnMean) / sigma;
    const double w = std::exp(-0.5 * z * z);

    out.push_back(RadiusSample{r, w});
    wSum += w;
  }

  if (wSum <= 0.0) {
    out.clear();
    out.push_back(RadiusSample{s.radiusUm, 1.0});
    return out;
  }

  for (auto& rs : out) rs.w /= wSum;
  return out;
}

void normalizePhase(const std::vector<double>& I,
                    const std::vector<double>& mu,
                    std::vector<double>& outP,
                    double& outPMax,
                    double& outG) {
  const int n = (int)I.size();
  outP.assign(I.begin(), I.end());
  outPMax = 0.0;
  outG = 0.0;

  if (n < 2) {
    if (n == 1) {
      outP[0] = 1.0 / (4.0 * kPi);
      outPMax = outP[0];
      outG = mu.empty() ? 0.0 : mu[0];
    }
    return;
  }

  const double dMu = (mu.back() - mu.front()) / (double)(n - 1);

  // Trapezoidal in μ.
  auto trap = [&](auto fn) {
    double sum = 0.5 * (fn(0) + fn(n - 1));
    for (int i = 1; i < n - 1; ++i) sum += fn(i);
    return sum * dMu;
  };

  const double intI = 2.0 * kPi * trap([&](int i) { return I[i]; });
  if (!(intI > 0.0) || !std::isfinite(intI)) {
    // Fallback: isotropic.
    const double pIso = 1.0 / (4.0 * kPi);
    for (int i = 0; i < n; ++i) outP[i] = pIso;
    outPMax = pIso;
    outG = 0.0;
    return;
  }

  // Physically normalized phase.
  for (int i = 0; i < n; ++i) {
    outP[i] = I[i] / intI;
    outPMax = std::max(outPMax, outP[i]);
  }

  // Asymmetry parameter g = <cosθ> = ∫ p(μ) μ dΩ.
  outG = 2.0 * kPi * trap([&](int i) { return outP[i] * mu[i]; });
  if (!std::isfinite(outG)) outG = 0.0;
}

} // namespace

SpectralMiePhaseResult generateSpectralMiePhase(const SpectralMiePhaseSettings& s) {
  SpectralMiePhaseResult out;

  const int N = std::clamp(s.angleSamples, 8, 4096);
  out.angleSamples = N;
  out.mu.resize((std::size_t)N);
  for (int i = 0; i < N; ++i) {
    const double t = (N == 1) ? 0.0 : (double)i / (double)(N - 1);
    out.mu[(std::size_t)i] = (float)(-1.0 + 2.0 * t);
  }

  // Convert to double mu array for integration.
  std::vector<double> muD((std::size_t)N);
  for (int i = 0; i < N; ++i) muD[(std::size_t)i] = (double)out.mu[(std::size_t)i];

  const std::complex<double> m(s.refractiveIndexN, s.refractiveIndexK);

  const auto radii = buildRadiusSamples(s);

  // Accumulate intensities per channel.
  std::array<std::vector<double>, 3> I{};
  for (int c = 0; c < 3; ++c) I[c].assign((std::size_t)N, 0.0);

  std::array<double, 3> qScaAcc{{0.0, 0.0, 0.0}};

  for (const auto& rs : radii) {
    const double rM = rs.rUm * 1e-6; // µm -> m

    for (int c = 0; c < 3; ++c) {
      const double lambdaM = std::max(1e-12, s.wavelengthsNm[(std::size_t)c]) * 1e-9;
      const double x = (2.0 * kPi * rM) / lambdaM;

      const auto coeff = mieCoefficients(m, x);
      qScaAcc[(std::size_t)c] += rs.w * coeff.qSca;

      for (int i = 0; i < N; ++i) {
        const double u = muD[(std::size_t)i];
        I[c][(std::size_t)i] += rs.w * mieUnpolarizedIntensity(coeff, u);
      }
    }
  }

  out.phaseR.resize((std::size_t)N);
  out.phaseG.resize((std::size_t)N);
  out.phaseB.resize((std::size_t)N);

  for (int c = 0; c < 3; ++c) {
    std::vector<double> Pnorm;
    double pMax = 0.0;
    double g = 0.0;
    normalizePhase(I[c], muD, Pnorm, pMax, g);

    out.pMax[(std::size_t)c] = (float)pMax;
    out.asymmetryG[(std::size_t)c] = (float)g;
    out.qSca[(std::size_t)c] = (float)qScaAcc[(std::size_t)c];

    const double invMax = (pMax > 0.0) ? (1.0 / pMax) : 1.0;

    for (int i = 0; i < N; ++i) {
      double v = Pnorm[(std::size_t)i];
      if (s.peakNormalize) v *= invMax;
      const float vf = (float)std::clamp(v, 0.0, 1e9);

      if (c == 0) out.phaseR[(std::size_t)i] = vf;
      else if (c == 1) out.phaseG[(std::size_t)i] = vf;
      else out.phaseB[(std::size_t)i] = vf;
    }
  }

  return out;
}

std::vector<std::uint8_t> spectralPhaseToRgba8(const SpectralMiePhaseResult& r) {
  const int N = r.angleSamples;
  std::vector<std::uint8_t> out;
  if (N <= 0) return out;

  out.resize((std::size_t)N * 4);
  for (int i = 0; i < N; ++i) {
    const float pr = (i < (int)r.phaseR.size()) ? r.phaseR[(std::size_t)i] : 0.0f;
    const float pg = (i < (int)r.phaseG.size()) ? r.phaseG[(std::size_t)i] : 0.0f;
    const float pb = (i < (int)r.phaseB.size()) ? r.phaseB[(std::size_t)i] : 0.0f;

    auto toU8 = [](float v) -> std::uint8_t {
      const float c = std::clamp(v, 0.0f, 1.0f);
      return (std::uint8_t)std::lround(c * 255.0f);
    };

    out[(std::size_t)i * 4 + 0] = toU8(pr);
    out[(std::size_t)i * 4 + 1] = toU8(pg);
    out[(std::size_t)i * 4 + 2] = toU8(pb);
    out[(std::size_t)i * 4 + 3] = 255;
  }

  return out;
}

} // namespace stellar::proc
