#pragma once

#include "stellar/proc/MiePhase.h"

namespace stellar::proc {

// Analytic multiple-scattering phase function approximation.
//
// Our spectral Mie LUT generator produces a *single-scattering* phase function
// p(\mu) (\mu = cos\theta) normalized such that \int_{4\pi} p d\Omega = 1.
//
// For a rotationally-symmetric phase function, a convenient basis is the
// Legendre expansion (aka zonal harmonics):
//   p(\mu) = (1 / 4\pi) \sum_{l=0..\infty} (2l+1) \chi_l P_l(\mu)
// where \chi_l are the (azimuthally integrated) Legendre moments.
//
// Repeated scattering corresponds to an n-fold spherical convolution of the
// phase function with itself. In the zonal/Legendre domain, convolution becomes
// multiplication, so the n-th scattering order has moments \chi_l^n.
//
// If we assume a geometric distribution over scattering order controlled by an
// effective single-scattering albedo \omega (0..1), we can sum orders in closed
// form and reconstruct an *effective* multiple-scattering phase function.
// This is a cheap, stable way to get a broader (more isotropic) lobe that is
// consistent with the single-scattering phase.

enum class MultipleScatteringMode : int {
  // Normalized mixture of orders n >= 1 (single + multiple).
  TotalOrders = 0,
  // Normalized mixture of orders n >= 2 (multiple only).
  MultipleOnly = 1,
};

struct MultipleScatteringPhaseSettings {
  // Maximum Legendre order L used for reconstruction (0..L).
  int legendreOrder{32};

  // Effective albedo \omega in [0, 0.999]. Must be < 1 for convergence.
  double scatteringAlbedo{0.92};

  MultipleScatteringMode mode{MultipleScatteringMode::MultipleOnly};

  // If true, peak-normalize each output channel so max(p)=1 (useful for LUTs).
  // If false, keep the physically normalized phase function.
  bool peakNormalize{true};

  // Set to true if the input phase samples were peak-normalized.
  // We can recover the physically normalized samples using the input pMax.
  bool inputWasPeakNormalized{true};
};

// Builds an effective multiple-scattering phase function from an existing
// spectral phase LUT.
//
// The returned SpectralMiePhaseResult uses the same \mu sampling as the input.
// qSca is copied through for convenience (it does not change under this model).
SpectralMiePhaseResult generateMultipleScatteringPhase(
    const SpectralMiePhaseResult& single,
    const MultipleScatteringPhaseSettings& s);

} // namespace stellar::proc
