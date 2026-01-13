#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace stellar::proc {

// How particle radii are distributed when generating the phase function.
//
// Mono:      single radius (fast)
// LogNormal: approximate a log-normal distribution by sampling multiple radii
//            in log-space and averaging the resulting phase functions.
enum class MieRadiusDistribution : int {
  Mono = 0,
  LogNormal = 1,
};

struct SpectralMiePhaseSettings {
  // Number of samples across μ = cos(θ) in [-1, 1]. This is also the LUT width.
  int angleSamples{256};

  // Wavelengths (nanometers) used for the spectral RGB channels.
  // Defaults are roughly: R=680nm, G=550nm, B=440nm.
  std::array<double, 3> wavelengthsNm{{680.0, 550.0, 440.0}};

  // Mean particle radius (micrometers).
  double radiusUm{0.20};

  // Optional particle size distribution.
  MieRadiusDistribution radiusDist{MieRadiusDistribution::Mono};

  // Log-normal sigma (dimensionless). Only used if radiusDist==LogNormal.
  double lognormalSigma{0.35};

  // Number of radii used to approximate the distribution. Only used if LogNormal.
  int radiusSampleCount{8};

  // Complex refractive index m = n + i*k (relative to the ambient medium).
  double refractiveIndexN{1.33};
  double refractiveIndexK{0.0};

  // If true, the returned phase curves are peak-normalized so max(p)=1 for each
  // channel (nice for driving real-time forward-scatter highlights).
  //
  // If false, the returned curves are physically normalized so the integral over
  // solid angle equals 1 (∫_{4π} p dΩ = 1).
  bool peakNormalize{true};
};

struct SpectralMiePhaseResult {
  int angleSamples{0};

  // μ values in [-1, 1] (size = angleSamples).
  std::vector<float> mu;

  // Phase curve samples per channel (size = angleSamples each).
  std::vector<float> phaseR;
  std::vector<float> phaseG;
  std::vector<float> phaseB;

  // Diagnostics per channel.
  std::array<float, 3> asymmetryG{{0.0f, 0.0f, 0.0f}}; // <cos θ>
  std::array<float, 3> qSca{{0.0f, 0.0f, 0.0f}};       // scattering efficiency (dimensionless)
  std::array<float, 3> pMax{{0.0f, 0.0f, 0.0f}};       // max of physically normalized p(μ)
};

// Generate a spectral Mie phase function (as a 1D LUT over μ = cos(θ)).
SpectralMiePhaseResult generateSpectralMiePhase(const SpectralMiePhaseSettings& s);

// Converts phase samples into an RGBA8 strip (width = angleSamples, height = 1).
// R,G,B map to phaseR/G/B. Alpha is set to 255.
std::vector<std::uint8_t> spectralPhaseToRgba8(const SpectralMiePhaseResult& r);

} // namespace stellar::proc
