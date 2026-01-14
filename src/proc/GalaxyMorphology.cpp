#include "stellar/proc/GalaxyMorphology.h"

#include "stellar/core/Hash.h"
#include "stellar/math/Math.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {

static double clamp01(double v) {
  if (v < 0.0) return 0.0;
  if (v > 1.0) return 1.0;
  return v;
}

double galaxyBarMask01(const GalaxyParams& gp, double xLy, double yLy) {
  // Treat strength==0 as disabled for convenience.
  if (gp.barStrength <= 0.0) return 0.0;

  const double a = std::max(1.0e-6, std::abs(gp.barLengthLy));
  const double b = std::max(1.0e-6, std::abs(gp.barWidthLy));
  const double ang = stellar::math::degToRad(gp.barAngleDeg);

  const double ca = std::cos(ang);
  const double sa = std::sin(ang);

  // Rotate (x,y) into bar local frame.
  const double xr = ca * xLy + sa * yLy;
  const double yr = -sa * xLy + ca * yLy;

  const double qx = xr / a;
  const double qy = yr / b;
  const double q = std::sqrt(qx * qx + qy * qy);
  if (q >= 1.0) return 0.0;

  const double power = std::max(1.0, gp.barPower);
  const double t = 1.0 - q;
  return clamp01(std::pow(t, power));
}

double galaxyRingMask01(const GalaxyParams& gp, double rxyLy) {
  if (gp.ringStrength <= 0.0) return 0.0;

  const double sigma = std::max(1.0e-6, std::abs(gp.ringWidthLy));
  const double d = (rxyLy - gp.ringRadiusLy) / sigma;
  double m = std::exp(-0.5 * d * d);

  const double p = std::max(1.0, gp.ringPower);
  if (p != 1.0) m = std::pow(m, p);

  return clamp01(m);
}

double galaxyThicknessHalfLy(const GalaxyParams& gp, double rxyLy) {
  const double base = std::max(1.0, gp.thicknessLy * 0.5);
  const double flareS = std::max(0.0, gp.flareStrength);
  if (flareS <= 0.0) return base;

  const double R = std::max(1.0, gp.radiusLy);
  const double t0 = std::clamp(std::max(0.0, rxyLy) / R, 0.0, 1.0);
  const double p = std::max(0.05, gp.flarePower);
  const double t = std::pow(t0, p);

  return base * (1.0 + flareS * t);
}

double galaxyWarpZLy(core::u64 seed, const GalaxyParams& gp, double xLy, double yLy) {
  const double amp0 = gp.warpAmplitudeLy;
  if (std::abs(amp0) <= 1.0e-9) return 0.0;

  const double r = std::sqrt(xLy * xLy + yLy * yLy);
  const double start = std::max(0.0, gp.warpStartRadiusLy);

  if (r <= start) return 0.0;

  const double R = std::max(start + 1.0, gp.radiusLy);
  const double t0 = std::clamp((r - start) / std::max(1.0, (R - start)), 0.0, 1.0);
  const double p = std::max(0.05, gp.warpPower);
  const double amp = amp0 * std::pow(t0, p);

  int lobes = gp.warpLobes;
  if (lobes == 0) lobes = 2;
  lobes = std::clamp(std::abs(lobes), 1, 8);

  const double theta = std::atan2(yLy, xLy);
  const double phase = stellar::math::degToRad(gp.warpPhaseDeg);

  double w = amp * std::sin(static_cast<double>(lobes) * theta + phase);

  // Optional low-frequency modulation so the warp isn't perfectly sinusoidal.
  const double ns = std::clamp(gp.warpNoiseStrength, 0.0, 1.0);
  const double nf = gp.warpNoiseFreq;
  if (ns > 0.0 && nf > 0.0) {
    const core::u64 nSeed = core::hashCombine(seed, core::fnv1a64("galaxy_warp_noise"));
    const double n = perlin2D(nSeed, xLy * nf, yLy * nf); // 0..1
    const double mul = (1.0 - ns) + (2.0 * ns) * n;       // [1-ns, 1+ns]
    w *= mul;
  }

  return w;
}

GalaxyMorphologySample sampleGalaxyMorphology(core::u64 seed, const GalaxyParams& gp, const math::Vec3d& posLy) {
  GalaxyMorphologySample out{};

  const double rxy = std::sqrt(posLy.x * posLy.x + posLy.y * posLy.y);

  out.thicknessHalfLy = galaxyThicknessHalfLy(gp, rxy);
  out.warpZLy = galaxyWarpZLy(seed, gp, posLy.x, posLy.y);

  out.bar01 = galaxyBarMask01(gp, posLy.x, posLy.y);
  out.ring01 = galaxyRingMask01(gp, rxy);

  const double bs = std::max(0.0, gp.barStrength);
  const double rs = std::max(0.0, gp.ringStrength);

  out.densityMul = (1.0 + bs * out.bar01) * (1.0 + rs * out.ring01);
  if (out.densityMul < 0.0) out.densityMul = 0.0;

  return out;
}

double galaxyMorphologyDensityMul(core::u64 seed, const GalaxyParams& gp, const math::Vec3d& posLy) {
  const auto s = sampleGalaxyMorphology(seed, gp, posLy);
  return s.densityMul;
}

} // namespace stellar::proc
