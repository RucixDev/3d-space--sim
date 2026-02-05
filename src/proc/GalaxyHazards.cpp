#include "stellar/proc/GalaxyHazards.h"

#include "stellar/core/Hash.h"
#include "stellar/proc/GalaxyRegions.h"
#include "stellar/proc/Noise.h"

#include <algorithm>
#include <cmath>

namespace stellar::proc {

static double clamp01(double v) {
  return std::clamp(v, 0.0, 1.0);
}

static double smoothstep(double e0, double e1, double x) {
  if (e0 == e1) return (x < e0) ? 0.0 : 1.0;
  const double t = clamp01((x - e0) / (e1 - e0));
  return t * t * (3.0 - 2.0 * t);
}

static double ampSum(int octaves, double gain) {
  // Matches proc::fbmPerlin* amp schedule: amp0=0.5 then *= gain.
  octaves = std::clamp(octaves, 1, 12);
  gain = std::clamp(gain, 0.0, 0.999);
  double amp = 0.5;
  double sum = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp;
    amp *= gain;
  }
  return std::max(1e-12, sum);
}

static double unit01FromHash(core::u64 h) {
  // High 53 bits => [0,1)
  const core::u64 v = (h >> 11) & ((1ull << 53) - 1ull);
  return static_cast<double>(v) / static_cast<double>(1ull << 53);
}

static double signedUnitFromHash(core::u64 h) {
  return 2.0 * unit01FromHash(h) - 1.0;
}

const char* galaxyHazardKindName(GalaxyHazardKind k) {
  switch (k) {
    case GalaxyHazardKind::Clear: return "Clear";
    case GalaxyHazardKind::Nebula: return "Nebula";
    case GalaxyHazardKind::Storm: return "Storm";
    case GalaxyHazardKind::Rift: return "Rift";
    default: return "Unknown";
  }
}

GalaxyHazardSample sampleGalaxyHazards(core::u64 universeSeed,
                                       const math::Vec3d& posLy,
                                       const GalaxyHazardsParams& params) {
  GalaxyHazardSample out{};

  GalaxyHazardsParams p = params;
  p.nebulaScaleLy = std::max(1.0, p.nebulaScaleLy);
  p.stormScaleLy = std::max(1.0, p.stormScaleLy);
  p.nebulaOctaves = std::clamp(p.nebulaOctaves, 1, 10);
  p.stormOctaves = std::clamp(p.stormOctaves, 1, 10);
  p.driftLyPerDay = std::max(0.0, p.driftLyPerDay);
  p.regionInfluence = std::clamp(p.regionInfluence, 0.0, 1.5);
  p.strength = std::max(0.0, p.strength);

  // Drift direction derived from the seed.
  const core::u64 driftSalt = core::fnv1a64("galaxy_hazard_drift");
  const core::u64 baseH = core::hashCombine(universeSeed, driftSalt);
  math::Vec3d driftDir{
      signedUnitFromHash(core::hashCombine(baseH, 1ull)),
      signedUnitFromHash(core::hashCombine(baseH, 2ull)),
      signedUnitFromHash(core::hashCombine(baseH, 3ull)),
  };
  driftDir = driftDir.normalized();

  const math::Vec3d drift = driftDir * (p.driftLyPerDay * p.timeDays);
  const math::Vec3d qLy = posLy + drift;

  // Seeds.
  const core::u64 nebSeed = core::hashCombine(universeSeed, core::fnv1a64("galaxy_hazard_nebula"));
  const core::u64 stormSeed = core::hashCombine(universeSeed, core::fnv1a64("galaxy_hazard_storm"));

  // Nebula channel: ridged fractal => “clouds with edges”.
  {
    const math::Vec3d n = qLy / p.nebulaScaleLy;
    const double sumAmp = ampSum(p.nebulaOctaves, 0.5);
    double v = ridgedFbmPerlin3D(nebSeed, n.x, n.y, n.z, p.nebulaOctaves, 2.0, 0.5);
    v = clamp01(v / sumAmp);

    // Patchify: encourage broad quiet space with occasional dense nebula.
    v = smoothstep(0.42, 0.78, v);
    out.nebula01 = v;
  }

  // Storm channel: higher frequency, mostly inside nebulas.
  {
    const math::Vec3d s = qLy / p.stormScaleLy;
    const double sumAmp = ampSum(p.stormOctaves, 0.5);
    double v = fbmPerlin3D(stormSeed, s.x, s.y, s.z, p.stormOctaves, 2.0, 0.5);
    v = clamp01(v / sumAmp);

    // Make storms more "event-like": few strong cells.
    v = smoothstep(0.60, 0.88, v);

    // Storms rarely exist in clean space.
    v *= (0.12 + 0.88 * out.nebula01);
    out.storm01 = clamp01(v);
  }

  // Region modulation (art-direction): bias hazards to feel consistent with
  // GalaxyRegions without making hazards identical to regions.
  GalaxyRegionSample reg{};
  bool haveReg = false;
  if (p.regionCellSizeLy > 0.0 && p.regionInfluence > 0.0) {
    reg = sampleGalaxyRegion(universeSeed, posLy, p.regionCellSizeLy);
    haveReg = true;

    double nebBias = 0.0;
    double stormBias = 0.0;

    switch (reg.kind) {
      case GalaxyRegionKind::Core:
        nebBias = -0.10;
        stormBias = -0.16;
        break;
      case GalaxyRegionKind::InnerDisc:
        nebBias = 0.00;
        stormBias = 0.00;
        break;
      case GalaxyRegionKind::OuterRim:
        nebBias = 0.06;
        stormBias = 0.10;
        break;
      case GalaxyRegionKind::Nebula:
        nebBias = 0.28;
        stormBias = 0.16;
        break;
      case GalaxyRegionKind::Cluster:
        nebBias = 0.10;
        stormBias = 0.06;
        break;
      case GalaxyRegionKind::Rift:
        nebBias = 0.14;
        stormBias = 0.32;
        break;
      default:
        break;
    }

    // Boundaries are more chaotic.
    stormBias += 0.18 * clamp01(reg.edge01);

    out.nebula01 = clamp01(out.nebula01 + p.regionInfluence * nebBias);
    out.storm01 = clamp01(out.storm01 + p.regionInfluence * stormBias);
  }

  // Combined hazard: nebula is pervasive, storms are spiky.
  out.hazard01 = clamp01(p.strength * clamp01(0.70 * out.nebula01 + 0.95 * out.storm01));

  out.sensorOcclusion01 = out.nebula01;
  out.navDisruption01 = out.storm01;

  // Classification.
  out.kind = GalaxyHazardKind::Clear;
  if (out.hazard01 > 0.22) {
    if (out.storm01 > 0.60) {
      out.kind = GalaxyHazardKind::Storm;
    } else if (out.nebula01 > 0.50) {
      out.kind = GalaxyHazardKind::Nebula;
    }
  }

  // If we're in a Rift region, expose that explicitly when hazards are present.
  if (haveReg && reg.kind == GalaxyRegionKind::Rift && out.hazard01 > 0.35) {
    out.kind = GalaxyHazardKind::Rift;
  }

  return out;
}

double sampleGalaxyHazardAvgOnSegment(core::u64 universeSeed,
                                      const math::Vec3d& aLy,
                                      const math::Vec3d& bLy,
                                      const GalaxyHazardsParams& params,
                                      int samples) {
  samples = std::clamp(samples, 1, 11);
  const math::Vec3d d = bLy - aLy;
  double sum = 0.0;
  for (int i = 0; i < samples; ++i) {
    const double t = (static_cast<double>(i) + 0.5) / static_cast<double>(samples);
    const math::Vec3d p = aLy + d * t;
    sum += sampleGalaxyHazards(universeSeed, p, params).hazard01;
  }
  return sum / static_cast<double>(samples);
}



double sampleGalaxyNavDisruptionAvgOnSegment(core::u64 universeSeed,
                                             const math::Vec3d& aLy,
                                             const math::Vec3d& bLy,
                                             const GalaxyHazardsParams& params,
                                             int samples) {
  samples = std::clamp(samples, 1, 11);
  const math::Vec3d d = bLy - aLy;
  double sum = 0.0;
  for (int i = 0; i < samples; ++i) {
    const double t = (static_cast<double>(i) + 0.5) / static_cast<double>(samples);
    const math::Vec3d p = aLy + d * t;
    sum += sampleGalaxyHazards(universeSeed, p, params).navDisruption01;
  }
  return std::clamp(sum / static_cast<double>(samples), 0.0, 1.0);
}

double sampleGalaxySensorOcclusionAvgOnSegment(core::u64 universeSeed,
                                               const math::Vec3d& aLy,
                                               const math::Vec3d& bLy,
                                               const GalaxyHazardsParams& params,
                                               int samples) {
  samples = std::clamp(samples, 1, 11);
  const math::Vec3d d = bLy - aLy;
  double sum = 0.0;
  for (int i = 0; i < samples; ++i) {
    const double t = (static_cast<double>(i) + 0.5) / static_cast<double>(samples);
    const math::Vec3d p = aLy + d * t;
    sum += sampleGalaxyHazards(universeSeed, p, params).sensorOcclusion01;
  }
  return std::clamp(sum / static_cast<double>(samples), 0.0, 1.0);
}

} // namespace stellar::proc
