#include "stellar/sim/Atmosphere.h"

#include "stellar/math/Math.h"
#include "stellar/sim/Units.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double safeExp(double x) {
  // Prevent extreme under/overflow while remaining deterministic.
  x = std::clamp(x, -80.0, 80.0);
  return std::exp(x);
}

AtmosphereModel atmosphereModelForPlanet(const Planet& p) {
  AtmosphereModel m{};

  // Base (Earth-ish) parameters per type at gRel ~= 1.
  // These are intentionally "gamey".
  bool has = true;
  double rho0 = 0.0;   // kg/m^3
  double H = 0.0;      // km
  double top = 0.0;    // km

  switch (p.type) {
    case PlanetType::Ocean:
      rho0 = 1.20;
      H = 8.5;
      top = 120.0;
      break;
    case PlanetType::Desert:
      rho0 = 0.90;
      H = 7.2;
      top = 90.0;
      break;
    case PlanetType::Ice:
      rho0 = 0.65;
      H = 6.5;
      top = 80.0;
      break;
    case PlanetType::GasGiant:
      // Gas giant "upper atmosphere" only. The near-surface region is not
      // represented; we just want a thick-ish drag well.
      rho0 = 3.50;
      H = 32.0;
      top = 900.0;
      break;
    case PlanetType::Rocky:
    default:
      // Treat most rocky worlds as airless. Some large rocky planets can have
      // a thin atmosphere.
      if (p.massEarth >= 0.35 && p.radiusEarth >= 0.55) {
        rho0 = 0.12;
        H = 6.0;
        top = 55.0;
      } else {
        has = false;
      }
      break;
  }

  if (!has) return m;

  const double r = std::max(0.05, p.radiusEarth);
  const double gRel = std::max(0.05, p.massEarth) / (r * r); // relative to Earth

  // Simple scaling: lower surface gravity => puffier atmosphere.
  const double gClamp = std::clamp(gRel, 0.25, 2.5);
  const double HScaled = H / gClamp;
  const double topScaled = top * (HScaled / std::max(1e-6, H));

  m.hasAtmosphere = true;
  m.seaLevelDensityKgM3 = rho0;
  m.scaleHeightKm = std::clamp(HScaled, 1.5, 120.0);
  m.topAltitudeKm = std::clamp(topScaled, 10.0, 4000.0);
  return m;
}

double atmosphereDensityKgM3(const AtmosphereModel& model, double altitudeKm) {
  if (!model.hasAtmosphere) return 0.0;
  if (!(altitudeKm >= 0.0)) altitudeKm = 0.0;
  if (altitudeKm >= model.topAltitudeKm) return 0.0;
  if (model.scaleHeightKm <= 1e-9) return 0.0;
  return model.seaLevelDensityKgM3 * safeExp(-altitudeKm / model.scaleHeightKm);
}

AtmosphereSample sampleSystemAtmosphere(const StarSystem& sys,
                                        double timeDays,
                                        const math::Vec3d& shipPosKm,
                                        const math::Vec3d& shipVelKmS,
                                        double shipMassKg,
                                        const AtmosphereParams& params) {
  AtmosphereSample best{};
  best.inAtmosphere = false;

  if (!params.includePlanets) return best;
  if (!(shipMassKg > 1e-6)) return best;
  if (sys.planets.empty()) return best;

  double bestQ = 0.0;

  for (std::size_t i = 0; i < sys.planets.size(); ++i) {
    const Planet& p = sys.planets[i];
    const AtmosphereModel model = atmosphereModelForPlanet(p);
    if (!model.hasAtmosphere) continue;

    const math::Vec3d pPos = planetPosKm(p, timeDays);
    const math::Vec3d pVel = planetVelKmS(p, timeDays);

    const double radiusKm = std::max(0.05, p.radiusEarth) * 6371.0;
    const math::Vec3d relPos = shipPosKm - pPos;
    const double distKm = relPos.length();
    const double altKm = distKm - radiusKm;
    if (altKm >= model.topAltitudeKm) continue;

    double rho = atmosphereDensityKgM3(model, altKm) * params.densityScale;
    if (!(rho > 1e-12)) continue;

    math::Vec3d vRel = shipVelKmS - pVel;
    double vRelMag = vRel.length();
    if (params.maxRelSpeedKmS > 0.0) {
      vRelMag = std::min(vRelMag, params.maxRelSpeedKmS);
    }
    if (!(vRelMag > 1e-9)) continue;

    // Dynamic pressure q = 0.5 * rho * v^2, with v in m/s.
    const double vMS = vRelMag * 1000.0;
    const double q = 0.5 * rho * (vMS * vMS);
    if (!(q > bestQ)) continue;

    bestQ = q;
    best.inAtmosphere = true;
    best.planetIndex = i;
    best.planetName = p.name;
    best.planetRadiusKm = radiusKm;
    best.altitudeKm = std::max(0.0, altKm);
    best.densityKgM3 = rho;
    best.dynamicPressurePa = q;
    best.relVelKmS = (vRel.lengthSq() > 1e-12) ? vRel : math::Vec3d{0, 0, 0};
    best.relSpeedKmS = vRelMag;

    // Drag acceleration magnitude.
    const double CdA = std::max(0.0, params.dragCd) * std::max(0.0, params.referenceAreaM2);
    double aMS2 = (shipMassKg > 0.0) ? (q * CdA / shipMassKg) : 0.0;
    double aKmS2 = aMS2 / 1000.0;

    if (params.maxDecelKmS2 > 0.0) {
      aKmS2 = std::min(aKmS2, params.maxDecelKmS2);
    }

    if (vRel.lengthSq() > 1e-12) {
      best.dragAccelKmS2 = (-vRel.normalized()) * aKmS2;
    } else {
      best.dragAccelKmS2 = {0, 0, 0};
    }

    if (params.applyHeating) {
      const double qKPa = q / 1000.0;
      best.heatingHeatPerSec = qKPa * params.heatPerKPaSec * params.heatingScale;
    } else {
      best.heatingHeatPerSec = 0.0;
    }
  }

  return best;
}

} // namespace stellar::sim
