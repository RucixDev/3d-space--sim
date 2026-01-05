#include "stellar/sim/Gravity.h"

#include <cmath>

namespace stellar::sim {

double muStarKm3S2(const Star& star) {
  return muFromSolarMass(star.massSol);
}

double muPlanetKm3S2(const Planet& planet) {
  return muFromEarthMass(planet.massEarth);
}

double radiusStarKm(const Star& star) {
  return radiusKmFromSolarRadius(star.radiusSol);
}

double radiusPlanetKm(const Planet& planet) {
  return radiusKmFromEarthRadius(planet.radiusEarth);
}

void gatherGravityBodies(const StarSystem& sys,
                         double timeDays,
                         std::vector<GravityBody>& out,
                         const GravityParams& params) {
  out.clear();

  if (params.includeStar) {
    GravityBody b{};
    b.kind = GravityBody::Kind::Star;
    b.id = sys.stub.id;
    b.name = "Star";
    b.posKm = {0, 0, 0};
    b.velKmS = {0, 0, 0};
    b.muKm3S2 = muStarKm3S2(sys.star);
    b.radiusKm = radiusStarKm(sys.star);
    out.push_back(std::move(b));
  }

  if (params.includePlanets) {
    out.reserve(out.size() + sys.planets.size());
    for (std::size_t i = 0; i < sys.planets.size(); ++i) {
      const auto& p = sys.planets[i];
      GravityBody b{};
      b.kind = GravityBody::Kind::Planet;
      b.id = (core::u64)i;
      b.name = p.name;
      b.posKm = planetPosKm(p, timeDays);
      b.velKmS = planetVelKmS(p, timeDays);
      b.muKm3S2 = muPlanetKm3S2(p);
      b.radiusKm = radiusPlanetKm(p);
      out.push_back(std::move(b));
    }
  }
}

math::Vec3d gravityAccelFromBodyKmS2(const math::Vec3d& bodyPosKm,
                                    double muKm3S2,
                                    const math::Vec3d& queryPosKm,
                                    double minRadiusKm) {
  // r = (query - body). Accel points toward the body.
  const math::Vec3d r = queryPosKm - bodyPosKm;
  const double r2 = r.lengthSq();
  if (r2 <= 1e-18 || muKm3S2 <= 0.0) return {0, 0, 0};

  const double dist = std::sqrt(r2);
  const double denom = std::max(dist, std::max(1e-6, minRadiusKm));
  const double invR3 = 1.0 / (denom * denom * denom);
  return r * (-muKm3S2 * invR3);
}

math::Vec3d systemGravityAccelKmS2(const StarSystem& sys,
                                  double timeDays,
                                  const math::Vec3d& posKm,
                                  const GravityParams& params) {
  math::Vec3d acc{0, 0, 0};

  const double softenScale = std::max(1.0, params.softeningRadiusScale);

  if (params.includeStar) {
    const double mu = muStarKm3S2(sys.star);
    const double rMin = radiusStarKm(sys.star) * softenScale;
    acc += gravityAccelFromBodyKmS2({0, 0, 0}, mu, posKm, rMin);
  }

  if (params.includePlanets) {
    for (const auto& p : sys.planets) {
      const double mu = muPlanetKm3S2(p);
      const double rMin = radiusPlanetKm(p) * softenScale;
      acc += gravityAccelFromBodyKmS2(planetPosKm(p, timeDays), mu, posKm, rMin);
    }
  }

  acc *= params.scale;

  if (params.maxAccelKmS2 > 1e-12) {
    const double m = acc.length();
    if (m > params.maxAccelKmS2) {
      acc *= (params.maxAccelKmS2 / m);
    }
  }

  return acc;
}

DominantGravityBody dominantGravityBody(const StarSystem& sys,
                                       double timeDays,
                                       const math::Vec3d& posKm,
                                       const GravityParams& params) {
  DominantGravityBody out{};

  const double softenScale = std::max(1.0, params.softeningRadiusScale);

  auto consider = [&](GravityBody::Kind kind,
                      core::u64 id,
                      const std::string& name,
                      const math::Vec3d& bPos,
                      const math::Vec3d& bVel,
                      double mu,
                      double radiusKm) {
    const double rMin = radiusKm * softenScale;
    const math::Vec3d a = gravityAccelFromBodyKmS2(bPos, mu, posKm, rMin) * params.scale;
    const double m = a.length();
    if (!out.valid || m > out.accelMagKmS2) {
      out.valid = true;
      out.body.kind = kind;
      out.body.id = id;
      out.body.name = name;
      out.body.posKm = bPos;
      out.body.velKmS = bVel;
      out.body.muKm3S2 = mu;
      out.body.radiusKm = radiusKm;
      out.accelKmS2 = a;
      out.accelMagKmS2 = m;
    }
  };

  if (params.includeStar) {
    consider(GravityBody::Kind::Star,
             sys.stub.id,
             "Star",
             {0, 0, 0},
             {0, 0, 0},
             muStarKm3S2(sys.star),
             radiusStarKm(sys.star));
  }

  if (params.includePlanets) {
    for (std::size_t i = 0; i < sys.planets.size(); ++i) {
      const auto& p = sys.planets[i];
      consider(GravityBody::Kind::Planet,
               (core::u64)i,
               p.name,
               planetPosKm(p, timeDays),
               planetVelKmS(p, timeDays),
               muPlanetKm3S2(p),
               radiusPlanetKm(p));
    }
  }

  return out;
}

} // namespace stellar::sim
