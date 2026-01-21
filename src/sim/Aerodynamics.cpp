#include "stellar/sim/Aerodynamics.h"

#include "stellar/math/Math.h"

#include <algorithm>
#include <cmath>

namespace stellar::sim {

static double smoothstep(double a, double b, double x) {
  if (!(b > a)) return (x >= b) ? 1.0 : 0.0;
  const double t = std::clamp((x - a) / (b - a), 0.0, 1.0);
  return t * t * (3.0 - 2.0 * t);
}

static math::Vec3d clampVecMag(const math::Vec3d& v, double maxMag) {
  if (!(maxMag > 0.0)) return v;
  const double m2 = v.lengthSq();
  if (m2 <= maxMag * maxMag) return v;
  const double m = std::sqrt(std::max(0.0, m2));
  if (!(m > 1e-12)) return {0, 0, 0};
  return v * (maxMag / m);
}

AerodynamicsSample computeAerodynamics(const AtmosphereSample& atmo,
                                      const math::Quatd& shipOrient,
                                      const math::Vec3d& shipAngVelLocalRadS,
                                      double shipMassKg,
                                      const AerodynamicsParams& params,
                                      math::Vec3d controlInputLocal) {
  AerodynamicsSample out{};
  out.active = false;

  if (!params.enabled) return out;
  if (!atmo.inAtmosphere) return out;
  if (!(shipMassKg > 1e-6)) return out;

  const double q = atmo.dynamicPressurePa;
  out.dynamicPressurePa = q;

  if (!(q >= params.minDynamicPressurePa)) return out;

  const math::Vec3d vRel = atmo.relVelKmS;
  const double vRelMag = vRel.length();
  if (!(vRelMag > 1e-9)) return out;
  const math::Vec3d vHat = vRel / vRelMag;

  // Compute body-frame wind direction (opposite the velocity).
  const math::Vec3d vLocal = shipOrient.conjugate().rotate(vRel);
  const double vLocalMag = vLocal.length();
  if (!(vLocalMag > 1e-9)) return out;
  const math::Vec3d windLocal = (-vLocal) / vLocalMag;

  // Convention:
  //   - Body axes: +X right, +Y up, +Z forward.
  //   - Relative wind points *toward* the ship.
  //   - alpha (AoA): rotation around +X (pitch plane), alpha>0 = nose up.
  //   - beta (sideslip): rotation around +Y (yaw plane), beta>0 = wind from left.
  const double alpha = std::atan2(windLocal.y, -windLocal.z);
  const double beta = std::atan2(windLocal.x, -windLocal.z);

  out.alphaRad = alpha;
  out.betaRad = beta;

  // --- Stall shaping (soft roll-off) ---
  const double stallStart = math::degToRad(std::max(0.1, params.stallAngleDeg));
  // A bit wider than onset so the response isn't razor sharp.
  const double stallEnd = std::max(stallStart + math::degToRad(8.0), stallStart * 1.8);
  const double aAbs = std::abs(alpha);
  const double stall01 = (aAbs <= stallStart) ? 0.0 : smoothstep(stallStart, stallEnd, aAbs);
  out.stall01 = stall01;

  // --- Lift coefficient ---
  double cl = params.liftSlopePerRad * alpha;
  const double clMax = std::max(0.0, params.clMax);
  if (clMax > 0.0) {
    cl = std::clamp(cl, -clMax, clMax);
  }

  // Roll-off in deep stall (keeps sign but reduces magnitude).
  cl *= (1.0 - 0.85 * stall01);

  out.cl = cl;

  const double sM2 = std::max(0.0, params.wingAreaM2);

  // --- Lift acceleration ---
  if (sM2 > 0.0 && std::abs(cl) > 1e-9) {
    // Lift direction starts from a ship-local axis (typically +Y) and is projected
    // onto the plane perpendicular to the flow.
    math::Vec3d liftAxisW = shipOrient.rotate(params.liftAxisLocal);
    if (liftAxisW.lengthSq() < 1e-18) {
      liftAxisW = shipOrient.rotate({0, 1, 0});
    }
    liftAxisW = liftAxisW.normalized();

    math::Vec3d liftDirW = liftAxisW - vHat * math::dot(liftAxisW, vHat);
    if (liftDirW.lengthSq() < 1e-12) {
      // Degenerate (axis aligned to flow). Fall back to a perpendicular vector.
      const math::Vec3d rightW = shipOrient.rotate({1, 0, 0}).normalized();
      liftDirW = math::cross(vHat, rightW);
    }

    if (liftDirW.lengthSq() > 1e-12) {
      liftDirW = liftDirW.normalized();

      // Lift: L = q * S * CL (Newtons).
      const double liftN = q * sM2 * cl;
      const double liftMS2 = liftN / shipMassKg;
      double liftKmS2 = liftMS2 / 1000.0;

      if (params.maxLiftAccelKmS2 > 1e-12) {
        liftKmS2 = std::clamp(liftKmS2, -params.maxLiftAccelKmS2, params.maxLiftAccelKmS2);
      }

      out.liftAccelKmS2 = liftDirW * liftKmS2;
    }
  }

  // --- Extra drag (induced + stall) ---
  double cdExtra = 0.0;
  if (sM2 > 0.0) {
    const double cdi = std::max(0.0, params.inducedDragK) * (cl * cl);
    const double cdStall = std::max(0.0, params.stallDragCd) * (stall01 * stall01);
    cdExtra = cdi + cdStall;
    out.cdExtra = cdExtra;

    if (cdExtra > 1e-9) {
      // Extra drag: D = q * S * CD_extra.
      const double dragN = q * sM2 * cdExtra;
      const double dragMS2 = dragN / shipMassKg;
      double dragKmS2 = dragMS2 / 1000.0;

      if (params.maxExtraDragAccelKmS2 > 1e-12) {
        dragKmS2 = std::min(dragKmS2, params.maxExtraDragAccelKmS2);
      }

      out.extraDragAccelKmS2 = (-vHat) * dragKmS2;
    }
  }

  // --- Angular acceleration (body-local) ---
  const double qRef = std::max(1.0, params.qRefPa);
  double qScale = q / qRef;
  qScale = std::clamp(qScale, 0.0, std::max(0.0, params.qScaleMax));

  math::Vec3d angAccel{0, 0, 0};

  // Aerodynamic rotational damping.
  angAccel.x += -shipAngVelLocalRadS.x * params.dampingPerS.x * qScale;
  angAccel.y += -shipAngVelLocalRadS.y * params.dampingPerS.y * qScale;
  angAccel.z += -shipAngVelLocalRadS.z * params.dampingPerS.z * qScale;

  // Weathervane alignment (reduces AoA + sideslip over time).
  if (params.alignToVelocity) {
    angAccel.x += (-alpha) * params.alignGain * qScale;
    angAccel.y += (-beta) * params.alignGain * qScale;
  }

  // Control surfaces: treat pilot torque input as aerodynamic control.
  if (params.controlSurfaces) {
    controlInputLocal.x = std::clamp(controlInputLocal.x, -1.0, 1.0);
    controlInputLocal.y = std::clamp(controlInputLocal.y, -1.0, 1.0);
    controlInputLocal.z = std::clamp(controlInputLocal.z, -1.0, 1.0);

    angAccel.x += controlInputLocal.x * params.controlAngAccelRadS2.x * qScale;
    angAccel.y += controlInputLocal.y * params.controlAngAccelRadS2.y * qScale;
    angAccel.z += controlInputLocal.z * params.controlAngAccelRadS2.z * qScale;
  }

  if (params.maxAngAccelRadS2 > 1e-12) {
    angAccel = clampVecMag(angAccel, params.maxAngAccelRadS2);
  }

  out.angAccelLocalRadS2 = angAccel;

  out.active = true;
  return out;
}

} // namespace stellar::sim
