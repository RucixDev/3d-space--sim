#include "OrbitAnalyzerWindow.h"

#include "stellar/math/Math.h"
#include "stellar/sim/OrbitalMechanics.h"
#include "stellar/sim/Units.h"

#include "imgui.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <string>

namespace stellar::game {

namespace {

struct RefBody {
  bool valid{false};
  stellar::sim::GravityBody body{};
};

static RefBody pickRefBody(const stellar::sim::StarSystem& sys,
                           double timeDays,
                           const stellar::sim::Ship& ship,
                           const stellar::sim::GravityParams& gParams,
                           int choice) {
  RefBody out{};

  // -1 = auto/dominant, 0 = star, 1..n = planet index+1
  if (choice == -1) {
    const auto dom = stellar::sim::dominantGravityBody(sys, timeDays, ship.positionKm(), gParams);
    if (dom.valid) {
      out.valid = true;
      out.body = dom.body;
    }
    return out;
  }

  if (choice == 0) {
    // Star at origin in local system frame.
    stellar::sim::GravityBody b;
    b.kind = stellar::sim::GravityBody::Kind::Star;
    b.id = sys.stub.id;
    b.name = sys.stub.name.empty() ? "Star" : (sys.stub.name + " (Star)");
    b.posKm = {0, 0, 0};
    b.velKmS = {0, 0, 0};
    b.muKm3S2 = stellar::sim::muStarKm3S2(sys.star);
    b.radiusKm = stellar::sim::radiusStarKm(sys.star);
    out.valid = true;
    out.body = b;
    return out;
  }

  const int planetIndex = choice - 1;
  if (planetIndex >= 0 && planetIndex < (int)sys.planets.size()) {
    const auto& p = sys.planets[(size_t)planetIndex];
    stellar::sim::GravityBody b;
    b.kind = stellar::sim::GravityBody::Kind::Planet;
    b.id = (stellar::core::u64)planetIndex;
    b.name = p.name;
    b.posKm = stellar::sim::planetPosKm(p, timeDays);
    b.velKmS = stellar::sim::planetVelKmS(p, timeDays);
    b.muKm3S2 = stellar::sim::muPlanetKm3S2(p);
    b.radiusKm = stellar::sim::radiusPlanetKm(p);
    out.valid = true;
    out.body = b;
    return out;
  }

  return out;
}

static void computeRtnBasis(const stellar::math::Vec3d& relPosKm,
                            const stellar::math::Vec3d& relVelKmS,
                            stellar::math::Vec3d& radial,
                            stellar::math::Vec3d& along,
                            stellar::math::Vec3d& normal) {
  radial = relPosKm;
  const double r2 = radial.lengthSq();
  if (r2 > 1e-18) {
    radial *= 1.0 / std::sqrt(r2);
  } else {
    radial = {1, 0, 0};
  }

  normal = stellar::math::cross(relPosKm, relVelKmS);
  const double n2 = normal.lengthSq();
  if (n2 > 1e-18) {
    normal *= 1.0 / std::sqrt(n2);
  } else {
    normal = {0, 0, 1};
  }

  along = stellar::math::cross(normal, radial);
  const double a2 = along.lengthSq();
  if (a2 > 1e-18) {
    along *= 1.0 / std::sqrt(a2);
  } else {
    along = {0, 1, 0};
  }
}


static stellar::math::Vec3d rotateAroundAxis(const stellar::math::Vec3d& v,
                                             const stellar::math::Vec3d& axisUnit,
                                             double angleRad) {
  // Rodrigues' rotation formula: v' = v c + (k×v) s + k (k·v) (1-c)
  const double c = std::cos(angleRad);
  const double s = std::sin(angleRad);
  const auto k = axisUnit;
  return (v * c) + (stellar::math::cross(k, v) * s) + (k * (stellar::math::dot(k, v) * (1.0 - c)));
}

static bool computePlaneAlignDvAtNode(const stellar::math::Vec3d& relPosKm,
                                      const stellar::math::Vec3d& relVelKmS,
                                      bool forcePrograde,
                                      stellar::math::Vec3d& outDvKmS) {
  outDvKmS = {0, 0, 0};

  const auto h = stellar::math::cross(relPosKm, relVelKmS);
  const double h2 = h.lengthSq();
  if (h2 < 1e-18) return false;

  const auto hHat = h * (1.0 / std::sqrt(h2));

  stellar::math::Vec3d targetNormal{0, 0, 1};
  if (!forcePrograde && hHat.z < 0.0) targetNormal = {0, 0, -1};

  const double dotN = stellar::math::clamp(stellar::math::dot(hHat, targetNormal), -1.0, 1.0);
  const double angle = std::acos(dotN);
  if (angle < 1e-9) {
    outDvKmS = {0, 0, 0};
    return true;
  }

  stellar::math::Vec3d rHat = relPosKm;
  const double r2 = rHat.lengthSq();
  if (r2 > 1e-18) {
    rHat *= 1.0 / std::sqrt(r2);
  } else {
    rHat = {1, 0, 0};
  }

  // Choose the sign that best aligns the rotated orbit normal to the target normal.
  const auto hPlus = rotateAroundAxis(hHat, rHat, angle);
  const auto hMinus = rotateAroundAxis(hHat, rHat, -angle);
  const double dp = stellar::math::dot(hPlus, targetNormal);
  const double dm = stellar::math::dot(hMinus, targetNormal);
  const double signedAngle = (dp >= dm) ? angle : -angle;

  const auto vDesired = rotateAroundAxis(relVelKmS, rHat, signedAngle);
  outDvKmS = vDesired - relVelKmS;
  return true;
}

static void writeNodeRtn(const OrbitAnalyzerBindings& b,
                         double nodeTimeSec,
                         const stellar::math::Vec3d& dvWorldKmS,
                         const stellar::math::Vec3d& radial,
                         const stellar::math::Vec3d& along,
                         const stellar::math::Vec3d& normal) {
  if (b.maneuverNodeEnabled) *b.maneuverNodeEnabled = true;
  if (b.maneuverNodeTimeSec) *b.maneuverNodeTimeSec = (float)std::max(0.0, nodeTimeSec);

  const double dvAlong = stellar::math::dot(dvWorldKmS, along);
  const double dvRad = stellar::math::dot(dvWorldKmS, radial);
  const double dvNorm = stellar::math::dot(dvWorldKmS, normal);

  if (b.dvAlongMS) *b.dvAlongMS = (float)(dvAlong * 1000.0);
  if (b.dvRadialMS) *b.dvRadialMS = (float)(dvRad * 1000.0);
  if (b.dvNormalMS) *b.dvNormalMS = (float)(dvNorm * 1000.0);
}

static std::string formatTime(double tSec) {
  if (!(tSec >= 0.0)) return "-";
  const int total = (int)std::round(tSec);
  const int s = total % 60;
  const int m = (total / 60) % 60;
  const int h = total / 3600;

  char buf[64];
  if (h > 0) {
    std::snprintf(buf, sizeof(buf), "%dh %dm %ds", h, m, s);
  } else if (m > 0) {
    std::snprintf(buf, sizeof(buf), "%dm %ds", m, s);
  } else {
    std::snprintf(buf, sizeof(buf), "%ds", s);
  }
  return std::string(buf);
}

static const char* orbitTypeLabel(stellar::sim::TwoBodyOrbit::Type t) {
  using T = stellar::sim::TwoBodyOrbit::Type;
  switch (t) {
    case T::Elliptic: return "Elliptic";
    case T::Hyperbolic: return "Hyperbolic";
    case T::Parabolic: return "Parabolic";
    default: return "Invalid";
  }
}

} // namespace

void drawOrbitAnalyzerWindow(OrbitAnalyzerWindowState& state,
                             const stellar::sim::StarSystem& sys,
                             double timeDays,
                             const stellar::sim::Ship& ship,
                             const stellar::sim::GravityParams& gravityParams,
                             OrbitAnalyzerBindings bindings) {
  if (!state.open) return;

  if (!ImGui::Begin("Orbit Analyzer", &state.open)) {
    ImGui::End();
    return;
  }

  // --- Reference body selection (optionally shared with trajectory preview) ---
  int choice = bindings.refBodyChoice ? *bindings.refBodyChoice : -1;
  choice = std::clamp(choice, -1, (int)sys.planets.size());

  if (ImGui::BeginCombo("Reference Body", (choice == -1) ? "Auto (dominant)" : (choice == 0) ? "Star" : sys.planets[(size_t)(choice - 1)].name.c_str())) {
    if (ImGui::Selectable("Auto (dominant)", choice == -1)) choice = -1;
    if (ImGui::Selectable("Star", choice == 0)) choice = 0;
    for (int i = 0; i < (int)sys.planets.size(); ++i) {
      const bool selected = (choice == i + 1);
      if (ImGui::Selectable(sys.planets[(size_t)i].name.c_str(), selected)) choice = i + 1;
    }
    ImGui::EndCombo();
  }

  if (bindings.refBodyChoice) *bindings.refBodyChoice = choice;

  ImGui::Checkbox("Use gravity scale", &state.useGravityScale);

  auto ref = pickRefBody(sys, timeDays, ship, gravityParams, choice);
  if (!ref.valid) {
    ImGui::TextDisabled("No reference body found (gravity bodies disabled?)");
    ImGui::End();
    return;
  }

  const double muBase = ref.body.muKm3S2;
  const double mu = state.useGravityScale ? (muBase * gravityParams.scale) : muBase;

  const stellar::math::Vec3d relPosKm = ship.positionKm() - ref.body.posKm;
  const stellar::math::Vec3d relVelKmS = ship.velocityKmS() - ref.body.velKmS;

  const auto orb = stellar::sim::solveTwoBodyOrbit(relPosKm, relVelKmS, mu);

  const double altKm = relPosKm.length() - ref.body.radiusKm;

  const auto el = stellar::sim::solveClassicalOrbitElements(relPosKm, relVelKmS, mu);

  ImGui::SeparatorText("Current orbit");
  ImGui::Text("Reference: %s", ref.body.name.empty() ? "(unnamed)" : ref.body.name.c_str());
  ImGui::Text("Type: %s", orbitTypeLabel(orb.type));
  ImGui::Text("Alt: %.0f km  |  r: %.0f km  |  |v|: %.3f km/s", altKm, orb.rKm, orb.vKmS);
  ImGui::Text("vr: %.3f km/s | vt: %.3f km/s", orb.radialSpeedKmS, orb.tangentialSpeedKmS);

  if (orb.type == stellar::sim::TwoBodyOrbit::Type::Elliptic) {
    const double periAlt = orb.periapsisKm - ref.body.radiusKm;
    const double apoAlt = orb.apoapsisKm - ref.body.radiusKm;

    ImGui::Text("a: %.0f km  e: %.4f", orb.semiMajorAxisKm, orb.eccentricity);
    ImGui::Text("Periapsis: %.0f km (alt %.0f km)", orb.periapsisKm, periAlt);
    ImGui::Text("Apoapsis:  %.0f km (alt %.0f km)", orb.apoapsisKm, apoAlt);
    ImGui::Text("Period: %s", formatTime(orb.periodSec).c_str());
    ImGui::Text("To periapsis: %s", formatTime(orb.timeToPeriapsisSec).c_str());
    ImGui::Text("To apoapsis:  %s", formatTime(orb.timeToApoapsisSec).c_str());
  } else if (orb.type == stellar::sim::TwoBodyOrbit::Type::Hyperbolic) {
    const double periAlt = orb.periapsisKm - ref.body.radiusKm;
    ImGui::Text("a: %.0f km  e: %.4f", orb.semiMajorAxisKm, orb.eccentricity);
    ImGui::Text("Periapsis: %.0f km (alt %.0f km)", orb.periapsisKm, periAlt);
    ImGui::Text("Time since periapsis: %s", formatTime(std::abs(orb.timeSincePeriapsisSec)).c_str());
  }


  if (el.valid) {
    ImGui::SeparatorText("Elements");
    ImGui::Text("Inclination: %.2f deg", stellar::math::radToDeg(el.inclinationRad));

    if (el.equatorial) {
      ImGui::TextDisabled("RAAN Ω: (equatorial)");
    } else {
      ImGui::Text("RAAN Ω: %.2f deg", stellar::math::radToDeg(el.raanRad));
    }

    if (el.circular) {
      ImGui::TextDisabled("Arg peri ω: (circular)");
    } else {
      ImGui::Text("Arg peri ω: %.2f deg", stellar::math::radToDeg(el.argPeriapsisRad));
    }

    ImGui::Text("True anomaly: %.2f deg", stellar::math::radToDeg(el.trueAnomalyRad));
  }

  ImGui::SeparatorText("Planner");

  const bool canWriteNode = bindings.maneuverNodeTimeSec && bindings.dvAlongMS && bindings.dvRadialMS && bindings.dvNormalMS;
  if (!canWriteNode) {
    ImGui::TextDisabled("Planner is not wired (missing maneuver node bindings)");
    ImGui::End();
    return;
  }

  // Targets
  ImGui::InputFloat("Target apo alt (km)", &state.targetApoAltKm, 100.0f, 1000.0f, "%.0f");
  ImGui::InputFloat("Target peri alt (km)", &state.targetPeriAltKm, 100.0f, 1000.0f, "%.0f");

  // Common basis vectors at the apses (two-body assumption)
  stellar::math::Vec3d nHat = orb.angularMomentumKm2S;
  if (nHat.lengthSq() > 1e-18) {
    nHat = nHat.normalized();
  } else {
    nHat = {0, 0, 1};
  }

  stellar::math::Vec3d periRad = relPosKm;
  if (orb.eccentricityVec.lengthSq() > 1e-18) {
    periRad = orb.eccentricityVec.normalized();
  } else if (periRad.lengthSq() > 1e-18) {
    periRad = periRad.normalized();
  } else {
    periRad = {1, 0, 0};
  }

  stellar::math::Vec3d apoRad = periRad * -1.0;

  // --- Circularize at apses ---
  if (orb.type == stellar::sim::TwoBodyOrbit::Type::Elliptic) {
    const double a = orb.semiMajorAxisKm;
    const double rp = orb.periapsisKm;
    const double ra = orb.apoapsisKm;

    const double vpCur = (rp > 0.0 && a != 0.0) ? std::sqrt(mu * (2.0 / rp - 1.0 / a)) : 0.0;
    const double vaCur = (ra > 0.0 && a != 0.0) ? std::sqrt(mu * (2.0 / ra - 1.0 / a)) : 0.0;

    const double vpCirc = (rp > 0.0) ? std::sqrt(mu / rp) : 0.0;
    const double vaCirc = (ra > 0.0) ? std::sqrt(mu / ra) : 0.0;

    stellar::math::Vec3d periAlong = stellar::math::cross(nHat, periRad).normalized();
    stellar::math::Vec3d apoAlong = stellar::math::cross(nHat, apoRad).normalized();

    if (ImGui::Button("Set node: circularize @ peri")) {
      const double dv = vpCirc - vpCur;
      writeNodeRtn(bindings, orb.timeToPeriapsisSec, periAlong * dv, periRad, periAlong, nHat);
    }
    ImGui::SameLine();
    if (ImGui::Button("Set node: circularize @ apo")) {
      const double dv = vaCirc - vaCur;
      writeNodeRtn(bindings, orb.timeToApoapsisSec, apoAlong * dv, apoRad, apoAlong, nHat);
    }

    ImGui::Spacing();

    // --- Adjust apoapsis (burn at periapsis) ---
    {
      const double raTarget = std::max(rp, ref.body.radiusKm + (double)state.targetApoAltKm);
      const double aNew = 0.5 * (rp + raTarget);
      const double vpNew = (rp > 0.0 && aNew != 0.0) ? std::sqrt(mu * (2.0 / rp - 1.0 / aNew)) : vpCur;
      const double dv = vpNew - vpCur;

      if (ImGui::Button("Set node: burn @ peri to target apo")) {
        writeNodeRtn(bindings, orb.timeToPeriapsisSec, periAlong * dv, periRad, periAlong, nHat);
      }
      ImGui::SameLine();
      ImGui::TextDisabled("(apo target %.0f km)", raTarget);
    }

    // --- Adjust periapsis (burn at apoapsis) ---
    {
      const double rpTarget = std::min(ra, ref.body.radiusKm + (double)state.targetPeriAltKm);
      const double aNew = 0.5 * (ra + rpTarget);
      const double vaNew = (ra > 0.0 && aNew != 0.0) ? std::sqrt(mu * (2.0 / ra - 1.0 / aNew)) : vaCur;
      const double dv = vaNew - vaCur;

      if (ImGui::Button("Set node: burn @ apo to target peri")) {
        writeNodeRtn(bindings, orb.timeToApoapsisSec, apoAlong * dv, apoRad, apoAlong, nHat);
      }
      ImGui::SameLine();
      ImGui::TextDisabled("(peri target %.0f km)", rpTarget);
    }

    ImGui::Spacing();

    // Escape at current radius (instantaneous)
    {
      const double vEsc = (orb.rKm > 0.0) ? std::sqrt(2.0 * mu / orb.rKm) : 0.0;
      const double dv = vEsc - orb.vKmS;
      stellar::math::Vec3d radial{}, along{}, normal{};
      computeRtnBasis(relPosKm, relVelKmS, radial, along, normal);

      if (ImGui::Button("Set node: escape now")) {
        writeNodeRtn(bindings, /*nodeTimeSec=*/0.0, along * dv, radial, along, normal);
      }
      ImGui::SameLine();
      ImGui::TextDisabled("(dv %.1f m/s)", dv * 1000.0);
    }

    ImGui::Spacing();
    ImGui::SeparatorText("Plane change (ref plane)");
    ImGui::Checkbox("Force prograde align", &state.planeAlignForcePrograde);

    if (el.valid && !el.equatorial) {
      auto wrap2pi = [](double a) {
        const double twoPi = 2.0 * stellar::math::kPi;
        a = std::fmod(a, twoPi);
        if (a < 0.0) a += twoPi;
        return a;
      };

      const double nuAsc = wrap2pi(-el.argPeriapsisRad);
      const double nuDesc = wrap2pi(stellar::math::kPi - el.argPeriapsisRad);

      const double dtAsc = stellar::sim::timeToTrueAnomalySecElliptic(el, nuAsc);
      const double dtDesc = stellar::sim::timeToTrueAnomalySecElliptic(el, nuDesc);

      stellar::math::Vec3d rAsc{}, vAsc{}, dvAsc{};
      stellar::math::Vec3d rDesc{}, vDesc{}, dvDesc{};

      const bool okAsc = (dtAsc >= 0.0) &&
                         stellar::sim::stateFromClassicalOrbitElementsAtTrueAnomaly(el, nuAsc, rAsc, vAsc) &&
                         computePlaneAlignDvAtNode(rAsc, vAsc, state.planeAlignForcePrograde, dvAsc);

      const bool okDesc = (dtDesc >= 0.0) &&
                          stellar::sim::stateFromClassicalOrbitElementsAtTrueAnomaly(el, nuDesc, rDesc, vDesc) &&
                          computePlaneAlignDvAtNode(rDesc, vDesc, state.planeAlignForcePrograde, dvDesc);

      if (okAsc) {
        const double altNodeKm = rAsc.length() - ref.body.radiusKm;
        ImGui::Text("Ascending node: in %s | alt %.0f km | dv %.1f m/s",
                    formatTime(dtAsc).c_str(),
                    altNodeKm,
                    dvAsc.length() * 1000.0);
      } else {
        ImGui::TextDisabled("Ascending node: (unavailable)");
      }

      ImGui::BeginDisabled(!okAsc);
      if (ImGui::Button("Set node: align plane @ AN")) {
        stellar::math::Vec3d radial{}, along{}, normal{};
        computeRtnBasis(rAsc, vAsc, radial, along, normal);
        writeNodeRtn(bindings, dtAsc, dvAsc, radial, along, normal);
      }
      ImGui::EndDisabled();

      if (okDesc) {
        const double altNodeKm = rDesc.length() - ref.body.radiusKm;
        ImGui::Text("Descending node: in %s | alt %.0f km | dv %.1f m/s",
                    formatTime(dtDesc).c_str(),
                    altNodeKm,
                    dvDesc.length() * 1000.0);
      } else {
        ImGui::TextDisabled("Descending node: (unavailable)");
      }

      ImGui::BeginDisabled(!okDesc);
      if (ImGui::Button("Set node: align plane @ DN")) {
        stellar::math::Vec3d radial{}, along{}, normal{};
        computeRtnBasis(rDesc, vDesc, radial, along, normal);
        writeNodeRtn(bindings, dtDesc, dvDesc, radial, along, normal);
      }
      ImGui::EndDisabled();

    } else {
      ImGui::TextDisabled("Plane align requires a non-equatorial bound orbit.");
    }

  } else {
    ImGui::TextDisabled("Planner actions require a bound (elliptic) orbit.");
  }

  ImGui::End();
}

} // namespace stellar::game
