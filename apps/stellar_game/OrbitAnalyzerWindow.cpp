#include "OrbitAnalyzerWindow.h"

#include "stellar/math/Math.h"
#include "stellar/sim/ManeuverProgramComputer.h"
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

static void writePlanToNodeControls(const OrbitAnalyzerBindings& bindings,
                                   const stellar::sim::ManeuverProgramResult& res) {
  if (!res.valid) return;
  stellar::math::Vec3d radial{}, along{}, normal{};
  computeRtnBasis(res.nodeRelPosKm, res.nodeRelVelKmS, radial, along, normal);
  writeNodeRtn(bindings, res.timeToNodeSec, res.plan.deltaVWorldKmS, radial, along, normal);
}

static const char* gravityBodyKindLabel(stellar::sim::GravityBody::Kind kind) {
  switch (kind) {
  case stellar::sim::GravityBody::Kind::Star:
    return "Star";
  case stellar::sim::GravityBody::Kind::Planet:
    return "Planet";
  }
  return "?";
}

static std::string bodyLabel(const stellar::sim::GravityBody& body) {
  if (!body.name.empty()) {
    return body.name;
  }
  return gravityBodyKindLabel(body.kind);
}

static void recomputeForecast(OrbitAnalyzerWindowState& state,
                              const stellar::sim::StarSystem& sys,
                              double timeDays,
                              const stellar::sim::Ship& ship,
                              const stellar::sim::GravityParams& gravityParams) {
  const double horizonSec = std::max(0.0, static_cast<double>(state.forecastHorizonMin) * 60.0);
  const double stepSec = std::max(0.05, static_cast<double>(state.forecastStepSec));
  const int maxSamples = std::clamp(state.forecastMaxSamples, 2, 50000);

  stellar::sim::TrajectoryPredictParams pred{};
  pred.horizonSec = horizonSec;
  pred.stepSec = stepSec;
  pred.maxSamples = maxSamples;
  pred.includeGravity = true;
  pred.gravity = gravityParams;
  if (!state.useGravityScale) {
    pred.gravity.scale = 1.0;
  }

  state.forecastSamples =
      stellar::sim::predictTrajectoryRK4(sys, timeDays, ship.positionKm(), ship.velocityKmS(), pred);
  state.forecastLastComputeTimeDays = timeDays;

  if (state.forecastSamples.size() < 2) {
    state.forecastValid = false;
    state.forecastAnalysis = {};
    state.forecastDominantTransitions.clear();
    state.forecastStatus = "Forecast: insufficient samples";
    return;
  }

  state.forecastAnalysis = stellar::sim::analyzeTrajectory(sys, timeDays, state.forecastSamples, pred.gravity);

  stellar::sim::DominantBodyTransitionParams domParams{};
  domParams.refineDepth = std::clamp(state.forecastRefineDepth, 0, 32);
  domParams.minSeparationSec = 0.0;

  state.forecastDominantTransitions =
      stellar::sim::detectDominantBodyTransitions(sys, timeDays, state.forecastSamples, pred.gravity, domParams);

  state.forecastValid = true;

  char buf[256];
  std::snprintf(buf,
                sizeof(buf),
                "Forecast: %zu samples (horizon %.0f min, step %.2f s)",
                state.forecastSamples.size(),
                static_cast<double>(state.forecastHorizonMin),
                stepSec);
  state.forecastStatus = buf;
}

} // namespace

void drawOrbitAnalyzerWindow(OrbitAnalyzerWindowState& state,
                             const OrbitAnalyzerBindings& bindings,
                             const stellar::sim::StarSystem& sys,
                             double timeDays,
                             const stellar::sim::Ship& ship,
                             const stellar::sim::GravityParams& gravityParams) {
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

  // Planner gravity scaling: pass a gravityScale multiplier so headless planners
  // produce numbers consistent with the game's effective gravity.
  const double plannerGravityScale = state.useGravityScale ? gravityParams.scale : 1.0;

  if (orb.type == stellar::sim::TwoBodyOrbit::Type::Elliptic) {
    // --- Circularize at apses ---
    if (ImGui::Button("Set node: circularize @ peri")) {
      const auto res = stellar::sim::planCircularize(ship,
                                                     timeDays,
                                                     ref.body,
                                                     stellar::sim::ManeuverProgramKind::CircularizeAtPeriapsis,
                                                     plannerGravityScale);
      writePlanToNodeControls(bindings, res);
    }
    ImGui::SameLine();
    if (ImGui::Button("Set node: circularize @ apo")) {
      const auto res = stellar::sim::planCircularize(ship,
                                                     timeDays,
                                                     ref.body,
                                                     stellar::sim::ManeuverProgramKind::CircularizeAtApoapsis,
                                                     plannerGravityScale);
      writePlanToNodeControls(bindings, res);
    }

    ImGui::Spacing();

    // --- Adjust apoapsis (burn at periapsis) ---
    {
      const double raTargetReq = ref.body.radiusKm + (double)state.targetApoAltKm;
      const double raTarget = std::max(orb.periapsisKm, raTargetReq);

      if (ImGui::Button("Set node: burn @ peri to target apo")) {
        const auto res = stellar::sim::planSetApoapsisAtPeriapsis(ship,
                                                                  timeDays,
                                                                  ref.body,
                                                                  raTarget,
                                                                  plannerGravityScale);
        writePlanToNodeControls(bindings, res);
      }
      ImGui::SameLine();
      ImGui::TextDisabled("(apo target %.0f km)", raTarget);
    }

    // --- Adjust periapsis (burn at apoapsis) ---
    {
      const double rpTargetReq = ref.body.radiusKm + (double)state.targetPeriAltKm;
      const double rpTarget = std::min(orb.apoapsisKm, rpTargetReq);

      if (ImGui::Button("Set node: burn @ apo to target peri")) {
        const auto res = stellar::sim::planSetPeriapsisAtApoapsis(ship,
                                                                  timeDays,
                                                                  ref.body,
                                                                  rpTarget,
                                                                  plannerGravityScale);
        writePlanToNodeControls(bindings, res);
      }
      ImGui::SameLine();
      ImGui::TextDisabled("(peri target %.0f km)", rpTarget);
    }

    ImGui::Spacing();

    // Escape at current radius (instantaneous)
    {
      const auto esc = stellar::sim::planEscapeNow(ship, timeDays, ref.body, plannerGravityScale);
      const double dvMs = esc.valid ? (esc.dvKmS * 1000.0) : 0.0;

      if (ImGui::Button("Set node: escape now")) {
        writePlanToNodeControls(bindings, esc);
      }
      ImGui::SameLine();
      if (esc.valid) {
        ImGui::TextDisabled("(dv %.1f m/s)", dvMs);
      } else {
        ImGui::TextDisabled("(unavailable)");
      }
    }

    ImGui::Spacing();
    ImGui::SeparatorText("Plane change (ref plane)");
    ImGui::Checkbox("Force prograde align", &state.planeAlignForcePrograde);

    // Plane align is valid for non-equatorial bound orbits.
    const auto asc = stellar::sim::planAlignPlaneAtAscendingNode(ship,
                                                                 timeDays,
                                                                 ref.body,
                                                                 state.planeAlignForcePrograde,
                                                                 plannerGravityScale);

    const auto desc = stellar::sim::planAlignPlaneAtDescendingNode(ship,
                                                                   timeDays,
                                                                   ref.body,
                                                                   state.planeAlignForcePrograde,
                                                                   plannerGravityScale);

    if (asc.valid) {
      ImGui::Text("Ascending node: in %s | alt %.0f km | dv %.1f m/s",
                  formatTime(asc.timeToNodeSec).c_str(),
                  asc.nodeAltitudeKm,
                  asc.dvKmS * 1000.0);
    } else {
      ImGui::TextDisabled("Ascending node: (unavailable)");
    }

    ImGui::BeginDisabled(!asc.valid);
    if (ImGui::Button("Set node: align plane @ AN")) {
      writePlanToNodeControls(bindings, asc);
    }
    ImGui::EndDisabled();

    if (desc.valid) {
      ImGui::Text("Descending node: in %s | alt %.0f km | dv %.1f m/s",
                  formatTime(desc.timeToNodeSec).c_str(),
                  desc.nodeAltitudeKm,
                  desc.dvKmS * 1000.0);
    } else {
      ImGui::TextDisabled("Descending node: (unavailable)");
    }

    ImGui::BeginDisabled(!desc.valid);
    if (ImGui::Button("Set node: align plane @ DN")) {
      writePlanToNodeControls(bindings, desc);
    }
    ImGui::EndDisabled();

  } else {
    ImGui::TextDisabled("Planner actions require a bound (elliptic) orbit.");
  }

  ImGui::Spacing();
  ImGui::SeparatorText("Forecast");

  bool doRecomputeForecast = ImGui::Button("Recompute##orbitForecast");
  ImGui::SameLine();
  ImGui::Checkbox("Auto##orbitForecast", &state.forecastAutoUpdate);
  ImGui::SameLine();
  ImGui::SetNextItemWidth(90.0f);
  ImGui::InputFloat("Interval (sec)##orbitForecast", &state.forecastAutoUpdateIntervalSec, 0.5f, 2.0f, "%.1f");

  ImGui::SetNextItemWidth(120.0f);
  ImGui::InputFloat("Horizon (min)##orbitForecast", &state.forecastHorizonMin, 5.0f, 30.0f, "%.0f");
  ImGui::SameLine();
  ImGui::SetNextItemWidth(120.0f);
  ImGui::InputFloat("Step (sec)##orbitForecast", &state.forecastStepSec, 1.0f, 10.0f, "%.1f");

  ImGui::SetNextItemWidth(120.0f);
  ImGui::InputInt("Max samples##orbitForecast", &state.forecastMaxSamples);
  ImGui::SameLine();
  ImGui::SetNextItemWidth(120.0f);
  ImGui::InputInt("Refine depth##orbitForecast", &state.forecastRefineDepth);

  state.forecastHorizonMin = std::max(0.0f, state.forecastHorizonMin);
  state.forecastStepSec = std::max(0.01f, state.forecastStepSec);
  state.forecastAutoUpdateIntervalSec = std::max(0.1f, state.forecastAutoUpdateIntervalSec);
  state.forecastMaxSamples = std::max(2, state.forecastMaxSamples);
  state.forecastRefineDepth = std::clamp(state.forecastRefineDepth, 0, 32);

  if (state.forecastAutoUpdate) {
    const double lastDays = state.forecastLastComputeTimeDays;
    const double dtSec = (lastDays < 0.0) ? 1e30 : (timeDays - lastDays) * stellar::sim::kSecondsPerDay;
    if (dtSec >= static_cast<double>(state.forecastAutoUpdateIntervalSec)) {
      doRecomputeForecast = true;
    }
  }

  if (doRecomputeForecast) {
    recomputeForecast(state, sys, timeDays, ship, gravityParams);
  }

  if (!state.forecastStatus.empty()) {
    ImGui::TextDisabled("%s", state.forecastStatus.c_str());
  }

  if (state.forecastValid) {
    const auto& analysis = state.forecastAnalysis;

    if (analysis.firstImpact.valid) {
      const auto& hit = analysis.firstImpact;
      const double distKm = (hit.posKm - hit.body.posKm).length();
      const double altKm = distKm - hit.body.radiusKm;
      ImGui::Text("Impact: %s in %s | alt %.1f km",
                  bodyLabel(hit.body).c_str(),
                  formatTime(hit.tSec).c_str(),
                  altKm);
    } else {
      ImGui::TextDisabled("Impact: (none within horizon)");
    }

    if (ImGui::BeginTable("orbitForecastApproaches", 4,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
      ImGui::TableSetupColumn("Body");
      ImGui::TableSetupColumn("Closest in");
      ImGui::TableSetupColumn("Altitude (km)");
      ImGui::TableSetupColumn("Distance (km)");
      ImGui::TableHeadersRow();

      for (const auto& a : analysis.closestApproachByBody) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        if (!a.valid) {
          ImGui::TextDisabled("-");
          ImGui::TableSetColumnIndex(1);
          ImGui::TextDisabled("-");
          ImGui::TableSetColumnIndex(2);
          ImGui::TextDisabled("-");
          ImGui::TableSetColumnIndex(3);
          ImGui::TextDisabled("-");
          continue;
        }

        ImGui::TextUnformatted(bodyLabel(a.body).c_str());
        ImGui::TableSetColumnIndex(1);
        ImGui::TextUnformatted(formatTime(a.tSec).c_str());
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%.1f", a.altitudeKm);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.1f", a.distanceKm);
      }

      ImGui::EndTable();
    }

    if (ImGui::BeginTable("orbitForecastDominant", 3,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
      ImGui::TableSetupColumn("t");
      ImGui::TableSetupColumn("From");
      ImGui::TableSetupColumn("To");
      ImGui::TableHeadersRow();

      if (state.forecastDominantTransitions.empty()) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::TextDisabled("-");
        ImGui::TableSetColumnIndex(1);
        ImGui::TextDisabled("-");
        ImGui::TableSetColumnIndex(2);
        ImGui::TextDisabled("-");
      } else {
        for (const auto& e : state.forecastDominantTransitions) {
          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::TextUnformatted(formatTime(e.tSec).c_str());
          ImGui::TableSetColumnIndex(1);
          ImGui::TextUnformatted(bodyLabel(e.from).c_str());
          ImGui::TableSetColumnIndex(2);
          ImGui::TextUnformatted(bodyLabel(e.to).c_str());
        }
      }

      ImGui::EndTable();
    }
  } else {
    ImGui::TextDisabled("Forecast: click Recompute to generate a short-horizon prediction.");
  }


  ImGui::End();
}

} // namespace stellar::game
