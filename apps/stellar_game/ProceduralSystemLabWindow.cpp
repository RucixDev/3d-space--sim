#include "ProceduralSystemLabWindow.h"

#include "stellar/core/Log.h"
#include "stellar/core/Random.h"
#include "stellar/core/Hash.h"
#include "stellar/math/Math.h"
#include "stellar/sim/System.h"
#include "stellar/sim/Signals.h"
#include "stellar/sim/ResourceField.h"
#include "stellar/sim/Mission.h"
#include "stellar/sim/Universe.h"
#include "stellar/econ/Commodity.h"
#include "stellar/ui/Format.h"
#include "stellar/render/ProceduralRings.h"
#include "stellar/proc/AsteroidBeltGenerator.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <cmath>
#include <vector>

#include <imgui.h>

namespace stellar::game {
namespace {

static const char* starClassStr(sim::StarClass c) {
  switch (c) {
    case sim::StarClass::O: return "O";
    case sim::StarClass::B: return "B";
    case sim::StarClass::A: return "A";
    case sim::StarClass::F: return "F";
    case sim::StarClass::G: return "G";
    case sim::StarClass::K: return "K";
    case sim::StarClass::M: return "M";
    default: return "?";
  }
}

static const char* planetTypeStr(sim::PlanetType t) {
  switch (t) {
    case sim::PlanetType::Rocky: return "Rocky";
    case sim::PlanetType::Desert: return "Desert";
    case sim::PlanetType::Ocean: return "Ocean";
    case sim::PlanetType::Ice: return "Ice";
    case sim::PlanetType::GasGiant: return "Gas Giant";
    default: return "Unknown";
  }
}

// Solar masses per Earth mass.
constexpr double kEarthMassSol = 5.9722e24 / 1.98847e30;

static double hillRadiusAU(const sim::Planet& host, const sim::Star& star) {
  const double mStar = std::max(0.08, star.massSol);
  const double mP = std::max(0.0, host.massEarth * kEarthMassSol);
  if (mP <= 0.0) { return 0.0; }
  return host.orbit.semiMajorAxisAU * std::cbrt(mP / (3.0 * mStar));
}

} // namespace

void drawProceduralSystemLabWindow(ProceduralSystemLabWindowState& state,
                                  sim::Universe& universe,
                                  sim::SystemId currentSystemId,
                                  double timeDays,
                                  float /*timeSec*/,
                                  const ToastFn& toast) {
  if (!state.open) {
    return;
  }

  if (state.followCurrentSystem) {
    state.systemId = currentSystemId;
  }

  ImGui::SetNextWindowSize(ImVec2(880, 720), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Procedural System Lab", &state.open)) {
    ImGui::End();
    return;
  }

  ImGui::Checkbox("Follow current system", &state.followCurrentSystem);
  ImGui::SameLine();
  ImGui::BeginDisabled(state.followCurrentSystem);
  ImGui::SetNextItemWidth(260.0f);
  ImGui::InputScalar("System ID", ImGuiDataType_U64, &state.systemId, nullptr, nullptr, "%llu");
  ImGui::EndDisabled();

  sim::SystemId sid = state.systemId;
  if (sid == 0) {
    ImGui::TextDisabled("No system selected.");
    ImGui::End();
    return;
  }

  const sim::StarSystem& sys = universe.getSystem(sid);

  // Header
  ImGui::SeparatorText("System");
  ImGui::Text("Name: %s", sys.stub.name.c_str());
  ImGui::Text("SystemId: %s", ui::toString(sys.stub.id).c_str());
  ImGui::SameLine();
  ImGui::TextDisabled("Seed: %s", ui::toString(sys.stub.seed).c_str());
  ImGui::Text("Star: class %s  M=%.2f Msun  R=%.2f Rsun  L=%.2f Lsun  Teff=%d K",
              starClassStr(sys.star.cls),
              sys.star.massSol,
              sys.star.radiusSol,
              sys.star.luminositySol,
              sys.star.temperatureK);

  if (state.showZones) {
    const double L = std::max(0.000001, sys.star.luminositySol);
    const double earthEqAU = std::sqrt(L); // flux ~= 1 at sqrt(L)
    const double innerHZ = 0.75 * earthEqAU;
    const double outerHZ = 1.75 * earthEqAU;
    const double frostAU = 2.7 * earthEqAU;

    ImGui::TextDisabled("Approx zones: HZ [%.2f, %.2f] AU   Frost line ~ %.2f AU   (scaled by sqrt(L))",
                        innerHZ, outerHZ, frostAU);
  }

  ImGui::SeparatorText("Bodies");
  ImGui::Checkbox("Planets", &state.showPlanets);
  ImGui::SameLine();
  ImGui::Checkbox("Moons", &state.showMoons);
  ImGui::SameLine();
  ImGui::Checkbox("Signals", &state.showSignals);
  ImGui::SameLine();
  ImGui::Checkbox("Rings", &state.showRings);
  ImGui::SameLine();
  ImGui::Checkbox("Belts", &state.showBelts);

  // Precompute moon counts per planet
  std::vector<int> moonCount(sys.planets.size(), 0);
  for (const auto& m : sys.moons) {
    if (m.parentPlanetIndex < moonCount.size()) {
      moonCount[m.parentPlanetIndex] += 1;
    }
  }

  if (state.showPlanets) {
    ImGui::SeparatorText("Planets");
    const ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable |
                                  ImGuiTableFlags_ScrollY;
    if (ImGui::BeginTable("##planet_table", 9, flags, ImVec2(0.0f, 260.0f))) {
      ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 24.0f);
      ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthStretch);
      ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed, 86.0f);
      ImGui::TableSetupColumn("a (AU)", ImGuiTableColumnFlags_WidthFixed, 74.0f);
      ImGui::TableSetupColumn("Period (d)", ImGuiTableColumnFlags_WidthFixed, 84.0f);
      ImGui::TableSetupColumn("e", ImGuiTableColumnFlags_WidthFixed, 44.0f);
      ImGui::TableSetupColumn("i (deg)", ImGuiTableColumnFlags_WidthFixed, 64.0f);
      ImGui::TableSetupColumn("R (Re)", ImGuiTableColumnFlags_WidthFixed, 64.0f);
      ImGui::TableSetupColumn("Moons", ImGuiTableColumnFlags_WidthFixed, 54.0f);
      ImGui::TableHeadersRow();

      for (std::size_t i = 0; i < sys.planets.size(); ++i) {
        const auto& p = sys.planets[i];
        ImGui::TableNextRow();

        ImGui::TableSetColumnIndex(0);
        ImGui::Text("%d", (int)i);

        ImGui::TableSetColumnIndex(1);
        ImGui::TextUnformatted(p.name.c_str());

        ImGui::TableSetColumnIndex(2);
        ImGui::TextUnformatted(planetTypeStr(p.type));

        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%.2f", p.orbit.semiMajorAxisAU);

        ImGui::TableSetColumnIndex(4);
        ImGui::Text("%.0f", p.orbit.periodDays);

        ImGui::TableSetColumnIndex(5);
        ImGui::Text("%.2f", p.orbit.eccentricity);

        ImGui::TableSetColumnIndex(6);
        ImGui::Text("%.1f", p.orbit.inclinationRad * 57.29577951308232);

        ImGui::TableSetColumnIndex(7);
        ImGui::Text("%.2f", p.radiusEarth);

        ImGui::TableSetColumnIndex(8);
        ImGui::Text("%d", moonCount[i]);
      }

      ImGui::EndTable();
    }
  }

  if (state.showBelts) {
    ImGui::SeparatorText("Minor Bodies");

    ImGui::TextDisabled("Procedural asteroid belts / debris disks with resonance gaps + trojan swarms.");

    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("Points", &state.beltPointCount);
    state.beltPointCount = std::clamp(state.beltPointCount, 0, 60000);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("Candidates/pt", &state.beltCandidatesPerPoint);
    state.beltCandidatesPerPoint = std::clamp(state.beltCandidatesPerPoint, 1, 64);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("Plot res", &state.beltRadialPlotRes);
    state.beltRadialPlotRes = std::clamp(state.beltRadialPlotRes, 32, 512);

    ImGui::Checkbox("Scatter", &state.beltShowScatter);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("Max scatter pts", &state.beltScatterMaxPoints);
    state.beltScatterMaxPoints = std::clamp(state.beltScatterMaxPoints, 0, 20000);
    ImGui::SameLine();
    ImGui::Checkbox("Color by density", &state.beltScatterColorByDensity);
    ImGui::SameLine();
    ImGui::Checkbox("Resonance rings", &state.beltShowResonanceRings);

    const auto beltPlan = proc::generateAsteroidBelts(universe.seed(), sys);

    if (beltPlan.belts.empty()) {
      ImGui::TextDisabled("(no belts)");
    } else {
      state.selectedBelt = std::clamp(state.selectedBelt, 0, (int)beltPlan.belts.size() - 1);

      // Belt list with selection.
      if (ImGui::BeginTable("##belt_table", 7,
                            ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable,
                            ImVec2(0.0f, 140.0f))) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 24.0f);
        ImGui::TableSetupColumn("Kind", ImGuiTableColumnFlags_WidthFixed, 108.0f);
        ImGui::TableSetupColumn("Inner (AU)", ImGuiTableColumnFlags_WidthFixed, 86.0f);
        ImGui::TableSetupColumn("Outer (AU)", ImGuiTableColumnFlags_WidthFixed, 86.0f);
        ImGui::TableSetupColumn("Thickness", ImGuiTableColumnFlags_WidthFixed, 86.0f);
        ImGui::TableSetupColumn("Control", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Res", ImGuiTableColumnFlags_WidthFixed, 42.0f);
        ImGui::TableHeadersRow();

        for (int i = 0; i < (int)beltPlan.belts.size(); ++i) {
          const auto& b = beltPlan.belts[(std::size_t)i];
          ImGui::TableNextRow();

          ImGui::TableSetColumnIndex(0);
          const bool sel = (i == state.selectedBelt);
          char buf[16];
          std::snprintf(buf, sizeof(buf), "%d", i);
          if (ImGui::Selectable(buf, sel, ImGuiSelectableFlags_SpanAllColumns)) {
            state.selectedBelt = i;
            state.beltCacheKey = 0; // invalidate
          }

          ImGui::TableSetColumnIndex(1);
          ImGui::TextUnformatted(proc::asteroidBeltKindName(b.kind));

          ImGui::TableSetColumnIndex(2);
          ImGui::Text("%.2f", b.aInnerAU);

          ImGui::TableSetColumnIndex(3);
          ImGui::Text("%.2f", b.aOuterAU);

          ImGui::TableSetColumnIndex(4);
          ImGui::Text("%.3f", b.thicknessAU);

          ImGui::TableSetColumnIndex(5);
          if (b.controllingPlanetIndex >= 0 && b.controllingPlanetIndex < (int)sys.planets.size()) {
            ImGui::Text("%s", sys.planets[(std::size_t)b.controllingPlanetIndex].name.c_str());
          } else {
            ImGui::TextDisabled("-");
          }

          ImGui::TableSetColumnIndex(6);
          ImGui::Text("%d", (int)b.resonances.size());
        }

        ImGui::EndTable();
      }

      const auto& b = beltPlan.belts[(std::size_t)state.selectedBelt];

      // Cache key includes the selected belt + sampling knobs.
      core::u64 key = 0;
      key = core::hashCombine(key, sys.stub.id);
      key = core::hashCombine(key, sys.stub.seed);
      key = core::hashCombine(key, b.id);
      key = core::hashCombine(key, (core::u64)state.selectedBelt);
      key = core::hashCombine(key, (core::u64)state.beltPointCount);
      key = core::hashCombine(key, (core::u64)state.beltCandidatesPerPoint);
      key = core::hashCombine(key, (core::u64)state.beltRadialPlotRes);
      key = core::hashCombine(key, (core::u64)state.beltScatterMaxPoints);

      if (key != state.beltCacheKey || state.beltCacheSelected != state.selectedBelt) {
        state.beltCacheKey = key;
        state.beltCacheSelected = state.selectedBelt;

        // --- Radial mean density plot ---
        state.beltRadialMean01.assign((std::size_t)state.beltRadialPlotRes, 0.0f);
        const int nTheta = 24;
        for (int i = 0; i < state.beltRadialPlotRes; ++i) {
          const double u = (state.beltRadialPlotRes <= 1) ? 0.0 : (double)i / (double)(state.beltRadialPlotRes - 1);
          const double aAU = b.aInnerAU + (b.aOuterAU - b.aInnerAU) * u;

          double sum = 0.0;
          for (int t = 0; t < nTheta; ++t) {
            const double th = (double)t / (double)nTheta * (2.0 * stellar::math::kPi);
            sum += proc::asteroidBeltDensity01(b, aAU, th);
          }
          const double mean = sum / (double)nTheta;
          state.beltRadialMean01[(std::size_t)i] = (float)std::clamp(mean, 0.0, 1.0);
        }

        // --- Scatter points ---
        state.beltScatterPosAU.clear();
        state.beltScatterDensity01.clear();
        const int nPts = std::min(state.beltPointCount, state.beltScatterMaxPoints);
        if (state.beltShowScatter && nPts > 0) {
          const auto pts = proc::sampleAsteroidBeltPoints(universe.seed(), sys, b, nPts, state.beltCandidatesPerPoint);
          state.beltScatterPosAU.reserve(pts.size());
          state.beltScatterDensity01.reserve(pts.size());

          for (const auto& p : pts) {
            // Store belt-local coordinates (x along basisX, y along basisY, z along basisZ).
            const float x = (float)math::dot(p.posAU, b.basisX);
            const float y = (float)math::dot(p.posAU, b.basisY);
            const float z = (float)math::dot(p.posAU, b.basisZ);
            state.beltScatterPosAU.push_back({x, y, z});
            state.beltScatterDensity01.push_back((float)p.density01);
          }
        }
      }

      // Summary line.
      ImGui::Text("Selected: %s  span [%.2f, %.2f] AU  control=%s",
                  proc::asteroidBeltKindName(b.kind),
                  b.aInnerAU,
                  b.aOuterAU,
                  (b.controllingPlanetIndex >= 0 && b.controllingPlanetIndex < (int)sys.planets.size())
                      ? sys.planets[(std::size_t)b.controllingPlanetIndex].name.c_str()
                      : "-"
      );

      if (!state.beltRadialMean01.empty()) {
        ImGui::PlotLines("Radial mean density",
                         state.beltRadialMean01.data(),
                         (int)state.beltRadialMean01.size(),
                         0,
                         nullptr,
                         0.0f,
                         1.0f,
                         ImVec2(0.0f, 78.0f));
      }

      if (!b.resonances.empty()) {
        ImGui::TextDisabled("Resonances:");
        for (const auto& f : b.resonances) {
          const bool ridge = (f.strength01 < 0.0);
          if (ridge) {
            ImGui::BulletText("%d:%d @ %.2f AU  ridge  width %.3f AU  amp %.2f",
                              f.m, f.n, f.aAU, f.halfWidthAU, std::abs(f.strength01));
          } else {
            ImGui::BulletText("%d:%d @ %.2f AU  gap    width %.3f AU  depth %.2f",
                              f.m, f.n, f.aAU, f.halfWidthAU, f.strength01);
          }
        }
      }

      if (state.beltShowScatter && !state.beltScatterPosAU.empty()) {
        ImGui::SeparatorText("Belt Scatter");

        // 2D plot in belt plane (x,z). y is thickness.
        const ImVec2 plotSize = ImVec2(0.0f, 240.0f);
        ImGui::BeginChild("##belt_scatter", plotSize, true, ImGuiWindowFlags_NoScrollbar);

        ImDrawList* dl = ImGui::GetWindowDrawList();
        const ImVec2 p0 = ImGui::GetCursorScreenPos();
        const ImVec2 avail = ImGui::GetContentRegionAvail();
        const float wpx = std::max(10.0f, avail.x);
        const float hpx = std::max(10.0f, avail.y);
        const ImVec2 p1 = ImVec2(p0.x + wpx, p0.y + hpx);

        // Background
        dl->AddRectFilled(p0, p1, ImGui::GetColorU32(ImVec4(0.05f, 0.05f, 0.07f, 1.0f)));
        dl->AddRect(p0, p1, ImGui::GetColorU32(ImGuiCol_Border));

        const float rAU = (float)std::max(0.1, b.aOuterAU) * 1.05f;
        auto toScreen = [&](float xAU, float zAU) -> ImVec2 {
          const float nx = (xAU / rAU) * 0.5f + 0.5f;
          const float ny = (zAU / rAU) * 0.5f + 0.5f;
          // Flip Y so +Z is up.
          return ImVec2(p0.x + nx * wpx, p1.y - ny * hpx);
        };

        // Resonance rings (gaps/ridges).
        if (state.beltShowResonanceRings && !b.resonances.empty()) {
          for (const auto& f : b.resonances) {
            const float rr = (float)f.aAU;
            const float radPx = (rr / rAU) * 0.5f * std::min(wpx, hpx);
            const bool ridge = (f.strength01 < 0.0);
            const ImU32 col = ridge ? ImGui::GetColorU32(ImVec4(0.30f, 0.75f, 0.40f, 0.60f))
                                    : ImGui::GetColorU32(ImVec4(0.85f, 0.35f, 0.25f, 0.55f));
            dl->AddCircle(ImVec2(p0.x + wpx * 0.5f, p0.y + hpx * 0.5f), radPx, col, 72, 1.0f);
          }
        }

        // Draw points.
        const float baseSize = 1.45f;
        for (std::size_t i = 0; i < state.beltScatterPosAU.size(); ++i) {
          const auto& v = state.beltScatterPosAU[i];
          const float x = v[0];
          const float z = v[2];
          const float dens = (i < state.beltScatterDensity01.size()) ? state.beltScatterDensity01[i] : 1.0f;

          ImU32 col = ImGui::GetColorU32(ImGuiCol_Text);
          if (state.beltScatterColorByDensity) {
            const float d = std::clamp(dens, 0.0f, 1.0f);
            col = ImGui::GetColorU32(ImVec4(0.20f + 0.80f * d, 0.25f + 0.75f * d, 0.35f + 0.65f * d, 1.0f));
          }

          const ImVec2 sp = toScreen(x, z);
          dl->AddCircleFilled(sp, baseSize, col, 6);
        }

        // Keep cursor/child size stable.
        ImGui::Dummy(avail);
        ImGui::EndChild();
      }
    }
  }


  if (state.showRings) {
    ImGui::SeparatorText("Rings");

    // Controls mirror the world ring logic in apps/stellar_game/main.cpp (chance by planet type).
    state.ringChanceMul = std::clamp(state.ringChanceMul, 0.0, 3.0);
    state.ringOpacity = std::clamp(state.ringOpacity, 0.0, 2.0);

    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputDouble("Chance mul", &state.ringChanceMul, 0.05, 0.25, "%.2f");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputDouble("Opacity", &state.ringOpacity, 0.05, 0.25, "%.2f");

    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("Preview W", &state.ringPreviewW);
    ImGui::SameLine();
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("Preview H", &state.ringPreviewH);
    state.ringPreviewW = std::clamp(state.ringPreviewW, 64, 1024);
    state.ringPreviewH = std::clamp(state.ringPreviewH, 32, 512);

    ImGui::Checkbox("Alpha-only", &state.ringPreviewAlphaOnly);
    ImGui::SameLine();
    ImGui::Checkbox("Use auto variant", &state.ringUseAutoVariant);

    // Table: predicted rings per planet.
    if (ImGui::BeginTable("##ring_planets", 7,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable,
                          ImVec2(0.0f, 150.0f))) {
      ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 24.0f);
      ImGui::TableSetupColumn("Planet", ImGuiTableColumnFlags_WidthStretch);
      ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed, 86.0f);
      ImGui::TableSetupColumn("Chance", ImGuiTableColumnFlags_WidthFixed, 72.0f);
      ImGui::TableSetupColumn("Has", ImGuiTableColumnFlags_WidthFixed, 36.0f);
      ImGui::TableSetupColumn("Var", ImGuiTableColumnFlags_WidthFixed, 36.0f);
      ImGui::TableSetupColumn("Seed", ImGuiTableColumnFlags_WidthStretch);
      ImGui::TableHeadersRow();

      auto baseRingChance = [](sim::PlanetType t) -> double {
        switch (t) {
          case sim::PlanetType::GasGiant: return 0.30;
          case sim::PlanetType::Ice:      return 0.22;
          case sim::PlanetType::Ocean:    return 0.18;
          case sim::PlanetType::Desert:   return 0.14;
          case sim::PlanetType::Rocky:    return 0.12;
          default:                        return 0.10;
        }
      };

      for (int i = 0; i < (int)sys.planets.size(); ++i) {
        const auto& p = sys.planets[(std::size_t)i];

        const core::u64 pSurfaceSeed = core::hashCombine(
            core::hashCombine(core::fnv1a64("planet_surface"), (core::u64)sys.stub.seed),
            (core::u64)i);

        const core::u64 ringSeedBase = core::hashCombine(pSurfaceSeed, core::fnv1a64("rings"));
        core::SplitMix64 rrng(ringSeedBase);

        const double chance = baseRingChance(p.type) * state.ringChanceMul;
        bool has = false;
        int var = 0;
        if (chance > 1e-6 && rrng.chance(chance)) {
          has = true;
          var = (int)rrng.range(0, 2);
        }

        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::Text("%d", i);
        ImGui::TableNextColumn();
        if (ImGui::Selectable(p.name.c_str(), state.ringPlanetIndex == i, ImGuiSelectableFlags_SpanAllColumns)) {
          state.ringPlanetIndex = i;
        }
        ImGui::TableNextColumn();
        ImGui::TextUnformatted(planetTypeStr(p.type));
        ImGui::TableNextColumn();
        ImGui::Text("%.2f", chance);
        ImGui::TableNextColumn();
        ImGui::TextUnformatted(has ? "Y" : "N");
        ImGui::TableNextColumn();
        ImGui::Text("%d", has ? var : -1);
        ImGui::TableNextColumn();
        ImGui::TextDisabled("%s", ui::toString(ringSeedBase).c_str());
      }

      ImGui::EndTable();
    }

    // Selected planet preview.
    if (!sys.planets.empty()) {
      state.ringPlanetIndex = std::clamp(state.ringPlanetIndex, 0, (int)sys.planets.size() - 1);
      const int i = state.ringPlanetIndex;
      const auto& p = sys.planets[(std::size_t)i];

      const core::u64 pSurfaceSeed = core::hashCombine(
          core::hashCombine(core::fnv1a64("planet_surface"), (core::u64)sys.stub.seed),
          (core::u64)i);
      const core::u64 ringSeedBase = core::hashCombine(pSurfaceSeed, core::fnv1a64("rings"));
      core::SplitMix64 rrng(ringSeedBase);

      auto baseRingChance = [](sim::PlanetType t) -> double {
        switch (t) {
          case sim::PlanetType::GasGiant: return 0.30;
          case sim::PlanetType::Ice:      return 0.22;
          case sim::PlanetType::Ocean:    return 0.18;
          case sim::PlanetType::Desert:   return 0.14;
          case sim::PlanetType::Rocky:    return 0.12;
          default:                        return 0.10;
        }
      };

      const double chance = baseRingChance(p.type) * state.ringChanceMul;
      bool has = false;
      int autoVar = 0;
      double autoAlphaMul = 0.0;
      if (chance > 1e-6 && rrng.chance(chance)) {
        has = true;
        autoVar = (int)rrng.range(0, 2);
        autoAlphaMul = state.ringOpacity * rrng.range(0.75, 1.10);
      }

      ImGui::Separator();
      ImGui::Text("Selected: #%d %s (%s)", i, p.name.c_str(), planetTypeStr(p.type));
      ImGui::Text("Has rings: %s  autoVar=%d  alphaMul=%.2f", has ? "yes" : "no", autoVar, autoAlphaMul);

      ImGui::SetNextItemWidth(120.0f);
      ImGui::SliderInt("Variant", &state.ringVariant, 0, 2);

      const int variant = (state.ringUseAutoVariant && has) ? autoVar : std::clamp(state.ringVariant, 0, 2);
      const core::u64 ringSeed = core::hashCombine(ringSeedBase, (core::u64)variant);

      // Build/refresh cached preview.
      if (state.ringPreviewSeed != ringSeed || state.ringPreviewWCache != state.ringPreviewW ||
          state.ringPreviewHCache != state.ringPreviewH || state.ringPreviewRGBA.empty()) {

        const render::RingImage img = render::generateRingTexture(ringSeed, state.ringPreviewW, state.ringPreviewH);
        state.ringPreviewRGBA = img.rgba;
        state.ringPreviewSeed = ringSeed;
        state.ringPreviewWCache = img.w;
        state.ringPreviewHCache = img.h;

        // Radial mean density (alpha averaged over u).
        state.ringRadialMean01.assign((std::size_t)img.h, 0.0f);
        for (int y = 0; y < img.h; ++y) {
          double sum = 0.0;
          for (int x = 0; x < img.w; ++x) {
            sum += (double)img.rgba[(std::size_t)((y * img.w + x) * 4 + 3)] / 255.0;
          }
          state.ringRadialMean01[(std::size_t)y] = (float)(sum / std::max(1.0, (double)img.w));
        }
      }

      ImGui::TextDisabled("ringSeed: %s", ui::toString(state.ringPreviewSeed).c_str());

      // Plot the radial profile first (cheap).
      if (!state.ringRadialMean01.empty()) {
        ImGui::PlotLines("Radial mean density", state.ringRadialMean01.data(),
                         (int)state.ringRadialMean01.size(), 0, nullptr, 0.0f, 1.0f, ImVec2(0.0f, 64.0f));
      }

      // Draw a pixel preview using ImDrawList rectangles (CPU-side, no GL texture needed).
      const int pw = state.ringPreviewWCache;
      const int ph = state.ringPreviewHCache;
      if (pw > 0 && ph > 0 && state.ringPreviewRGBA.size() == (std::size_t)(pw * ph * 4)) {
        const float canvasW = std::min(560.0f, ImGui::GetContentRegionAvail().x);
        const float canvasH = canvasW * ((float)ph / (float)pw);
        const ImVec2 p0 = ImGui::GetCursorScreenPos();
        const ImVec2 p1 = ImVec2(p0.x + canvasW, p0.y + canvasH);
        ImGui::InvisibleButton("##ring_canvas", ImVec2(canvasW, canvasH));

        ImDrawList* dl = ImGui::GetWindowDrawList();
        dl->AddRectFilled(p0, p1, IM_COL32(10, 10, 10, 255));

        const float sx = canvasW / (float)pw;
        const float sy = canvasH / (float)ph;

        for (int y = 0; y < ph; ++y) {
          for (int x = 0; x < pw; ++x) {
            const std::size_t idx = (std::size_t)((y * pw + x) * 4);
            std::uint8_t R = state.ringPreviewRGBA[idx + 0];
            std::uint8_t G = state.ringPreviewRGBA[idx + 1];
            std::uint8_t B = state.ringPreviewRGBA[idx + 2];
            std::uint8_t A = state.ringPreviewRGBA[idx + 3];

            if (state.ringPreviewAlphaOnly) {
              R = G = B = A;
              A = 255;
            }

            const ImU32 col = IM_COL32(R, G, B, A);
            const ImVec2 q0(p0.x + (float)x * sx, p0.y + (float)y * sy);
            const ImVec2 q1(q0.x + sx, q0.y + sy);
            dl->AddRectFilled(q0, q1, col);
          }
        }

        dl->AddRect(p0, p1, IM_COL32(200, 200, 200, 120));
      }
    }
  }

  if (state.showMoons) {
    ImGui::SeparatorText("Moons");
    if (sys.moons.empty()) {
      ImGui::TextDisabled("(no moons in this system)");
    } else {
      const ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable |
                                    ImGuiTableFlags_ScrollY;
      if (ImGui::BeginTable("##moon_table", 12, flags, ImVec2(0.0f, 0.0f))) {
        ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableSetupColumn("Parent", ImGuiTableColumnFlags_WidthFixed, 150.0f);
        ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed, 86.0f);
        ImGui::TableSetupColumn("a (AU)", ImGuiTableColumnFlags_WidthFixed, 74.0f);
        ImGui::TableSetupColumn("Period (d)", ImGuiTableColumnFlags_WidthFixed, 84.0f);
        ImGui::TableSetupColumn("P ratio", ImGuiTableColumnFlags_WidthFixed, 92.0f);
        ImGui::TableSetupColumn("dRH", ImGuiTableColumnFlags_WidthFixed, 56.0f);
        ImGui::TableSetupColumn("e", ImGuiTableColumnFlags_WidthFixed, 44.0f);
        ImGui::TableSetupColumn("i (deg)", ImGuiTableColumnFlags_WidthFixed, 64.0f);
        ImGui::TableSetupColumn("R (Re)", ImGuiTableColumnFlags_WidthFixed, 64.0f);
        ImGui::TableSetupColumn("Hill%", ImGuiTableColumnFlags_WidthFixed, 62.0f);
        ImGui::TableSetupColumn("ID", ImGuiTableColumnFlags_WidthFixed, 110.0f);
        ImGui::TableHeadersRow();

        struct Res { int num; int den; double r; };
        const std::array<Res, 6> kRes = {Res{4,3,4.0/3.0}, Res{3,2,1.5}, Res{5,3,5.0/3.0},
                                        Res{2,1,2.0}, Res{5,2,2.5}, Res{3,1,3.0}};

        auto mutualHillAU = [&](double a1, double a2, double m1, double m2, double hostMassEarth) {
          const double M = std::max(1.0e-6, hostMassEarth);
          const double mu = std::max(1.0e-12, (m1 + m2) / (3.0 * M));
          const double aBar = 0.5 * (a1 + a2);
          return std::cbrt(mu) * std::max(1.0e-12, aBar);
        };

        for (std::size_t mi = 0; mi < sys.moons.size(); ++mi) {
          const auto& m = sys.moons[mi];
          const sim::Moon* prev = nullptr;
          if (mi > 0 && sys.moons[mi - 1].parentPlanetIndex == m.parentPlanetIndex) {
            prev = &sys.moons[mi - 1];
          }

          const sim::Planet* parent = nullptr;
          if (m.parentPlanetIndex < sys.planets.size()) {
            parent = &sys.planets[m.parentPlanetIndex];
          }

          const double h = parent ? hillRadiusAU(*parent, sys.star) : 0.0;
          const double frac = (h > 0.0) ? (m.orbit.semiMajorAxisAU / h) : 0.0;

          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::TextUnformatted(m.name.c_str());

          ImGui::TableSetColumnIndex(1);
          if (parent) {
            ImGui::Text("%d: %s", (int)m.parentPlanetIndex, parent->name.c_str());
          } else {
            ImGui::TextDisabled("(invalid)");
          }

          ImGui::TableSetColumnIndex(2);
          ImGui::TextUnformatted(planetTypeStr(m.type));

          ImGui::TableSetColumnIndex(3);
          ImGui::Text("%.4f", m.orbit.semiMajorAxisAU);

          ImGui::TableSetColumnIndex(4);
          ImGui::Text("%.2f", m.orbit.periodDays);

          ImGui::TableSetColumnIndex(5);
          if (prev) {
            const double pr = m.orbit.periodDays / std::max(1.0e-9, prev->orbit.periodDays);
            Res best = kRes[0];
            double bestErr = 1.0e9;
            for (const auto& rr : kRes) {
              const double err = std::abs(pr - rr.r) / rr.r;
              if (err < bestErr) {
                bestErr = err;
                best = rr;
              }
            }
            if (bestErr < 0.03) {
              ImGui::Text("%.2f (%d:%d)", pr, best.num, best.den);
            } else {
              ImGui::Text("%.2f", pr);
            }
          } else {
            ImGui::TextDisabled("-");
          }

          ImGui::TableSetColumnIndex(6);
          if (prev && parent) {
            const double rH = mutualHillAU(prev->orbit.semiMajorAxisAU, m.orbit.semiMajorAxisAU, prev->massEarth, m.massEarth, parent->massEarth);
            const double sep = (m.orbit.semiMajorAxisAU - prev->orbit.semiMajorAxisAU) / std::max(1.0e-12, rH);
            ImGui::Text("%.1f", sep);
          } else {
            ImGui::TextDisabled("-");
          }

          ImGui::TableSetColumnIndex(7);
          ImGui::Text("%.2f", m.orbit.eccentricity);

          ImGui::TableSetColumnIndex(8);
          ImGui::Text("%.1f", m.orbit.inclinationRad * 57.29577951308232);

          ImGui::TableSetColumnIndex(9);
          ImGui::Text("%.2f", m.radiusEarth);

          ImGui::TableSetColumnIndex(10);
          ImGui::Text("%.0f%%", frac * 100.0);

          ImGui::TableSetColumnIndex(11);
          ImGui::Text("%s", ui::toString(m.id).c_str());
        }

        ImGui::EndTable();
      }
    }
  }

  if (state.showSignals) {
    ImGui::SeparatorText("Signals");

    ImGui::TextDisabled("Deterministic daily procedural content (resource fields, derelicts, distress, traffic).");

    ImGui::Checkbox("Use current timeDays", &state.useCurrentTime);
    ImGui::SameLine();
    if (state.useCurrentTime) {
      ImGui::Text("t = %.2f d", timeDays);
    } else {
      ImGui::SetNextItemWidth(180.0f);
      ImGui::InputDouble("Override timeDays", &state.timeDaysOverride, 1.0, 10.0, "%.2f");
    }

    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputInt("Resource field count", &state.resourceFieldCount);
    state.resourceFieldCount = std::clamp(state.resourceFieldCount, 0, 12);

    ImGui::Checkbox("Daily derelict", &state.includeDailyDerelict);
    ImGui::SameLine();
    ImGui::Checkbox("Distress", &state.includeDistress);
    ImGui::SameLine();
    ImGui::Checkbox("Traffic convoys", &state.includeTrafficConvoys);

    const double tDays = state.useCurrentTime ? timeDays : state.timeDaysOverride;

    sim::SignalGenParams params{};
    params.resourceFieldCount = state.resourceFieldCount;
    params.includeDailyDerelict = state.includeDailyDerelict;
    params.includeDistress = state.includeDistress;
    params.includeTrafficConvoys = state.includeTrafficConvoys;
    params.distressPerDay = 2;
    params.distressTtlDays = 1.0;

    const std::vector<sim::Mission> activeMissions{};
    const std::vector<core::u64> resolvedSignalIds{};
    const sim::SystemSignalPlan plan =
        sim::generateSystemSignals(universe.seed(), sys, tDays, activeMissions, resolvedSignalIds, params, nullptr);

    ImGui::Text("Sites: %d   Resource fields: %d   Asteroids: %d",
                (int)plan.sites.size(),
                (int)plan.resourceFields.fields.size(),
                (int)plan.resourceFields.asteroids.size());

    // --- Resource Fields ---
    ImGui::SeparatorText("Resource Fields");

    if (plan.resourceFields.fields.empty()) {
      ImGui::TextDisabled("(no resource fields for current parameters)");
    } else {
      const auto& fields = plan.resourceFields.fields;
      const auto& asteroids = plan.resourceFields.asteroids;

      // Keep selection in-bounds.
      if (state.selectedResourceField >= (int)fields.size()) {
        state.selectedResourceField = -1;
      }

      const ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_Resizable |
                                    ImGuiTableFlags_ScrollY;
      if (ImGui::BeginTable("##rf_table", 11, flags, ImVec2(0.0f, 220.0f))) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 24.0f);
        ImGui::TableSetupColumn("Kind", ImGuiTableColumnFlags_WidthFixed, 90.0f);
        ImGui::TableSetupColumn("Layout", ImGuiTableColumnFlags_WidthFixed, 70.0f);
        ImGui::TableSetupColumn("Struct", ImGuiTableColumnFlags_WidthFixed, 68.0f);
        ImGui::TableSetupColumn("Major (km)", ImGuiTableColumnFlags_WidthFixed, 88.0f);
        ImGui::TableSetupColumn("Minor (km)", ImGuiTableColumnFlags_WidthFixed, 88.0f);
        ImGui::TableSetupColumn("Arc (deg)", ImGuiTableColumnFlags_WidthFixed, 72.0f);
        ImGui::TableSetupColumn("Rich", ImGuiTableColumnFlags_WidthFixed, 56.0f);
        ImGui::TableSetupColumn("Primary", ImGuiTableColumnFlags_WidthFixed, 76.0f);
        ImGui::TableSetupColumn("Secondary", ImGuiTableColumnFlags_WidthFixed, 86.0f);
        ImGui::TableSetupColumn("ID", ImGuiTableColumnFlags_WidthStretch);
        ImGui::TableHeadersRow();

        for (int i = 0; i < (int)fields.size(); ++i) {
          const auto& f = fields[i];
          ImGui::TableNextRow();

          ImGui::TableSetColumnIndex(0);
          const bool selected = (state.selectedResourceField == i);
          char idxLabel[16];
          std::snprintf(idxLabel, sizeof(idxLabel), "%d", i);
          if (ImGui::Selectable(idxLabel, selected, ImGuiSelectableFlags_SpanAllColumns)) {
            state.selectedResourceField = i;
          }

          ImGui::TableSetColumnIndex(1);
          ImGui::TextUnformatted(sim::resourceFieldKindName(f.kind));

          ImGui::TableSetColumnIndex(2);
          ImGui::TextUnformatted(sim::resourceFieldLayoutName(f.layout));

          ImGui::TableSetColumnIndex(3);
          {
            int nHot = 0, nGap = 0, nStreak = 0, nSpokes = 0;
            for (const auto& ft : plan.resourceFields.features) {
              if (ft.fieldId != f.id) continue;
              switch (ft.kind) {
                case sim::ResourceFieldFeatureKind::Hotspot: ++nHot; break;
                case sim::ResourceFieldFeatureKind::Gap: ++nGap; break;
                case sim::ResourceFieldFeatureKind::Streak: ++nStreak; break;
                case sim::ResourceFieldFeatureKind::Spokes: ++nSpokes; break;
              }
            }

            char buf[32];
            if (f.layout == sim::ResourceFieldLayout::Sheet) {
              std::snprintf(buf, sizeof(buf), "St%d", nStreak);
            } else if (f.layout == sim::ResourceFieldLayout::Cluster) {
              std::snprintf(buf, sizeof(buf), "H%d", nHot);
            } else {
              std::snprintf(buf, sizeof(buf), "H%d G%d%s", nHot, nGap, (nSpokes > 0) ? "+" : "");
            }
            ImGui::TextUnformatted(buf);
          }

          ImGui::TableSetColumnIndex(4);
          ImGui::Text("%.0f", f.majorRadiusKm);

          ImGui::TableSetColumnIndex(5);
          ImGui::Text("%.0f", f.minorRadiusKm);

          ImGui::TableSetColumnIndex(6);
          ImGui::Text("%.0f", f.arcRad * 57.29577951308232);

          ImGui::TableSetColumnIndex(7);
          ImGui::Text("%.2f", f.richness);

          ImGui::TableSetColumnIndex(8);
          const auto primCode = econ::commodityCode(f.primary);
          ImGui::TextUnformatted(primCode.data(), primCode.data() + primCode.size());

          ImGui::TableSetColumnIndex(9);
          const auto secCode = econ::commodityCode(f.secondary);
          ImGui::TextUnformatted(secCode.data(), secCode.data() + secCode.size());

          ImGui::TableSetColumnIndex(10);
          ImGui::Text("%s", ui::toString(f.id).c_str());
        }

        ImGui::EndTable();
      }

      if (state.selectedResourceField >= 0 && state.selectedResourceField < (int)fields.size()) {
        const auto& f = fields[state.selectedResourceField];

        ImGui::SeparatorText("Selected Field");
        ImGui::Text("Kind: %s   Layout: %s", sim::resourceFieldKindName(f.kind), sim::resourceFieldLayoutName(f.layout));
        ImGui::Text("Primary: %s   Secondary: %s",
                    econ::commodityDef(f.primary).name,
                    econ::commodityDef(f.secondary).name);
        ImGui::Text("Radii: major %.0f km   minor %.0f km   arc %.0f deg",
                    f.majorRadiusKm, f.minorRadiusKm, f.arcRad * 57.29577951308232);

        // Structural features for this field.
        const auto feats = sim::filterFeaturesForField(plan.resourceFields.features, f.id);
        {
          int nHot = 0, nGap = 0, nStreak = 0, nSpokes = 0;
          for (const auto& ft : feats) {
            switch (ft.kind) {
              case sim::ResourceFieldFeatureKind::Hotspot: ++nHot; break;
              case sim::ResourceFieldFeatureKind::Gap: ++nGap; break;
              case sim::ResourceFieldFeatureKind::Streak: ++nStreak; break;
              case sim::ResourceFieldFeatureKind::Spokes: ++nSpokes; break;
            }
          }
          ImGui::Text("Structure: hotspots %d   gaps %d   streaks %d   spokes %d", nHot, nGap, nStreak, nSpokes);
        }

        if (!feats.empty()) {
          ImGui::TextDisabled("Features:");
          for (const auto& ft : feats) {
            const double aDeg = ft.angleRad * 57.29577951308232;
            const double wDeg = ft.width * 57.29577951308232;
            switch (ft.kind) {
              case sim::ResourceFieldFeatureKind::Gap:
                ImGui::BulletText("Gap: center %.1f deg  sigma %.1f deg  depth %.2f", aDeg, wDeg, ft.strength01);
                break;
              case sim::ResourceFieldFeatureKind::Hotspot:
                if (f.layout == sim::ResourceFieldLayout::Cluster) {
                  ImGui::BulletText("Hotspot: local(%.2f, %.2f, %.2f)  sigma %.2f  strength %.2f",
                                    ft.localPos.x, ft.localPos.y, ft.localPos.z,
                                    ft.width, ft.strength01);
                } else {
                  ImGui::BulletText("Hotspot: center %.1f deg  sigma %.1f deg  strength %.2f", aDeg, wDeg, ft.strength01);
                }
                break;
              case sim::ResourceFieldFeatureKind::Streak:
                ImGui::BulletText("Streak: dir %.1f deg  halfWidth %.0f km  strength %.2f", aDeg, ft.width, ft.strength01);
                break;
              case sim::ResourceFieldFeatureKind::Spokes:
                ImGui::BulletText("Spokes: freq %.0f  phase %.1f deg  amp %.2f", ft.param, aDeg, ft.strength01);
                break;
            }
          }
        }

        // Yield breakdown for this field.
        std::array<int, econ::kCommodityCount> byCommodity{};
        double totalUnits = 0.0;
        double totalDensity = 0.0;
        int nAst = 0;
        for (const auto& a : asteroids) {
          if (a.fieldId != f.id) continue;
          const std::size_t idx = static_cast<std::size_t>(a.yield);
          if (idx < byCommodity.size()) {
            byCommodity[idx] += 1;
          }
          totalUnits += a.baseUnits;
          totalDensity += a.density01;
          ++nAst;
        }
        ImGui::Text("Asteroids: %d   Total baseUnits (sum): %.0f   Mean density: %.2f",
                    nAst, totalUnits, (nAst > 0) ? (totalDensity / (double)nAst) : 0.0);

        ImGui::TextDisabled("Yield mix:");
        for (std::size_t ci = 0; ci < byCommodity.size(); ++ci) {
          const int c = byCommodity[ci];
          if (c == 0) continue;
          const auto id = static_cast<econ::CommodityId>(ci);
          ImGui::BulletText("%s: %d", econ::commodityDef(id).name, c);
        }

        if (ImGui::Button("Copy selected field debug")) {
          std::string s;
          s += "ResourceField ";
          s += ui::toString(f.id);
          s += "\nKind=";
          s += sim::resourceFieldKindName(f.kind);
          s += " Layout=";
          s += sim::resourceFieldLayoutName(f.layout);
          ui::appendf(s, "\nMajorKm=%.0f", f.majorRadiusKm);
          ui::appendf(s, " MinorKm=%.0f", f.minorRadiusKm);
          ui::appendf(s, " ArcDeg=%.2f", f.arcRad * 57.29577951308232);
          s += "\nPrimary=";
          s += econ::commodityDef(f.primary).name;
          s += " Secondary=";
          s += econ::commodityDef(f.secondary).name;
          ImGui::SetClipboardText(s.c_str());
          toast("Copied field.", 1.2);
        }

        ImGui::SameLine();
        if (ImGui::Button("Copy field JSON")) {
          std::string js;
          ui::JsonWriter jw(js);
          jw.beginObject();
          jw.member("id", (core::u64)f.id);
          jw.member("kind", sim::resourceFieldKindName(f.kind));
          jw.member("layout", sim::resourceFieldLayoutName(f.layout));
          jw.member("majorRadiusKm", f.majorRadiusKm);
          jw.member("minorRadiusKm", f.minorRadiusKm);
          jw.member("arcDeg", f.arcRad * 57.29577951308232);
          jw.member("primary", econ::commodityDef(f.primary).name);
          jw.member("secondary", econ::commodityDef(f.secondary).name);

          // Export features that belong to this field.
          const auto feats = sim::filterFeaturesForField(plan.resourceFields.features, f.id);
          jw.key("features");
          jw.beginArray();
          for (const auto& ft : feats) {
            jw.beginObject();
            jw.member("kind", sim::resourceFieldFeatureKindName(ft.kind));
            jw.member("strength01", ft.strength01);
            jw.member("width", ft.width);
            jw.member("param", ft.param);
            jw.key("localPos");
            jw.beginArray();
            jw.value(ft.localPos.x);
            jw.value(ft.localPos.y);
            jw.value(ft.localPos.z);
            jw.endArray();
            jw.endObject();
          }
          jw.endArray();
          jw.endObject();

          ImGui::SetClipboardText(js.c_str());
          toast("Copied field JSON.", 1.2);
        }

        if (state.showAsteroidScatter) {
          ImGui::SeparatorText("Asteroid Scatter (field-local X/Z)");
          ImGui::SetNextItemWidth(160.0f);
          ImGui::InputInt("Max points", &state.scatterMaxPoints);
          state.scatterMaxPoints = std::clamp(state.scatterMaxPoints, 32, 4096);

          ImGui::Checkbox("Density heatmap", &state.scatterDensityHeatmap);
          ImGui::SameLine();
          ImGui::Checkbox("Color by density", &state.scatterColorByDensity);
          if (state.scatterDensityHeatmap) {
            ImGui::SameLine();
            ImGui::SetNextItemWidth(90.0f);
            ImGui::InputInt("Heatmap res", &state.scatterHeatmapRes);
            state.scatterHeatmapRes = std::clamp(state.scatterHeatmapRes, 12, 96);
          }

          const double extentKm = (f.layout == sim::ResourceFieldLayout::Torus)
                                      ? (std::max(1.0, f.majorRadiusKm + f.minorRadiusKm))
                                      : std::max(1.0, f.majorRadiusKm);

          const ImVec2 canvasSize(360.0f, 360.0f);
          ImGui::BeginChild("##rf_scatter", canvasSize, true,
                            ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);

          const ImVec2 p0 = ImGui::GetCursorScreenPos();
          const ImVec2 sz = ImGui::GetContentRegionAvail();
          const float w = std::max(32.0f, sz.x);
          const float h = std::max(32.0f, sz.y);
          const ImVec2 p1(p0.x + w, p0.y + h);

          ImDrawList* dl = ImGui::GetWindowDrawList();
          dl->AddRect(p0, p1, IM_COL32(120, 120, 120, 255));

          // Density heatmap (field-local X/Z slice at y=0).
          if (state.scatterDensityHeatmap) {
            const int N = std::clamp(state.scatterHeatmapRes, 12, 96);
            const float cellW = w / (float)N;
            const float cellH = h / (float)N;
            for (int iy = 0; iy < N; ++iy) {
              for (int ix = 0; ix < N; ++ix) {
                const double nx = (((double)ix + 0.5) / (double)N) * 2.0 - 1.0;
                const double nz = (((double)iy + 0.5) / (double)N) * 2.0 - 1.0;
                const double xKm = nx * extentKm;
                const double zKm = nz * extentKm;

                const math::Vec3d wp = f.posKm + f.basisX * xKm + f.basisZ * zKm;
                const double d01 = sim::resourceFieldDensity01(f, plan.resourceFields.features, wp);

                const int c = (int)std::clamp(30.0 + 200.0 * d01, 0.0, 255.0);
                const int a = 65;
                const ImVec2 q0(p0.x + (float)ix * cellW, p0.y + (float)iy * cellH);
                const ImVec2 q1(p0.x + (float)(ix + 1) * cellW, p0.y + (float)(iy + 1) * cellH);
                dl->AddRectFilled(q0, q1, IM_COL32(c, c, c, a));
              }
            }
          }

          const ImVec2 center((p0.x + p1.x) * 0.5f, (p0.y + p1.y) * 0.5f);
          dl->AddLine(ImVec2(center.x, p0.y), ImVec2(center.x, p1.y), IM_COL32(70, 70, 70, 255));
          dl->AddLine(ImVec2(p0.x, center.y), ImVec2(p1.x, center.y), IM_COL32(70, 70, 70, 255));

          const float scale = (std::min(w, h) * 0.5f) / (float)(extentKm * 1.05);

          // Outline helper.
          if (f.layout == sim::ResourceFieldLayout::Sheet || f.layout == sim::ResourceFieldLayout::Torus) {
            const float r0 = (float)(f.majorRadiusKm * scale);
            dl->AddCircle(center, r0, IM_COL32(90, 90, 90, 255), 96, 1.0f);
          }

          int plotted = 0;
          for (const auto& a : asteroids) {
            if (a.fieldId != f.id) continue;
            if (plotted >= state.scatterMaxPoints) break;

            const math::Vec3d d = a.posKm - f.posKm;
            const float x = (float)math::dot(d, f.basisX);
            const float z = (float)math::dot(d, f.basisZ);
            const ImVec2 pt(center.x + x * scale, center.y + z * scale);
            const double d01 = std::clamp(a.density01, 0.0, 1.0);
            const float pr = std::clamp(1.2f + (float)(a.radiusKm / 8200.0) * 2.8f, 1.2f, 4.2f);

            ImU32 col = IM_COL32(200, 200, 200, 180);
            if (state.scatterColorByDensity) {
              const int c = (int)std::clamp(80.0 + 175.0 * d01, 0.0, 255.0);
              col = IM_COL32(c, c, c, 205);
            }

            dl->AddCircleFilled(pt, pr, col);

            ++plotted;
          }

          ImGui::EndChild();
        }
      }
    }
  }

  if (ImGui::Button("Copy system debug to clipboard")) {
    std::string s;
    s += "System ";
    s += sys.stub.name;
    s += "\nId=";
    s += ui::toString(sys.stub.id);
    s += " Seed=";
    s += ui::toString(sys.stub.seed);
    s += "\nStar class=";
    s += starClassStr(sys.star.cls);
    ui::appendf(s, " M=%.2f", sys.star.massSol);
    ui::appendf(s, " R=%.2f", sys.star.radiusSol);
    ui::appendf(s, " L=%.2f", sys.star.luminositySol);
    ui::appendf(s, "\nPlanets=%zu", sys.planets.size());
    ui::appendf(s, " Moons=%zu", sys.moons.size());
    ImGui::SetClipboardText(s.c_str());
    toast("Copied.", 1.2);
  }

  ImGui::SameLine();
  if (ImGui::Button("Copy system JSON")) {
    std::string js;
    ui::JsonWriter jw(js);
    jw.beginObject();
    jw.member("id", (core::u64)sys.stub.id);
    jw.member("seed", (core::u64)sys.stub.seed);
    jw.member("name", sys.stub.name);

    jw.key("star");
    jw.beginObject();
    jw.member("class", starClassStr(sys.star.cls));
    jw.member("massSol", sys.star.massSol);
    jw.member("radiusSol", sys.star.radiusSol);
    jw.member("luminositySol", sys.star.luminositySol);
    jw.member("temperatureK", (int)sys.star.temperatureK);
    jw.endObject();

    jw.key("planets");
    jw.beginArray();
    for (const auto& p : sys.planets) {
      jw.beginObject();
      jw.member("name", p.name);
      jw.member("type", planetTypeStr(p.type));
      jw.member("radiusEarth", p.radiusEarth);
      jw.member("massEarth", p.massEarth);
      jw.key("orbit");
      jw.beginObject();
      jw.member("semiMajorAxisAU", p.orbit.semiMajorAxisAU);
      jw.member("eccentricity", p.orbit.eccentricity);
      jw.member("inclinationRad", p.orbit.inclinationRad);
      jw.member("periodDays", p.orbit.periodDays);
      jw.endObject();
      jw.endObject();
    }
    jw.endArray();

    jw.key("moons");
    jw.beginArray();
    for (const auto& m : sys.moons) {
      jw.beginObject();
      jw.member("id", (core::u64)m.id);
      jw.member("name", m.name);
      jw.member("type", planetTypeStr(m.type));
      jw.member("radiusEarth", m.radiusEarth);
      jw.member("massEarth", m.massEarth);
      jw.member("parentPlanetIndex", (core::u32)m.parentPlanetIndex);
      jw.key("orbit");
      jw.beginObject();
      jw.member("semiMajorAxisAU", m.orbit.semiMajorAxisAU);
      jw.member("eccentricity", m.orbit.eccentricity);
      jw.member("inclinationRad", m.orbit.inclinationRad);
      jw.member("periodDays", m.orbit.periodDays);
      jw.endObject();
      jw.endObject();
    }
    jw.endArray();

    jw.endObject();

    ImGui::SetClipboardText(js.c_str());
    toast("Copied system JSON.", 1.2);
  }

  ImGui::End();
}

} // namespace stellar::game
