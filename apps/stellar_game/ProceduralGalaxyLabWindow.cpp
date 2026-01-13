#include "ProceduralGalaxyLabWindow.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include "stellar/math/Math.h"

#include <imgui.h>

#include <algorithm>
#include <chrono>
#include <cmath>

namespace stellar::game {

using Clock = std::chrono::high_resolution_clock;

static ImU32 rgba(float r, float g, float b, float a = 1.0f) {
  return ImGui::GetColorU32(ImVec4(r, g, b, a));
}

static ImU32 colorForStarClass(sim::StarClass c) {
  // Rough temperature-ish palette; intentionally stylized.
  switch (c) {
    case sim::StarClass::O: return rgba(0.55f, 0.70f, 1.00f);
    case sim::StarClass::B: return rgba(0.65f, 0.75f, 1.00f);
    case sim::StarClass::A: return rgba(0.80f, 0.85f, 1.00f);
    case sim::StarClass::F: return rgba(1.00f, 0.95f, 0.80f);
    case sim::StarClass::G: return rgba(1.00f, 0.90f, 0.60f);
    case sim::StarClass::K: return rgba(1.00f, 0.75f, 0.45f);
    case sim::StarClass::M: return rgba(1.00f, 0.55f, 0.40f);
    default: return rgba(1.0f, 1.0f, 1.0f);
  }
}

static ImU32 colorForFaction(core::u32 id, const std::vector<sim::Faction>& factions) {
  if (id == 0) return rgba(0.75f, 0.75f, 0.75f);
  for (const auto& f : factions) {
    if (f.id == id) {
      return rgba(static_cast<float>(f.colorRgb.x),
                  static_cast<float>(f.colorRgb.y),
                  static_cast<float>(f.colorRgb.z));
    }
  }
  return rgba(0.75f, 0.75f, 0.75f);
}

static bool dragDouble(const char* label, double& v, double step, double minV, double maxV, const char* fmt = "%.3f") {
  double tmp = v;
  const bool changed = ImGui::DragScalar(label, ImGuiDataType_Double, &tmp, step, &minV, &maxV, fmt);
  if (changed) v = tmp;
  return changed;
}

static void rebuildPreview(ProceduralGalaxyLabWindowState& st) {
  const auto t0 = Clock::now();

  const int fc = std::max(1, st.factionCount);
  st.factions = sim::generateFactions(st.seed, static_cast<core::u32>(fc));

  proc::GalaxyGenerator gen(st.seed, st.params);

  st.stubs.clear();
  st.stubs.reserve(4096);

  const double s = std::max(1.0, st.params.sectorSizeLy);
  const double r = std::max(1.0, st.viewRadiusLy);
  const double zHalf = std::max(0.0, st.zHalfLy);

  const math::Vec3d c = st.centerLy;

  const int minX = static_cast<int>(std::floor((c.x - r) / s));
  const int maxX = static_cast<int>(std::floor((c.x + r) / s));
  const int minY = static_cast<int>(std::floor((c.y - r) / s));
  const int maxY = static_cast<int>(std::floor((c.y + r) / s));
  const int minZ = static_cast<int>(std::floor((c.z - zHalf) / s));
  const int maxZ = static_cast<int>(std::floor((c.z + zHalf) / s));

  const double r2 = r * r;

  for (int z = minZ; z <= maxZ; ++z) {
    for (int y = minY; y <= maxY; ++y) {
      for (int x = minX; x <= maxX; ++x) {
        proc::Sector sec = gen.generateSector({x, y, z}, st.factions);
        for (const auto& stub : sec.systems) {
          if (std::abs(stub.posLy.z - c.z) > zHalf) continue;
          const double dx = stub.posLy.x - c.x;
          const double dy = stub.posLy.y - c.y;
          if (dx * dx + dy * dy > r2) continue;

          st.stubs.push_back(stub);
          if (st.stubs.size() >= st.maxStubs) goto done;
        }
      }
    }
  }

done:
  st.lastGenMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
  st.dirty = false;
}

void drawProceduralGalaxyLabWindow(ProceduralGalaxyLabWindowState& st, float timeSec, const ToastFn& toast) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(1180, 760), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Procedural Galaxy Lab", &st.open)) {
    ImGui::End();
    return;
  }

  bool dirty = false;

  if (ImGui::BeginTable("##galaxy_lab_layout", 2, ImGuiTableFlags_Resizable | ImGuiTableFlags_BordersInnerV)) {
    ImGui::TableNextColumn();

    // ----- Controls -----
    {
      ImGui::TextUnformatted("Preview generator params in a top-down XY slice.");
      ImGui::TextDisabled("Tip: set Spiral Arms > 0 and Arm Strength > 0 to enable.");
      ImGui::Separator();

      dirty |= ImGui::InputScalar("Seed", ImGuiDataType_U64, &st.seed);

      ImGui::SameLine();
      if (ImGui::Button("Randomize")) {
        const core::u64 t = static_cast<core::u64>(timeSec * 1'000'000.0f);
        const core::u64 h = core::hashCombine(st.seed, core::hashCombine(core::fnv1a64("galaxy_lab"), t));
        core::SplitMix64 rng(h);
        st.seed = rng.nextU64();
        dirty = true;
        if (toast) toast("Galaxy seed randomized", 2.0);
      }

      dirty |= ImGui::SliderInt("Factions", &st.factionCount, 1, 24);

      ImGui::Separator();
      ImGui::TextUnformatted("Galaxy Base");

      dirty |= dragDouble("Sector Size (ly)", st.params.sectorSizeLy, 0.1, 1.0, 200.0, "%.1f");
      dirty |= dragDouble("Disc Radius (ly)", st.params.radiusLy, 50.0, 500.0, 250'000.0, "%.0f");
      dirty |= dragDouble("Disc Thickness (ly)", st.params.thicknessLy, 10.0, 10.0, 25'000.0, "%.0f");
      dirty |= dragDouble("Radial Scale (ly)", st.params.radialScaleLengthLy, 50.0, 100.0, 250'000.0, "%.0f");
      dirty |= dragDouble("Mean Systems/Sector", st.params.baseMeanSystemsPerSector, 0.05, 0.0, 100.0, "%.2f");

      ImGui::Separator();
      ImGui::TextUnformatted("Spiral Arms (log spiral)");

      dirty |= ImGui::SliderInt("Arm Count", &st.params.spiralArmCount, 0, 8);
      dirty |= dragDouble("Arm Strength", st.params.spiralArmStrength, 0.02, 0.0, 5.0, "%.2f");
      dirty |= dragDouble("Pitch (deg)", st.params.spiralPitchDeg, 0.1, 1.0, 45.0, "%.1f");
      dirty |= dragDouble("Arm Width (deg)", st.params.spiralArmWidthDeg, 0.1, 1.0, 90.0, "%.1f");
      dirty |= dragDouble("Arm Phase (deg)", st.params.spiralArmPhaseDeg, 0.25, -360.0, 360.0, "%.2f");
      dirty |= dragDouble("Arm Noise Strength", st.params.spiralArmNoiseStrength, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Arm Noise Freq", st.params.spiralArmNoiseFreq, 0.0001, 0.0, 0.05, "%.4f");

      ImGui::Separator();
      ImGui::TextUnformatted("Density Noise (clumpiness)");

      dirty |= dragDouble("Density Noise Strength", st.params.densityNoiseStrength, 0.01, 0.0, 1.0, "%.2f");
      dirty |= dragDouble("Density Noise Freq", st.params.densityNoiseFreq, 0.0001, 0.0, 0.05, "%.4f");

      ImGui::Separator();
      ImGui::TextUnformatted("View");

      dirty |= dragDouble("Center X (ly)", st.centerLy.x, 5.0, -500'000.0, 500'000.0, "%.0f");
      dirty |= dragDouble("Center Y (ly)", st.centerLy.y, 5.0, -500'000.0, 500'000.0, "%.0f");
      dirty |= dragDouble("Center Z (ly)", st.centerLy.z, 5.0, -500'000.0, 500'000.0, "%.0f");
      dirty |= dragDouble("View Radius (ly)", st.viewRadiusLy, 5.0, 10.0, 50'000.0, "%.0f");
      dirty |= dragDouble("Z Half-Range (ly)", st.zHalfLy, 5.0, 0.0, 25'000.0, "%.0f");

      dirty |= ImGui::Checkbox("Auto regenerate", &st.autoRegenerate);
      dirty |= ImGui::Checkbox("Color by faction", &st.colorByFaction);
      ImGui::SameLine();
      dirty |= ImGui::Checkbox("Legend", &st.showLegend);

      ImGui::Checkbox("Arm guides", &st.showArmGuides);

      if (ImGui::Button("Regenerate")) {
        st.dirty = true;
      }

      ImGui::SameLine();
      if (ImGui::Button("Reset Params")) {
        st.params = proc::GalaxyParams{};
        st.dirty = true;
        if (toast) toast("Galaxy params reset", 2.0);
      }

      ImGui::Separator();
      ImGui::Text("Preview: %zu systems", st.stubs.size());
      ImGui::Text("Last gen: %.2f ms", st.lastGenMs);
      if (st.stubs.size() >= st.maxStubs) {
        ImGui::TextDisabled("(hit maxStubs cap: %zu)", st.maxStubs);
      }
    }

    ImGui::TableNextColumn();

    // ----- Map preview -----
    if (dirty) st.dirty = true;
    if (st.dirty && st.autoRegenerate) {
      rebuildPreview(st);
    }

    const ImVec2 canvasSize = ImGui::GetContentRegionAvail();
    ImGui::BeginChild("##galaxy_canvas", canvasSize, true, ImGuiWindowFlags_NoScrollWithMouse | ImGuiWindowFlags_NoScrollbar);

    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const ImVec2 sz = ImGui::GetContentRegionAvail();
    const ImVec2 p1 = ImVec2(p0.x + sz.x, p0.y + sz.y);

    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->AddRectFilled(p0, p1, rgba(0.06f, 0.07f, 0.10f));
    dl->AddRect(p0, p1, rgba(0.30f, 0.30f, 0.35f));

    const ImVec2 centerPx = ImVec2((p0.x + p1.x) * 0.5f, (p0.y + p1.y) * 0.5f);

    const double radiusLy = std::max(1.0, st.viewRadiusLy);
    const float scale = static_cast<float>(std::min(sz.x, sz.y) / (2.0 * radiusLy));

    const auto worldToScreen = [&](const math::Vec3d& w) -> ImVec2 {
      const double dx = w.x - st.centerLy.x;
      const double dy = w.y - st.centerLy.y;
      return ImVec2(centerPx.x + static_cast<float>(dx * scale),
                   centerPx.y - static_cast<float>(dy * scale));
    };

    // Axes + boundary.
    dl->AddLine(ImVec2(p0.x, centerPx.y), ImVec2(p1.x, centerPx.y), rgba(0.22f, 0.22f, 0.25f));
    dl->AddLine(ImVec2(centerPx.x, p0.y), ImVec2(centerPx.x, p1.y), rgba(0.22f, 0.22f, 0.25f));
    dl->AddCircle(centerPx, static_cast<float>(radiusLy * scale), rgba(0.25f, 0.25f, 0.28f));

    // Optional arm guides (only meaningful if you're looking near the galactic center).
    if (st.showArmGuides && st.params.spiralArmCount > 0 && st.params.spiralArmStrength > 0.0) {
      const int arms = st.params.spiralArmCount;
      const double pitchDeg = std::clamp(st.params.spiralPitchDeg, 1.0, 89.0);
      const double pitchRad = pitchDeg * (stellar::math::kPi / 180.0);
      const double k = 1.0 / std::tan(pitchRad);
      const double phaseRad = st.params.spiralArmPhaseDeg * (stellar::math::kPi / 180.0);
      const double twoPi = 2.0 * stellar::math::kPi;
      const double rRef = std::max(1.0, st.params.radiusLy * 0.02);

      // Draw each arm as a polyline. We sample r in log space.
      const int steps = 220;
      for (int a = 0; a < arms; ++a) {
        ImVec2 prev{};
        bool havePrev = false;
        for (int i = 0; i < steps; ++i) {
          const double t = static_cast<double>(i) / static_cast<double>(steps - 1);
          const double r = rRef * std::exp(std::log(std::max(1.0, radiusLy / rRef)) * t);
          const double lnTerm = std::log(std::max(1.0e-6, r / rRef));
          const double theta = k * lnTerm + phaseRad + twoPi * (static_cast<double>(a) / static_cast<double>(arms));

          const math::Vec3d w{r * std::cos(theta), r * std::sin(theta), st.centerLy.z};
          const ImVec2 p = worldToScreen(w);
          if (havePrev) {
            dl->AddLine(prev, p, rgba(0.16f, 0.19f, 0.28f, 0.65f));
          }
          prev = p;
          havePrev = true;
        }
      }
    }

    // Draw points.
    for (const auto& stub : st.stubs) {
      const ImVec2 p = worldToScreen(stub.posLy);
      if (p.x < p0.x || p.x > p1.x || p.y < p0.y || p.y > p1.y) continue;

      const ImU32 col = st.colorByFaction ? colorForFaction(stub.factionId, st.factions)
                                          : colorForStarClass(stub.primaryClass);
      dl->AddCircleFilled(p, 2.0f, col);
    }

    // Hover tooltip.
    if (ImGui::IsMouseHoveringRect(p0, p1) && !st.stubs.empty()) {
      const ImVec2 m = ImGui::GetIO().MousePos;
      int bestIdx = -1;
      float bestD2 = 6.0f * 6.0f;

      for (int i = 0; i < static_cast<int>(st.stubs.size()); ++i) {
        const ImVec2 p = worldToScreen(st.stubs[static_cast<std::size_t>(i)].posLy);
        const float dx = m.x - p.x;
        const float dy = m.y - p.y;
        const float d2 = dx * dx + dy * dy;
        if (d2 < bestD2) {
          bestD2 = d2;
          bestIdx = i;
        }
      }

      if (bestIdx >= 0) {
        const auto& stub = st.stubs[static_cast<std::size_t>(bestIdx)];
        ImGui::BeginTooltip();
        ImGui::TextUnformatted(stub.name.c_str());
        ImGui::Separator();
        ImGui::Text("SystemId: 0x%llx", static_cast<unsigned long long>(stub.id));
        ImGui::Text("Pos (ly): [%.1f, %.1f, %.1f]", stub.posLy.x, stub.posLy.y, stub.posLy.z);
        ImGui::Text("Star: %c", static_cast<char>(stub.primaryClass));
        ImGui::Text("Planets: %u", stub.planetCount);
        ImGui::Text("Stations: %u", stub.stationCount);
        ImGui::Text("Faction: %u", stub.factionId);
        ImGui::TextDisabled("Click to copy name");
        ImGui::EndTooltip();

        if (ImGui::IsMouseClicked(0)) {
          ImGui::SetClipboardText(stub.name.c_str());
          if (toast) toast("Copied system name to clipboard", 1.5);
        }
      }
    }

    // Legend.
    if (st.showLegend) {
      ImGui::SetCursorScreenPos(ImVec2(p0.x + 12.0f, p0.y + 12.0f));
      ImGui::BeginChild("##legend", ImVec2(210, 160), true);
      ImGui::TextUnformatted(st.colorByFaction ? "Legend (Faction)" : "Legend (Star Class)");
      ImGui::Separator();

      if (!st.colorByFaction) {
        const sim::StarClass classes[] = {sim::StarClass::O, sim::StarClass::B, sim::StarClass::A, sim::StarClass::F, sim::StarClass::G, sim::StarClass::K, sim::StarClass::M};
        for (auto c : classes) {
          ImGui::ColorButton("##c", ImGui::ColorConvertU32ToFloat4(colorForStarClass(c)), ImGuiColorEditFlags_NoTooltip, ImVec2(16, 16));
          ImGui::SameLine();
          ImGui::Text("%c", static_cast<char>(c));
        }
      } else {
        for (const auto& f : st.factions) {
          ImGui::ColorButton("##f", ImVec4(static_cast<float>(f.colorRgb.x), static_cast<float>(f.colorRgb.y), static_cast<float>(f.colorRgb.z), 1.0f), ImGuiColorEditFlags_NoTooltip, ImVec2(16, 16));
          ImGui::SameLine();
          ImGui::Text("%u %s", f.id, f.name.c_str());
        }
      }

      ImGui::EndChild();
    }

    ImGui::EndChild();

    ImGui::EndTable();
  }

  ImGui::End();
}

} // namespace stellar::game
