#include "ProceduralGalaxyLabWindow.h"

#include "stellar/core/Hash.h"
#include "stellar/core/Random.h"

#include "stellar/math/Math.h"

#include <imgui.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <unordered_map>

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

static ImU32 colorForRegionKind(proc::GalaxyRegionKind k) {
  // Intentionally stylized palette (not physically based).
  switch (k) {
    case proc::GalaxyRegionKind::Core: return rgba(1.00f, 0.86f, 0.40f);
    case proc::GalaxyRegionKind::InnerDisc: return rgba(0.55f, 0.90f, 0.55f);
    case proc::GalaxyRegionKind::OuterRim: return rgba(0.55f, 0.70f, 1.00f);
    case proc::GalaxyRegionKind::Nebula: return rgba(0.95f, 0.55f, 0.95f);
    case proc::GalaxyRegionKind::Cluster: return rgba(0.85f, 0.85f, 0.85f);
    case proc::GalaxyRegionKind::Rift: return rgba(0.90f, 0.30f, 0.30f);
    default: return rgba(0.75f, 0.75f, 0.75f);
  }
}

static ImU32 colorForHyperlane(const proc::HyperlaneEdge& e) {
  const float risk = static_cast<float>(std::clamp(e.risk01, 0.0, 1.0));
  const float bw = static_cast<float>(std::clamp(e.bandwidth01, 0.0, 1.0));

  // Low risk: bluish, high risk: reddish; bandwidth controls alpha.
  const float r = (1.0f - risk) * 0.30f + risk * 0.95f;
  const float g = (1.0f - risk) * 0.65f + risk * 0.35f;
  const float b = (1.0f - risk) * 0.95f + risk * 0.30f;
  const float a = 0.10f + 0.55f * bw;
  return rgba(r, g, b, a);
}

static bool dragDouble(const char* label, double& v, double step, double minV, double maxV, const char* fmt = "%.3f") {
  double tmp = v;
  const bool changed = ImGui::DragScalar(label, ImGuiDataType_Double, &tmp, step, &minV, &maxV, fmt);
  if (changed) v = tmp;
  return changed;
}

static void rebuildHyperlanes(ProceduralGalaxyLabWindowState& st);

static void rebuildPreview(ProceduralGalaxyLabWindowState& st) {
  const auto t0 = Clock::now();

  const int fc = std::max(1, st.factionCount);
  st.factions = sim::generateFactions(st.seed, static_cast<core::u32>(fc));

  proc::GalaxyGenerator gen(st.seed, st.params);

  st.stubs.clear();
  st.stubs.reserve(4096);

  st.stubRegionKind.clear();
  st.stubRegionId.clear();
  st.stubRegionKind.reserve(4096);
  st.stubRegionId.reserve(4096);

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

          // Cache region data for optional visualization.
          {
            const auto reg = proc::sampleGalaxyRegion(st.seed, stub.posLy, st.regionCellSizeLy);
            st.stubRegionKind.push_back(reg.kind);
            st.stubRegionId.push_back(reg.regionId);
          }

          if (st.stubs.size() >= st.maxStubs) goto done;
        }
      }
    }
  }

done:
  st.lastGenMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
  rebuildHyperlanes(st);
  st.dirty = false;
}

static void rebuildHyperlanes(ProceduralGalaxyLabWindowState& st) {
  st.hyperlanes.clear();
  st.lastLaneMs = 0.0;

  if (!st.showHyperlanes) {
    return;
  }
  if (st.stubs.size() < 2) {
    return;
  }

  const auto t0 = Clock::now();

  // NOTE: Hyperlane generation is deterministic and depends only on
  // (seed, stub list, params). It can be recomputed without regenerating stubs.
  proc::HyperlaneParams p = st.hyperlaneParams;

  // Keep lane region sampling resolution reasonably aligned with the region
  // visualization cell size unless the user explicitly disables it.
  if (p.regionCellSizeLy > 0.0) {
    p.regionCellSizeLy = std::max(1.0, p.regionCellSizeLy);
  }

  st.hyperlanes = proc::generateHyperlaneNetwork(st.seed, st.stubs, p).edges;
  st.lastLaneMs = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

void drawProceduralGalaxyLabWindow(ProceduralGalaxyLabWindowState& st, float timeSec, const ToastFn& toast) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(1180, 760), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Procedural Galaxy Lab", &st.open)) {
    ImGui::End();
    return;
  }

  bool dirty = false;
  bool laneDirty = false;

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
      dirty |= ImGui::Checkbox("Color by region", &st.colorByRegion);
      ImGui::SameLine();
      dirty |= ImGui::Checkbox("Legend", &st.showLegend);

      laneDirty |= ImGui::Checkbox("Hyperlanes", &st.showHyperlanes);

      if (st.colorByRegion) {
        ImGui::Indent();
        const bool regChanged = dragDouble("Region Cell Size (ly)", st.regionCellSizeLy, 10.0, 100.0, 5000.0, "%.0f");
        dirty |= regChanged;
        if (regChanged && st.hyperlaneParams.regionCellSizeLy > 0.0) {
          st.hyperlaneParams.regionCellSizeLy = st.regionCellSizeLy;
          laneDirty = true;
        }
        ImGui::Unindent();
      }

      if (st.showHyperlanes) {
        ImGui::Indent();
        laneDirty |= dragDouble("Lane Max Dist (ly)", st.hyperlaneParams.maxNeighborDistanceLy, 0.5, 2.0, 200.0, "%.1f");
        laneDirty |= ImGui::SliderInt("Lane Neighbor K", &st.hyperlaneParams.neighborK, 1, 12);
        laneDirty |= ImGui::Checkbox("Lane Force Connected", &st.hyperlaneParams.forceConnected);
        laneDirty |= ImGui::SliderInt("Lane Min Degree", &st.hyperlaneParams.minDegree, 0, 6);
        laneDirty |= dragDouble("Lane Extra Edge Chance", st.hyperlaneParams.extraEdgeChance, 0.01, 0.0, 1.0, "%.2f");
        laneDirty |= dragDouble("Lane Region Cell (ly)", st.hyperlaneParams.regionCellSizeLy, 10.0, 0.0, 5000.0, "%.0f");
        ImGui::Unindent();
      }

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
      if (st.showHyperlanes) {
        ImGui::Text("Hyperlanes: %zu", st.hyperlanes.size());
        ImGui::Text("Last lane gen: %.2f ms", st.lastLaneMs);
      }
      if (st.stubs.size() >= st.maxStubs) {
        ImGui::TextDisabled("(hit maxStubs cap: %zu)", st.maxStubs);
      }
    }

    ImGui::TableNextColumn();

    // ----- Map preview -----
    if (dirty) st.dirty = true;
    if (st.dirty && st.autoRegenerate) {
      rebuildPreview(st);
    } else if (laneDirty) {
      rebuildHyperlanes(st);
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

    // Hyperlane overlay (drawn behind points).
    if (st.showHyperlanes && !st.hyperlanes.empty()) {
      std::unordered_map<sim::SystemId, math::Vec3d> idToPos;
      idToPos.reserve(st.stubs.size() * 2);
      for (const auto& stub : st.stubs) {
        idToPos.emplace(stub.id, stub.posLy);
      }

      for (const auto& e : st.hyperlanes) {
        const auto itA = idToPos.find(e.a);
        const auto itB = idToPos.find(e.b);
        if (itA == idToPos.end() || itB == idToPos.end()) continue;

        const ImVec2 a = worldToScreen(itA->second);
        const ImVec2 b = worldToScreen(itB->second);

        // Quick reject: if both endpoints are on the same outside side, skip.
        if ((a.x < p0.x && b.x < p0.x) || (a.x > p1.x && b.x > p1.x) ||
            (a.y < p0.y && b.y < p0.y) || (a.y > p1.y && b.y > p1.y)) {
          continue;
        }

        const float bw = static_cast<float>(std::clamp(e.bandwidth01, 0.0, 1.0));
        const float thickness = 1.0f + 1.6f * bw;
        dl->AddLine(a, b, colorForHyperlane(e), thickness);
      }
    }


    // Draw points.
    for (std::size_t i = 0; i < st.stubs.size(); ++i) {
      const auto& stub = st.stubs[i];

      const ImVec2 p = worldToScreen(stub.posLy);
      if (p.x < p0.x || p.x > p1.x || p.y < p0.y || p.y > p1.y) continue;

      ImU32 col = colorForStarClass(stub.primaryClass);
      if (st.colorByFaction) {
        col = colorForFaction(stub.factionId, st.factions);
      } else if (st.colorByRegion && i < st.stubRegionKind.size()) {
        col = colorForRegionKind(st.stubRegionKind[i]);
      }
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


        if (st.showHyperlanes && !st.hyperlanes.empty()) {
          int degree = 0;
          double bestBw = 0.0;
          double bestRisk = 0.0;
          for (const auto& e : st.hyperlanes) {
            if (e.a == stub.id || e.b == stub.id) {
              degree += 1;
              bestBw = std::max(bestBw, e.bandwidth01);
              bestRisk = std::max(bestRisk, e.risk01);
            }
          }

          ImGui::Separator();
          ImGui::Text("Hyperlanes: %d", degree);
          ImGui::Text("Max lane BW: %.2f", bestBw);
          ImGui::Text("Max lane Risk: %.2f", bestRisk);
        }

        if (st.colorByRegion) {
          const auto reg = proc::sampleGalaxyRegion(st.seed, stub.posLy, st.regionCellSizeLy);
          ImGui::Separator();
          ImGui::Text("Region: %s", reg.name.c_str());
          ImGui::Text("Kind: %s", proc::galaxyRegionKindName(reg.kind));
          ImGui::Text("Edge: %.2f", reg.edge01);
        }
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

      const bool legendFaction = st.colorByFaction;
      const bool legendRegion = (!legendFaction && st.colorByRegion);
      ImGui::TextUnformatted(legendFaction ? "Legend (Faction)" : (legendRegion ? "Legend (Region)" : "Legend (Star Class)"));
      ImGui::Separator();

      if (legendRegion) {
        const proc::GalaxyRegionKind kinds[] = {
            proc::GalaxyRegionKind::Core,
            proc::GalaxyRegionKind::InnerDisc,
            proc::GalaxyRegionKind::OuterRim,
            proc::GalaxyRegionKind::Nebula,
            proc::GalaxyRegionKind::Cluster,
            proc::GalaxyRegionKind::Rift,
        };
        for (auto k : kinds) {
          ImGui::ColorButton("##r", ImGui::ColorConvertU32ToFloat4(colorForRegionKind(k)), ImGuiColorEditFlags_NoTooltip, ImVec2(16, 16));
          ImGui::SameLine();
          ImGui::TextUnformatted(proc::galaxyRegionKindName(k));
        }
      } else if (!legendFaction) {
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
