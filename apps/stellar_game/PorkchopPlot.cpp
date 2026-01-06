#include "PorkchopPlot.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

namespace stellar::game {
namespace {

static double metricValue(const stellar::sim::LambertPorkchopCell& c, PorkchopMetric m) {
  switch (m) {
    case PorkchopMetric::Score: return c.score;
    case PorkchopMetric::DepartDv: return c.dvDepartMagKmS;
    case PorkchopMetric::ArrivalRelSpeed: return c.arriveRelSpeedKmS;
    case PorkchopMetric::TotalDv: return c.totalDvKmS;
    default: return c.score;
  }
}

// A lightweight color ramp: best (low) -> green, worst (high) -> red.
static ImU32 heatColor(float t01, float alpha = 1.0f) {
  t01 = std::clamp(t01, 0.0f, 1.0f);
  // Hue 0.33 (~green) down to 0.0 (red).
  const float hue = (1.0f - t01) * 0.33f;
  ImVec4 c = ImColor::HSV(hue, 0.95f, 0.95f);
  c.w = alpha;
  return ImGui::GetColorU32(c);
}

struct PlotCell {
  bool ok{false};
  double departAfterSec{0.0};
  double tofSec{0.0};

  double dvDepartMagKmS{0.0};
  double arriveRelSpeedKmS{0.0};
  double totalDvKmS{0.0};
  double score{0.0};
};

} // namespace

PorkchopPlotPick drawPorkchopPlot(const sim::LambertPorkchopResult& res,
                                  const std::vector<sim::LambertTransferMetrics>& best,
                                  const PorkchopPlotOptions& opt,
                                  ImVec2 size) {
  PorkchopPlotPick pick{};

  const int depSteps = res.departSteps;
  const int tofSteps = res.tofSteps;
  const std::size_t expected = (depSteps > 0 && tofSteps > 0) ? (std::size_t)depSteps * (std::size_t)tofSteps : 0u;

  if (depSteps <= 0 || tofSteps <= 0) {
    ImGui::TextDisabled("No porkchop grid.");
    return pick;
  }
  if (res.grid.size() < expected) {
    ImGui::TextDisabled("Porkchop grid not complete (computed %d / %d cells).",
                        (int)res.grid.size(), (int)expected);
    ImGui::TextDisabled("Tip: wait for the search to finish, or increase Cells/frame.");
    return pick;
  }

  // Resolve plot size.
  if (size.x <= 0.0f) size.x = ImGui::GetContentRegionAvail().x;
  if (size.y <= 0.0f) size.y = 260.0f;
  if (size.x < 80.0f) size.x = 80.0f;
  if (size.y < 80.0f) size.y = 80.0f;

  // Downsample bins.
  const int cols = std::max(1, std::min(depSteps, std::max(1, opt.maxCols)));
  const int rows = std::max(1, std::min(tofSteps, std::max(1, opt.maxRows)));
  const int skipDep = std::max(1, (depSteps + cols - 1) / cols);
  const int skipTof = std::max(1, (tofSteps + rows - 1) / rows);
  const int pDep = (depSteps + skipDep - 1) / skipDep;
  const int pTof = (tofSteps + skipTof - 1) / skipTof;

  std::vector<PlotCell> binned((std::size_t)pDep * (std::size_t)pTof);
  std::vector<double> binnedMetric((std::size_t)pDep * (std::size_t)pTof, DBL_MAX);

  // Bin by keeping the best (min) metric in each block.
  for (int d = 0; d < depSteps; ++d) {
    const int bd = d / skipDep;
    for (int t = 0; t < tofSteps; ++t) {
      const int bt = t / skipTof;
      const std::size_t bi = (std::size_t)bd * (std::size_t)pTof + (std::size_t)bt;
      const std::size_t gi = (std::size_t)d * (std::size_t)tofSteps + (std::size_t)t;
      const auto& c = res.grid[gi];
      if (!c.ok) continue;
      const double mv = metricValue(c, opt.metric);
      if (!(mv < binnedMetric[bi])) continue;

      binnedMetric[bi] = mv;
      PlotCell pc;
      pc.ok = true;
      pc.departAfterSec = c.departAfterSec;
      pc.tofSec = c.tofSec;
      pc.dvDepartMagKmS = c.dvDepartMagKmS;
      pc.arriveRelSpeedKmS = c.arriveRelSpeedKmS;
      pc.totalDvKmS = c.totalDvKmS;
      pc.score = c.score;
      binned[bi] = pc;
    }
  }

  // Range for color mapping.
  double vMin = DBL_MAX;
  double vMax = -DBL_MAX;
  for (std::size_t i = 0; i < binned.size(); ++i) {
    if (!binned[i].ok) continue;
    const double mv = binnedMetric[i];
    if (!std::isfinite(mv)) continue;
    vMin = std::min(vMin, mv);
    vMax = std::max(vMax, mv);
  }
  if (!(vMin < vMax)) {
    // Degenerate: either only one value or no ok cells.
    vMin = 0.0;
    vMax = 1.0;
  }

  const ImVec2 p0 = ImGui::GetCursorScreenPos();
  const ImVec2 p1 = ImVec2(p0.x + size.x, p0.y + size.y);
  const float cellW = size.x / (float)pDep;
  const float cellH = size.y / (float)pTof;

  ImDrawList* dl = ImGui::GetWindowDrawList();
  dl->AddRectFilled(p0, p1, ImGui::GetColorU32(ImGuiCol_FrameBg));

  // Heatmap cells.
  for (int bd = 0; bd < pDep; ++bd) {
    for (int bt = 0; bt < pTof; ++bt) {
      const std::size_t bi = (std::size_t)bd * (std::size_t)pTof + (std::size_t)bt;
      const PlotCell& c = binned[bi];

      const ImVec2 a = ImVec2(p0.x + (float)bd * cellW, p0.y + (float)bt * cellH);
      const ImVec2 b = ImVec2(a.x + cellW, a.y + cellH);

      if (!c.ok) {
        // Unreachable / invalid area.
        dl->AddRectFilled(a, b, ImGui::GetColorU32(ImVec4(0, 0, 0, 0.22f)));
        continue;
      }

      const double mv = binnedMetric[bi];
      float t01 = 0.0f;
      if (std::isfinite(mv) && vMax > vMin) {
        t01 = (float)((mv - vMin) / (vMax - vMin));
      }
      t01 = std::clamp(t01, 0.0f, 1.0f);
      if (opt.gamma > 0.01f && opt.gamma != 1.0f) {
        t01 = std::pow(t01, opt.gamma);
      }

      dl->AddRectFilled(a, b, heatColor(t01, 0.92f));
    }
  }

  // Border.
  dl->AddRect(p0, p1, ImGui::GetColorU32(ImGuiCol_Border));

  // Mark best candidates.
  if (opt.showBestMarkers && !best.empty()) {
    // Infer axis extents from the full grid (uniform spacing, row-major).
    const double depMin = res.grid[0].departAfterSec;
    const double depMax = res.grid[(std::size_t)(depSteps - 1) * (std::size_t)tofSteps].departAfterSec;
    const double tofMin = res.grid[0].tofSec;
    const double tofMax = res.grid[(std::size_t)tofSteps - 1].tofSec;

    const double depDen = (depMax - depMin);
    const double tofDen = (tofMax - tofMin);

    for (std::size_t i = 0; i < best.size(); ++i) {
      const auto& c = best[i];
      if (!c.ok) continue;
      const double u = (depDen != 0.0) ? (c.departAfterSec - depMin) / depDen : 0.0;
      const double v = (tofDen != 0.0) ? (c.tofSec - tofMin) / tofDen : 0.0;

      const float x = p0.x + (float)std::clamp(u, 0.0, 1.0) * size.x;
      const float y = p0.y + (float)std::clamp(v, 0.0, 1.0) * size.y;

      const float r = 3.2f;
      const ImU32 col = ImGui::GetColorU32(ImVec4(1, 1, 1, 0.95f));
      dl->AddCircle(ImVec2(x, y), r, col, 12, 1.2f);
    }
  }

  // Hover + click handling.
  ImGui::InvisibleButton("##porkchop_plot", size);
  const bool hovered = ImGui::IsItemHovered();
  if (hovered && cellW > 0.0f && cellH > 0.0f) {
    const ImVec2 mp = ImGui::GetIO().MousePos;
    const int bd = std::clamp((int)((mp.x - p0.x) / cellW), 0, pDep - 1);
    const int bt = std::clamp((int)((mp.y - p0.y) / cellH), 0, pTof - 1);
    const std::size_t bi = (std::size_t)bd * (std::size_t)pTof + (std::size_t)bt;

    const PlotCell& c = binned[bi];
    if (c.ok) {
      ImGui::BeginTooltip();
      ImGui::Text("Depart: +%.2f h", c.departAfterSec / 3600.0);
      ImGui::Text("TOF: %.2f h", c.tofSec / 3600.0);
      ImGui::Separator();
      ImGui::Text("Δv dep: %.1f m/s", c.dvDepartMagKmS * 1000.0);
      ImGui::Text("Arr rel: %.1f m/s", c.arriveRelSpeedKmS * 1000.0);
      ImGui::Text("Total Δv: %.1f m/s", c.totalDvKmS * 1000.0);
      ImGui::Text("Score: %.5f", c.score);
      ImGui::Separator();
      ImGui::TextDisabled("Click to apply this cell.");
      ImGui::EndTooltip();
    }

    if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
      if (c.ok) {
        pick.clicked = true;
        pick.departAfterSec = c.departAfterSec;
        pick.tofSec = c.tofSec;
      }
    }
  }

  // Axis labels.
  if (opt.showAxesLabels) {
    const double depMin = res.grid[0].departAfterSec / 3600.0;
    const double depMax = res.grid[(std::size_t)(depSteps - 1) * (std::size_t)tofSteps].departAfterSec / 3600.0;
    const double tofMin = res.grid[0].tofSec / 3600.0;
    const double tofMax = res.grid[(std::size_t)tofSteps - 1].tofSec / 3600.0;

    ImGui::TextDisabled("Departure window: +%.2f h … +%.2f h", depMin, depMax);
    ImGui::TextDisabled("TOF window: %.2f h … %.2f h", tofMin, tofMax);
  }

  return pick;
}

bool exportPorkchopCsv(const char* path, const sim::LambertPorkchopResult& res, std::string* err) {
  if (!path || !*path) {
    if (err) *err = "Invalid path.";
    return false;
  }

  const int depSteps = res.departSteps;
  const int tofSteps = res.tofSteps;
  const std::size_t expected = (depSteps > 0 && tofSteps > 0) ? (std::size_t)depSteps * (std::size_t)tofSteps : 0u;
  if (expected == 0u) {
    if (err) *err = "No grid.";
    return false;
  }
  if (res.grid.size() < expected) {
    if (err) *err = "Grid not complete.";
    return false;
  }

  std::ofstream f(path, std::ios::out | std::ios::trunc);
  if (!f.is_open()) {
    if (err) *err = "Failed to open file.";
    return false;
  }

  f << "# departSteps," << depSteps << "\n";
  f << "# tofSteps," << tofSteps << "\n";
  f << "departAfterSec,tofSec,ok,dvDepartMagKmS,arriveRelSpeedKmS,totalDvKmS,score\n";

  for (std::size_t i = 0; i < expected; ++i) {
    const auto& c = res.grid[i];
    f << c.departAfterSec << ','
      << c.tofSec << ','
      << (c.ok ? 1 : 0) << ','
      << c.dvDepartMagKmS << ','
      << c.arriveRelSpeedKmS << ','
      << c.totalDvKmS << ','
      << c.score << '\n';
  }

  return true;
}

} // namespace stellar::game
