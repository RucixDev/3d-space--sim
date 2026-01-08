#include "ActionWheel.h"

#include <algorithm>
#include <cmath>

namespace stellar::game {

void openActionWheel(ActionWheelState& st) {
  st.open = true;
  st.justOpened = true;
  st.requestClose = false;
  st.hovered = -1;
  st.page = 0;
  st.openT = 0.0f;
}

void closeActionWheel(ActionWheelState& st) {
  st.open = false;
  st.justOpened = false;
  st.requestClose = false;
  st.hovered = -1;
  st.openT = 0.0f;
}

static float wrapAnglePos(float a) {
  const float twoPi = IM_PI * 2.0f;
  while (a < 0.0f) a += twoPi;
  while (a >= twoPi) a -= twoPi;
  return a;
}

static int computeHovered(const ActionWheelState& st, int count, ImVec2 mouse) {
  if (count <= 0) return -1;

  const float dx = mouse.x - st.center.x;
  const float dy = mouse.y - st.center.y;
  const float r2 = dx * dx + dy * dy;
  const float r = std::sqrt(std::max(0.0f, r2));
  if (r < st.deadzoneRadius || r > st.outerRadius + 6.0f) return -1;

  // We want segment 0 to be at the top (12 o'clock).
  const float start = -IM_PI * 0.5f;
  float ang = std::atan2(dy, dx);
  ang = wrapAnglePos(ang - start);

  const float slice = (IM_PI * 2.0f) / (float)count;
  int idx = (int)(ang / slice);
  idx = std::clamp(idx, 0, count - 1);
  return idx;
}

static ImVec4 mulAlpha(ImVec4 c, float a) {
  c.w *= a;
  return c;
}

static ImU32 colU32(ImVec4 c) {
  return ImGui::ColorConvertFloat4ToU32(c);
}

int drawActionWheelOverlay(ActionWheelState& st,
                          const std::vector<ActionWheelItem>& items,
                          const char* pageName,
                          int pageIndex,
                          int pageCount) {
  if (!st.open) return -1;

  ImGuiIO& io = ImGui::GetIO();

  // Smooth fade-in.
  st.openT = std::min(1.0f, st.openT + io.DeltaTime * 12.0f);
  const float fade = st.openT;

  // Fullscreen overlay window in the foreground.
  ImGuiViewport* vp = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(vp->WorkPos);
  ImGui::SetNextWindowSize(vp->WorkSize);
  ImGui::SetNextWindowViewport(vp->ID);

  ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoDocking |
                           ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings |
                           ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoFocusOnAppearing |
                           ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse;

  ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
  ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
  ImGui::PushStyleColor(ImGuiCol_WindowBg, mulAlpha(ImGui::GetStyleColorVec4(ImGuiCol_WindowBg), 0.0f));

  int confirm = -1;

  if (ImGui::Begin("##action_wheel_overlay", nullptr, flags)) {
    if (st.justOpened) {
      st.center = vp->GetWorkCenter();
      st.justOpened = false;
    }

    // Capture click anywhere so underlying windows don't receive it.
    ImGui::SetCursorScreenPos(vp->WorkPos);
    ImGui::InvisibleButton("##aw_capture", vp->WorkSize,
                           ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight);
    const bool click = ImGui::IsMouseReleased(ImGuiMouseButton_Left);

    ImDrawList* dl = ImGui::GetWindowDrawList();

    // Dim background.
    const ImVec2 p0 = vp->WorkPos;
    const ImVec2 p1 = ImVec2(vp->WorkPos.x + vp->WorkSize.x, vp->WorkPos.y + vp->WorkSize.y);
    ImVec4 dim = ImGui::GetStyleColorVec4(ImGuiCol_ModalWindowDimBg);
    if (dim.w <= 0.01f) dim = ImVec4(0.0f, 0.0f, 0.0f, 0.42f);
    dim.w *= 0.85f;
    dl->AddRectFilled(p0, p1, colU32(mulAlpha(dim, fade)));

    const int n = (int)items.size();
    st.hovered = computeHovered(st, n, io.MousePos);

    const float start = -IM_PI * 0.5f;
    const float slice = (n > 0) ? (IM_PI * 2.0f) / (float)n : IM_PI * 2.0f;

    // Colors from style.
    ImVec4 colBase = ImGui::GetStyleColorVec4(ImGuiCol_FrameBg);
    ImVec4 colHover = ImGui::GetStyleColorVec4(ImGuiCol_ButtonHovered);
    ImVec4 colDisabled = colBase;
    colDisabled.w *= 0.45f;
    ImVec4 colOutline = ImGui::GetStyleColorVec4(ImGuiCol_Border);

    const int arcSeg = 28;

    // Ring sectors.
    for (int i = 0; i < n; ++i) {
      const ActionWheelItem& it = items[(std::size_t)i];
      const bool enabled = it.enabled && (bool)it.action;
      const bool hovered = (i == st.hovered);

      ImVec4 col = enabled ? colBase : colDisabled;
      if (hovered) col = enabled ? colHover : colDisabled;
      col = mulAlpha(col, 0.95f * fade);

      const float a0 = start + slice * (float)i;
      const float a1 = start + slice * (float)(i + 1);

      dl->PathClear();
      dl->PathArcTo(st.center, st.outerRadius, a0, a1, arcSeg);
      dl->PathArcTo(st.center, st.innerRadius, a1, a0, arcSeg);
      dl->PathFillConvex(colU32(col));

      dl->PathClear();
      dl->PathArcTo(st.center, st.outerRadius, a0, a1, arcSeg);
      dl->PathArcTo(st.center, st.innerRadius, a1, a0, arcSeg);
      dl->PathStroke(colU32(mulAlpha(colOutline, 0.75f * fade)), true, 1.0f);
    }

    // Center disk.
    {
      ImVec4 centerCol = ImGui::GetStyleColorVec4(ImGuiCol_WindowBg);
      centerCol.w *= 0.88f;
      dl->AddCircleFilled(st.center, st.innerRadius - 6.0f, colU32(mulAlpha(centerCol, fade)), 32);
      dl->AddCircle(st.center, st.innerRadius - 6.0f, colU32(mulAlpha(colOutline, 0.75f * fade)), 32, 1.0f);
    }

    // Labels.
    for (int i = 0; i < n; ++i) {
      const ActionWheelItem& it = items[(std::size_t)i];
      const bool enabled = it.enabled && (bool)it.action;
      const bool hovered = (i == st.hovered);

      const float mid = start + slice * ((float)i + 0.5f);
      const ImVec2 pos = ImVec2(st.center.x + std::cos(mid) * st.labelRadius,
                                st.center.y + std::sin(mid) * st.labelRadius);

      ImVec4 textCol = enabled ? ImGui::GetStyleColorVec4(ImGuiCol_Text)
                               : ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled);
      if (hovered && enabled) {
        textCol.w = std::min(1.0f, textCol.w * 1.0f);
      }
      textCol = mulAlpha(textCol, fade);

      std::string line = it.label;
      if (!it.shortcut.empty()) {
        line += "\n";
        line += it.shortcut;
      }

      const ImVec2 sz = ImGui::CalcTextSize(line.c_str());
      const ImVec2 tpos = ImVec2(pos.x - sz.x * 0.5f, pos.y - sz.y * 0.5f);
      dl->AddText(tpos, colU32(textCol), line.c_str());
    }

    // Center info.
    {
      std::string title = "Action Wheel";
      std::string sub;

      if (pageCount > 0) {
        if (pageName && pageName[0]) {
          sub = std::string(pageName) + "  (" + std::to_string(pageIndex + 1) + "/" + std::to_string(pageCount) + ")";
        } else {
          sub = "Page " + std::to_string(pageIndex + 1) + "/" + std::to_string(pageCount);
        }
      }

      if (st.hovered >= 0 && st.hovered < (int)items.size()) {
        const auto& it = items[(std::size_t)st.hovered];
        title = it.label;
        if (!it.detail.empty()) sub = it.detail;
      } else {
        if (sub.empty()) sub = "Scroll: switch page • Release: run";
      }

      ImVec4 titleCol = mulAlpha(ImGui::GetStyleColorVec4(ImGuiCol_Text), fade);
      ImVec4 subCol = mulAlpha(ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled), fade);

      const ImVec2 tSz = ImGui::CalcTextSize(title.c_str());
      const ImVec2 tPos = ImVec2(st.center.x - tSz.x * 0.5f, st.center.y - tSz.y * 0.65f);
      dl->AddText(tPos, colU32(titleCol), title.c_str());

      if (!sub.empty()) {
        const ImVec2 sSz = ImGui::CalcTextSize(sub.c_str());
        const ImVec2 sPos = ImVec2(st.center.x - sSz.x * 0.5f, st.center.y + 6.0f);
        dl->AddText(sPos, colU32(subCol), sub.c_str());
      }
    }

    // Confirm via click.
    if (click && st.hovered >= 0 && st.hovered < (int)items.size()) {
      const auto& it = items[(std::size_t)st.hovered];
      if (it.enabled && (bool)it.action) {
        confirm = st.hovered;
      }
    }
  }
  ImGui::End();

  ImGui::PopStyleColor();
  ImGui::PopStyleVar(2);

  return confirm;
}

} // namespace stellar::game
