#include "PhotoModeWindow.h"

#include <imgui.h>

#include <algorithm>
#include <cctype>
#include <string>

namespace stellar::game {

static bool anyNonSpace(const char* s) {
  if (!s) return false;
  for (const unsigned char* p = (const unsigned char*)s; *p; ++p) {
    if (std::isspace(*p) == 0) return true;
  }
  return false;
}

void drawPhotoModeWindow(PhotoModeWindowState& st, const PhotoModeContext& ctx) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(520.0f, 380.0f), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Photo Mode", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextDisabled(
      "Take screenshots directly from the renderer. \n"
      "World-only captures happen before UI draws; UI captures include HUD/windows.");

  ImGui::Separator();

  ImGui::Text("Framebuffer: %d x %d", ctx.framebufferW, ctx.framebufferH);
  ImGui::SameLine();
  ImGui::TextDisabled("(%s)", ctx.paused ? "paused" : "running");

  ImGui::Checkbox("Pause simulation while open", &st.pauseWhileOpen);
  ImGui::Checkbox("Capture includes UI", &st.captureIncludeUi);
  ImGui::Checkbox("Append timestamp", &st.timestamp);
  ImGui::Checkbox("Copy saved path to clipboard", &st.copyPathToClipboard);

  ImGui::Separator();

  ImGui::TextDisabled("Output");

  ImGui::SetNextItemWidth(-FLT_MIN);
  if (st.focusDir) {
    ImGui::SetKeyboardFocusHere();
    st.focusDir = false;
  }
  ImGui::InputText("##shot_dir", st.outDir, (int)sizeof(st.outDir));
  ImGui::SameLine();
  ImGui::TextDisabled("Dir");

  ImGui::SetNextItemWidth(-FLT_MIN);
  ImGui::InputText("##shot_base", st.baseName, (int)sizeof(st.baseName));
  ImGui::SameLine();
  ImGui::TextDisabled("Base name");

  ImGui::Separator();

  const bool dirOk = anyNonSpace(st.outDir);
  const bool baseOk = anyNonSpace(st.baseName);

  if (!dirOk || !baseOk || ctx.framebufferW <= 0 || ctx.framebufferH <= 0) {
    ImGui::BeginDisabled();
  }
  if (ImGui::Button(st.captureIncludeUi ? "Take screenshot (with UI)" : "Take screenshot (world only)", ImVec2(-FLT_MIN, 0.0f))) {
    if (ctx.requestScreenshot) {
      ctx.requestScreenshot(st.captureIncludeUi,
                            std::string(st.outDir),
                            std::string(st.baseName),
                            st.timestamp,
                            st.copyPathToClipboard);
    }
  }
  if (!dirOk || !baseOk || ctx.framebufferW <= 0 || ctx.framebufferH <= 0) {
    ImGui::EndDisabled();
  }

  if (!dirOk) ImGui::TextDisabled("(Set an output directory)");
  if (!baseOk) ImGui::TextDisabled("(Set a base filename)");

  ImGui::Separator();

  if (!st.lastSavedPath.empty()) {
    ImGui::TextDisabled("Last saved:");
    ImGui::PushTextWrapPos();
    ImGui::TextUnformatted(st.lastSavedPath.c_str());
    ImGui::PopTextWrapPos();

    if (ImGui::Button("Copy path")) {
      ImGui::SetClipboardText(st.lastSavedPath.c_str());
      if (ctx.toast) ctx.toast("Copied screenshot path.", 1.4);
    }
  } else {
    ImGui::TextDisabled("No screenshot captured yet.");
  }

  ImGui::End();
}

} // namespace stellar::game
