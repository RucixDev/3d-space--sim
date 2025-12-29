#pragma once

// Small compatibility helpers for Dear ImGui API differences and safer overloads.
// We deliberately keep this header tiny and dependency-light.

#include <imgui.h>

#include <string_view>

#ifndef IM_PI
// Dear ImGui defines IM_PI in imgui_internal.h, but user code should not depend on that header.
// Provide it here for legacy code paths that used IM_PI directly.
#define IM_PI 3.14159265358979323846f
#endif

namespace ImGui {

// Some codebases expect an ImGui::GetMouseDelta() helper.
// Dear ImGui exposes mouse delta via ImGuiIO::MouseDelta.
inline ImVec2 GetMouseDelta() {
  return GetIO().MouseDelta;
}

// Convenience overload for std::string_view so callers can avoid allocating std::string.
// Uses the (text, text_end) overload, which does not require null-termination.
inline void TextUnformatted(std::string_view sv) {
  TextUnformatted(sv.data(), sv.data() + sv.size());
}

} // namespace ImGui
