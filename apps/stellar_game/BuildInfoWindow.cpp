#include "BuildInfoWindow.h"

#include "stellar/core/Log.h"

#include <imgui.h>

#include <SDL.h>
#include <SDL_opengl.h>

#include <algorithm>
#include <cstdint>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>

namespace stellar::game {

namespace {

static const char* yesNo(bool v) { return v ? "Yes" : "No"; }

static std::string compilerString() {
#if defined(__clang__)
  std::ostringstream ss;
  ss << "Clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
  return ss.str();
#elif defined(_MSC_VER)
  std::ostringstream ss;
  ss << "MSVC " << _MSC_VER;
  return ss.str();
#elif defined(__GNUC__)
  std::ostringstream ss;
  ss << "GCC " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
  return ss.str();
#else
  return "Unknown compiler";
#endif
}

static std::string platformString() {
#if defined(_WIN32)
  return "Windows";
#elif defined(__APPLE__)
  return "macOS";
#elif defined(__linux__)
  return "Linux";
#elif defined(__EMSCRIPTEN__)
  return "Emscripten";
#else
  return "Unknown platform";
#endif
}

static std::string archString() {
#if defined(_M_X64) || defined(__x86_64__)
  return "x86_64";
#elif defined(_M_IX86) || defined(__i386__)
  return "x86";
#elif defined(__aarch64__) || defined(_M_ARM64)
  return "arm64";
#elif defined(__arm__) || defined(_M_ARM)
  return "arm";
#else
  return "unknown";
#endif
}

static std::string buildConfigString() {
#if defined(NDEBUG)
  return "Release (NDEBUG)";
#elif defined(_DEBUG)
  return "Debug (_DEBUG)";
#else
  return "Debug";
#endif
}

static std::string cppStdString() {
  // __cplusplus uses yyyymmL for C++11+.
  std::ostringstream ss;
  ss << __cplusplus;
  return ss.str();
}

static std::string safeGlString(GLenum name) {
  const GLubyte* s = glGetString(name);
  if (!s) return std::string();
  return reinterpret_cast<const char*>(s);
}

static std::string buildInfoClipboardText(bool includeGl) {
  std::ostringstream out;

  out << "Stellar Forge - Build Info\n";
  out << "Platform: " << platformString() << " (" << archString() << ")\n";
  out << "Compiler: " << compilerString() << "\n";
  out << "C++: __cplusplus=" << cppStdString() << "\n";
  out << "Config: " << buildConfigString() << "\n";
  out << "Threads: " << std::max(1u, std::thread::hardware_concurrency()) << "\n";

  out << "\nFeatures\n";
  out << "  STELLAR_ENABLE_RENDER=" << (int)STELLAR_ENABLE_RENDER << "\n";
  out << "  STELLAR_ENABLE_IMGUI=" << (int)STELLAR_ENABLE_IMGUI << "\n";

  out << "\nDependencies\n";
  out << "  SDL: " << SDL_MAJOR_VERSION << "." << SDL_MINOR_VERSION << "." << SDL_PATCHLEVEL << "\n";
#ifdef IMGUI_VERSION
  out << "  ImGui: " << IMGUI_VERSION << "\n";
#endif
#ifdef IMGUI_HAS_DOCK
  out << "  ImGui Docking: yes\n";
#else
  out << "  ImGui Docking: no\n";
#endif

  if (includeGl) {
    const std::string glVendor = safeGlString(GL_VENDOR);
    const std::string glRenderer = safeGlString(GL_RENDERER);
    const std::string glVersion = safeGlString(GL_VERSION);
    if (!glVendor.empty() || !glRenderer.empty() || !glVersion.empty()) {
      out << "\nOpenGL\n";
      if (!glVendor.empty()) out << "  Vendor: " << glVendor << "\n";
      if (!glRenderer.empty()) out << "  Renderer: " << glRenderer << "\n";
      if (!glVersion.empty()) out << "  Version: " << glVersion << "\n";
    }
  }

  return out.str();
}

} // namespace

void drawBuildInfoWindow(BuildInfoWindowState& st, const ToastFn& toast) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(560.0f, 420.0f), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Build Info", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextUnformatted("Stellar Forge");
  ImGui::SameLine();
  ImGui::TextDisabled("(build + dependency snapshot)");

  ImGui::Separator();

  if (ImGui::Button("Copy summary to clipboard")) {
    const std::string clip = buildInfoClipboardText(/*includeGl=*/st.showRuntimeGlInfo);
    ImGui::SetClipboardText(clip.c_str());
    if (toast) toast("Build info copied.", 1.4);
  }
  ImGui::SameLine();
  ImGui::Checkbox("Include OpenGL", &st.showRuntimeGlInfo);

  ImGui::Separator();

  if (ImGui::CollapsingHeader("Compile-time", ImGuiTreeNodeFlags_DefaultOpen)) {
    st.showCompileTimeInfo = true;
    ImGui::Text("Platform: %s", platformString().c_str());
    ImGui::Text("Arch: %s", archString().c_str());
    ImGui::Text("Compiler: %s", compilerString().c_str());
    ImGui::Text("C++: __cplusplus=%s", cppStdString().c_str());
    ImGui::Text("Config: %s", buildConfigString().c_str());
    ImGui::Text("Build time: %s %s", __DATE__, __TIME__);
    ImGui::Text("Hardware threads: %u", std::max(1u, std::thread::hardware_concurrency()));
  }

  if (ImGui::CollapsingHeader("Feature flags", ImGuiTreeNodeFlags_DefaultOpen)) {
    st.showFeatureFlags = true;
    ImGui::BulletText("STELLAR_ENABLE_RENDER: %d", (int)STELLAR_ENABLE_RENDER);
    ImGui::BulletText("STELLAR_ENABLE_IMGUI: %d", (int)STELLAR_ENABLE_IMGUI);
#ifdef IMGUI_HAS_DOCK
    ImGui::BulletText("IMGUI_HAS_DOCK: %s", yesNo(true));
#else
    ImGui::BulletText("IMGUI_HAS_DOCK: %s", yesNo(false));
#endif
  }

  if (ImGui::CollapsingHeader("Dependency versions", ImGuiTreeNodeFlags_DefaultOpen)) {
    st.showDependencyVersions = true;
    ImGui::BulletText("SDL: %d.%d.%d", SDL_MAJOR_VERSION, SDL_MINOR_VERSION, SDL_PATCHLEVEL);
#ifdef IMGUI_VERSION
    ImGui::BulletText("ImGui: %s", IMGUI_VERSION);
#endif
  }

  if (st.showRuntimeGlInfo) {
    if (ImGui::CollapsingHeader("OpenGL runtime", ImGuiTreeNodeFlags_DefaultOpen)) {
      const std::string glVendor = safeGlString(GL_VENDOR);
      const std::string glRenderer = safeGlString(GL_RENDERER);
      const std::string glVersion = safeGlString(GL_VERSION);
      if (glVendor.empty() && glRenderer.empty() && glVersion.empty()) {
        ImGui::TextDisabled("OpenGL context not active yet.");
      } else {
        if (!glVendor.empty()) ImGui::Text("Vendor: %s", glVendor.c_str());
        if (!glRenderer.empty()) ImGui::Text("Renderer: %s", glRenderer.c_str());
        if (!glVersion.empty()) ImGui::Text("Version: %s", glVersion.c_str());
      }
    }
  }

  ImGui::End();
}

} // namespace stellar::game
