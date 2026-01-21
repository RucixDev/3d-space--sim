#include "ProceduralLabWindow.h"

#include "Screenshot.h"

#include "imgui.h"

#include <SDL_opengl.h>

#include <algorithm>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <random>
#include <vector>

namespace stellar::game {

namespace fs = std::filesystem;

namespace {

core::u64 randomU64() {
  static std::random_device rd;
  static std::mt19937_64 rng((core::u64(rd()) << 32) ^ core::u64(rd()));
  return (core::u64)rng();
}

int opInputCount(render::ProcNodeOp op) {
  using render::ProcNodeOp;
  switch (op) {
    case ProcNodeOp::Constant:
    case ProcNodeOp::UvX:
    case ProcNodeOp::UvY:
    case ProcNodeOp::Noise2D:
    case ProcNodeOp::Fbm2D:
    case ProcNodeOp::Perlin2D:
    case ProcNodeOp::FbmPerlin2D:
    case ProcNodeOp::RidgedFbmPerlin2D:
    case ProcNodeOp::Voronoi2D:
      return 0;

    case ProcNodeOp::Abs:
    case ProcNodeOp::Invert:
    case ProcNodeOp::Fract:
    case ProcNodeOp::Clamp01:
    case ProcNodeOp::Smoothstep:
    case ProcNodeOp::Pow:
    case ProcNodeOp::Sine:
    case ProcNodeOp::Warp:
    case ProcNodeOp::Pan:
      return 1;

    case ProcNodeOp::Add:
    case ProcNodeOp::Sub:
    case ProcNodeOp::Mul:
    case ProcNodeOp::Div:
    case ProcNodeOp::Min:
    case ProcNodeOp::Max:
      return 2;

    default:
      return 0;
  }
}

void setOpDefaults(render::ProcNode& n) {
  using render::ProcNodeOp;
  n.in0 = -1;
  n.in1 = -1;
  n.in2 = -1;
  n.seed = 0;
  n.p0 = n.p1 = n.p2 = n.p3 = 0.0f;

  switch (n.op) {
    case ProcNodeOp::Constant:
      n.p0 = 0.5f;
      break;
    case ProcNodeOp::Noise2D:
      n.p0 = 8.0f;   // freq
      n.p1 = 0.0f;   // speed
      n.p2 = 1.0f;   // amp
      n.p3 = 0.0f;   // bias
      break;
    case ProcNodeOp::Perlin2D:
      n.p0 = 8.0f;   // freq
      n.p1 = 0.0f;   // speed
      n.p2 = 1.0f;   // amp
      n.p3 = 0.0f;   // bias
      break;
    case ProcNodeOp::Fbm2D:
      n.p0 = 6.0f;   // freq
      n.p1 = 5.0f;   // octaves
      n.p2 = 2.0f;   // lacunarity
      n.p3 = 0.5f;   // gain
      break;
    case ProcNodeOp::FbmPerlin2D:
    case ProcNodeOp::RidgedFbmPerlin2D:
      n.p0 = 6.0f;   // freq
      n.p1 = 5.0f;   // octaves
      n.p2 = 2.0f;   // lacunarity
      n.p3 = 0.5f;   // gain
      break;
    case ProcNodeOp::Voronoi2D:
      n.p0 = 12.0f;  // freq
      n.p1 = 1.0f;   // jitter
      break;
    case ProcNodeOp::Smoothstep:
      n.p0 = 0.3f;
      n.p1 = 0.7f;
      break;
    case ProcNodeOp::Pow:
      n.p0 = 1.5f;
      break;
    case ProcNodeOp::Sine:
      n.p0 = 12.0f; // frequency
      n.p1 = 0.0f;  // phase
      break;
    case ProcNodeOp::Warp:
      n.p0 = 0.25f; // strength
      n.p1 = 2.0f;  // freq
      n.p2 = 4.0f;  // octaves
      n.p3 = 0.55f; // gain
      break;
    case ProcNodeOp::Pan:
      n.p0 = 0.05f; // vx
      n.p1 = 0.00f; // vy
      break;
    default:
      break;
  }
}

bool readTextureRgba8(const render::Texture2D& tex, int w, int h, std::vector<unsigned char>& outPixels, std::string* outErr) {
  if (w <= 0 || h <= 0 || tex.handle() == 0) {
    if (outErr) *outErr = "readTextureRgba8: invalid texture/size";
    return false;
  }

  outPixels.resize((std::size_t)w * (std::size_t)h * 4u);

  glBindTexture(GL_TEXTURE_2D, tex.handle());
  glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, outPixels.data());
  glBindTexture(GL_TEXTURE_2D, 0);

  return true;
}

void applyPreset(ProceduralLabWindowState& st) {
  st.graph = render::makeProceduralGraphPreset(st.preset, st.seed);
  st.usePalette = st.graph.usePalette;
  st.dirty = true;
}

bool drawOpCombo(render::ProcNode& n) {
  using render::ProcNodeOp;
  bool changed = false;

  const char* preview = render::procNodeOpName(n.op);
  if (ImGui::BeginCombo("Op", preview)) {
    // Keep this list in a sensible order for humans.
    const ProcNodeOp ops[] = {
        ProcNodeOp::Constant,
        ProcNodeOp::UvX,
        ProcNodeOp::UvY,
        ProcNodeOp::Noise2D,
        ProcNodeOp::Fbm2D,
        ProcNodeOp::Perlin2D,
        ProcNodeOp::FbmPerlin2D,
        ProcNodeOp::RidgedFbmPerlin2D,
        ProcNodeOp::Voronoi2D,
        ProcNodeOp::Warp,
        ProcNodeOp::Pan,
        ProcNodeOp::Add,
        ProcNodeOp::Sub,
        ProcNodeOp::Mul,
        ProcNodeOp::Div,
        ProcNodeOp::Min,
        ProcNodeOp::Max,
        ProcNodeOp::Abs,
        ProcNodeOp::Invert,
        ProcNodeOp::Fract,
        ProcNodeOp::Clamp01,
        ProcNodeOp::Smoothstep,
        ProcNodeOp::Pow,
        ProcNodeOp::Sine,
    };

    for (auto op : ops) {
      const bool isSelected = (n.op == op);
      if (ImGui::Selectable(render::procNodeOpName(op), isSelected)) {
        n.op = op;
        setOpDefaults(n);
        changed = true;
      }
      if (isSelected) ImGui::SetItemDefaultFocus();
    }
    ImGui::EndCombo();
  }

  return changed;
}

bool drawNodeInputs(render::ProcGraph& g, int nodeIndex) {
  auto& n = g.nodes[(std::size_t)nodeIndex];
  const int inputs = opInputCount(n.op);
  bool changed = false;

  auto inputCombo = [&](const char* label, int& idx) {
    const char* preview = (idx < 0) ? "None" : render::procNodeOpName(g.nodes[(std::size_t)idx].op);
    if (ImGui::BeginCombo(label, preview)) {
      if (ImGui::Selectable("None", idx < 0)) {
        idx = -1;
        changed = true;
      }

      // Only allow connections to earlier nodes to keep mental model simple (and avoid recursion).
      for (int j = 0; j < nodeIndex; ++j) {
        const char* opName = render::procNodeOpName(g.nodes[(std::size_t)j].op);
        char buf[128];
        snprintf(buf, sizeof(buf), "#%d %s", j, opName);
        const bool isSel = (idx == j);
        if (ImGui::Selectable(buf, isSel)) {
          idx = j;
          changed = true;
        }
      }

      ImGui::EndCombo();
    }
  };

  if (inputs >= 1) inputCombo("In0", n.in0);
  if (inputs >= 2) inputCombo("In1", n.in1);
  if (inputs >= 3) inputCombo("In2", n.in2);

  return changed;
}

bool drawNodeParams(render::ProcNode& n) {
  using render::ProcNodeOp;
  bool changed = false;

  switch (n.op) {
    case ProcNodeOp::Constant:
      changed |= ImGui::SliderFloat("Value", &n.p0, 0.0f, 1.0f);
      break;

    case ProcNodeOp::Noise2D:
      changed |= ImGui::SliderFloat("Freq", &n.p0, 0.0f, 64.0f);
      changed |= ImGui::SliderFloat("Speed", &n.p1, -2.0f, 2.0f);
      changed |= ImGui::SliderFloat("Amp", &n.p2, 0.0f, 4.0f);
      changed |= ImGui::SliderFloat("Bias", &n.p3, -1.0f, 1.0f);
      break;

    case ProcNodeOp::Perlin2D:
      changed |= ImGui::SliderFloat("Freq", &n.p0, 0.0f, 64.0f);
      changed |= ImGui::SliderFloat("Speed", &n.p1, -2.0f, 2.0f);
      changed |= ImGui::SliderFloat("Amp", &n.p2, 0.0f, 4.0f);
      changed |= ImGui::SliderFloat("Bias", &n.p3, -1.0f, 1.0f);
      break;

    case ProcNodeOp::Fbm2D:
      changed |= ImGui::SliderFloat("Freq", &n.p0, 0.0f, 64.0f);
      changed |= ImGui::SliderFloat("Octaves", &n.p1, 1.0f, 8.0f);
      changed |= ImGui::SliderFloat("Lacunarity", &n.p2, 1.0f, 4.0f);
      changed |= ImGui::SliderFloat("Gain", &n.p3, 0.0f, 0.95f);
      break;

    case ProcNodeOp::FbmPerlin2D:
    case ProcNodeOp::RidgedFbmPerlin2D:
      changed |= ImGui::SliderFloat("Freq", &n.p0, 0.0f, 64.0f);
      changed |= ImGui::SliderFloat("Octaves", &n.p1, 1.0f, 8.0f);
      changed |= ImGui::SliderFloat("Lacunarity", &n.p2, 1.0f, 4.0f);
      changed |= ImGui::SliderFloat("Gain", &n.p3, 0.0f, 0.95f);
      break;

    case ProcNodeOp::Voronoi2D:
      changed |= ImGui::SliderFloat("Freq", &n.p0, 0.0f, 64.0f);
      changed |= ImGui::SliderFloat("Jitter", &n.p1, 0.0f, 1.0f);
      break;

    case ProcNodeOp::Smoothstep:
      changed |= ImGui::SliderFloat("Edge0", &n.p0, 0.0f, 1.0f);
      changed |= ImGui::SliderFloat("Edge1", &n.p1, 0.0f, 1.0f);
      break;

    case ProcNodeOp::Pow:
      changed |= ImGui::SliderFloat("Exp", &n.p0, 0.05f, 6.0f);
      break;

    case ProcNodeOp::Sine:
      changed |= ImGui::SliderFloat("Freq", &n.p0, 0.0f, 120.0f);
      changed |= ImGui::SliderFloat("Phase", &n.p1, -6.28318f, 6.28318f);
      break;

    case ProcNodeOp::Warp:
      changed |= ImGui::SliderFloat("Strength", &n.p0, 0.0f, 2.0f);
      changed |= ImGui::SliderFloat("Warp Freq", &n.p1, 0.0f, 32.0f);
      changed |= ImGui::SliderFloat("Warp Oct", &n.p2, 1.0f, 8.0f);
      changed |= ImGui::SliderFloat("Warp Gain", &n.p3, 0.0f, 0.95f);
      break;

    case ProcNodeOp::Pan:
      changed |= ImGui::SliderFloat("Vel X", &n.p0, -1.0f, 1.0f);
      changed |= ImGui::SliderFloat("Vel Y", &n.p1, -1.0f, 1.0f);
      break;

    default:
      break;
  }

  return changed;
}

} // namespace

void drawProceduralLabWindow(ProceduralLabWindowState& st, float timeSec, const ToastFn& toast) {
  if (!st.open) return;

  if (!ImGui::Begin("Procedural Lab", &st.open)) {
    ImGui::End();
    return;
  }

  // Keep graph state consistent with quick toggles.
  st.graph.seed = st.seed;
  st.graph.usePalette = st.usePalette;

  ImGui::TextUnformatted("GPU Procedural Texture Graph (bakes to a Texture2D via generated GLSL).");
  ImGui::Separator();

  // Preset selection
  {
    const char* preview = render::procGraphPresetName(st.preset);
    if (ImGui::BeginCombo("Preset", preview)) {
      const render::ProcGraphPreset presets[] = {
          render::ProcGraphPreset::Nebula,
          render::ProcGraphPreset::Marble,
          render::ProcGraphPreset::Lava,
          render::ProcGraphPreset::AlienCircuit,
          render::ProcGraphPreset::Rocky,
      };
      for (auto p : presets) {
        const bool sel = (st.preset == p);
        if (ImGui::Selectable(render::procGraphPresetName(p), sel)) {
          st.preset = p;
          if (st.lockToPreset) applyPreset(st);
        }
        if (sel) ImGui::SetItemDefaultFocus();
      }
      ImGui::EndCombo();
    }

    ImGui::SameLine();
    ImGui::Checkbox("Lock preset", &st.lockToPreset);

    ImGui::SameLine();
    if (ImGui::Button("Apply preset")) applyPreset(st);
  }

  // Seed / resolution / bake toggles
  {
    unsigned long long seedULL = (unsigned long long)st.seed;
    if (ImGui::InputScalar("Seed", ImGuiDataType_U64, &seedULL)) {
      st.seed = (core::u64)seedULL;
      if (st.lockToPreset) applyPreset(st);
      st.dirty = true;
    }

    ImGui::SameLine();
    if (ImGui::Button("Randomize")) {
      st.seed = randomU64();
      if (st.lockToPreset) applyPreset(st);
      st.dirty = true;
    }

    if (ImGui::SliderInt("Resolution", &st.resolution, 64, 2048)) st.dirty = true;

    if (ImGui::Checkbox("Generate mips", &st.bakeGenerateMips)) st.dirty = true;
    if (ImGui::SliderFloat("Dither", &st.bakeDitherStrength, 0.0f, 2.0f, "%.2f")) st.dirty = true;
    if (ImGui::Checkbox("Pack height in alpha", &st.bakePackHeightInAlpha)) st.dirty = true;

    // Apply quality options to the baker (affects the next bake).
    st.baker.setGenerateMips(st.bakeGenerateMips);
    st.baker.setDitherStrength(st.bakeDitherStrength);
    st.baker.setPackHeightInAlpha(st.bakePackHeightInAlpha);

    ImGui::Checkbox("Auto-bake", &st.autoBake);
    ImGui::SameLine();
    const bool wantBake = ImGui::Button("Bake now") || (st.autoBake && st.dirty);
    if (wantBake) {
      std::string err;
      if (!st.baker.bake(st.graph, st.resolution, st.resolution, timeSec, &err)) {
        st.lastError = err;
        st.dirty = false; // don't spam compile attempts on failure
      } else {
        st.lastError.clear();
        st.dirty = false;
      }
    }
  }

  // Preview / export
  ImGui::Separator();

  if (st.baker.isReady()) {
    const float preview = 320.0f;
    ImGui::Image((ImTextureID)(intptr_t)st.baker.texture().handle(), ImVec2(preview, preview), ImVec2(0, 1), ImVec2(1, 0));

    const auto& stats = st.baker.stats();
    if (st.bakeGenerateMips) {
      ImGui::Text("Shader: %s (build %.2f ms) | Draw %.2f ms | Mips %.2f ms",
                  stats.shaderRebuilt ? "rebuilt" : "cached",
                  stats.shaderBuildMs,
                  stats.drawMs,
                  stats.mipsGenerated ? stats.mipGenMs : 0.0);
    } else {
      ImGui::Text("Shader: %s (build %.2f ms) | Draw %.2f ms | Mips: off",
                  stats.shaderRebuilt ? "rebuilt" : "cached",
                  stats.shaderBuildMs,
                  stats.drawMs);
    }
  } else {
    ImGui::TextUnformatted("(No baked texture yet)");
  }

  if (!st.lastError.empty()) {
    ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "Error: %s", st.lastError.c_str());
  }

  ImGui::Separator();
  ImGui::InputText("Export path", st.exportPath, sizeof(st.exportPath));
  ImGui::SameLine();
  ImGui::Checkbox("Flip Y", &st.exportFlipY);

  if (ImGui::Button("Export PNG")) {
    if (!st.baker.isReady()) {
      toast("Nothing to export yet (bake a texture first).", 3.0);
    } else {
      std::vector<unsigned char> pixels;
      std::string err;
      if (!readTextureRgba8(st.baker.texture(), st.resolution, st.resolution, pixels, &err)) {
        toast(std::string("Export failed: ") + err, 4.0);
      } else {
        try {
          const fs::path outPath(st.exportPath);
          if (outPath.has_parent_path()) fs::create_directories(outPath.parent_path());

          std::string pngErr;
          const int stride = st.resolution * 4;
          if (!writePixelsToPng(outPath.string(), st.resolution, st.resolution, 4, pixels.data(), stride, st.exportFlipY, &pngErr)) {
            toast(std::string("Export failed: ") + pngErr, 4.0);
          } else {
            toast(std::string("Exported ") + outPath.string(), 3.0);
          }
        } catch (const std::exception& e) {
          toast(std::string("Export failed: ") + e.what(), 4.0);
        }
      }
    }
  }

  // Graph save/load (asset pipeline)
  ImGui::Separator();
  ImGui::TextUnformatted("Graph I/O (..../*.procgraph)");
  ImGui::InputText("Graph file", st.graphPath, sizeof(st.graphPath));

  if (ImGui::Button("Save graph")) {
    std::string err;
    if (!render::saveProcGraphToFile(st.graph, st.graphPath, &err)) {
      toast(std::string("Save failed: ") + err, 4.0);
    } else {
      toast(std::string("Saved ") + st.graphPath, 3.0);
    }
  }
  ImGui::SameLine();
  if (ImGui::Button("Load graph")) {
    render::ProcGraph loaded;
    std::string err;
    if (!render::loadProcGraphFromFile(st.graphPath, loaded, &err)) {
      toast(std::string("Load failed: ") + err, 4.0);
    } else {
      st.graph = std::move(loaded);
      st.seed = st.graph.seed;
      st.usePalette = st.graph.usePalette;
      st.lockToPreset = false;
      st.dirty = true;
      toast(std::string("Loaded ") + st.graphPath, 3.0);
    }
  }

  // Palette editor
  ImGui::Separator();
  if (ImGui::Checkbox("Use palette", &st.usePalette)) st.dirty = true;

  if (st.usePalette) {
    int palCount = st.graph.paletteCount;
    if (ImGui::SliderInt("Palette stops", &palCount, 2, render::kProcGraphMaxPaletteStops)) {
      st.graph.paletteCount = palCount;
      st.dirty = true;
    }

    for (int i = 0; i < st.graph.paletteCount; ++i) {
      auto& s = st.graph.palette[(std::size_t)i];
      ImGui::PushID(i);

      float col[3] = {s.r, s.g, s.b};
      if (ImGui::ColorEdit3("Color", col, ImGuiColorEditFlags_NoInputs)) {
        s.r = col[0];
        s.g = col[1];
        s.b = col[2];
        st.dirty = true;
      }

      if (ImGui::SliderFloat("Pos", &s.pos, 0.0f, 1.0f)) st.dirty = true;

      ImGui::Separator();
      ImGui::PopID();
    }

    if (ImGui::Button("Evenly distribute")) {
      const int n = std::max(2, st.graph.paletteCount);
      for (int i = 0; i < n; ++i) st.graph.palette[(std::size_t)i].pos = (float)i / (float)(n - 1);
      st.dirty = true;
    }
  }

  // Graph editor
  ImGui::Separator();
  ImGui::Text("Nodes (%d / %d)", (int)st.graph.nodes.size(), render::kProcGraphMaxNodes);

  if (st.graph.nodes.size() >= (std::size_t)render::kProcGraphMaxNodes) {
    ImGui::TextColored(ImVec4(1, 0.8f, 0.2f, 1), "Node limit reached (shader supports up to %d).", render::kProcGraphMaxNodes);
  }

  ImGui::BeginDisabled(st.lockToPreset);

  if (ImGui::Button("Add node")) {
    if (st.graph.nodes.size() < (std::size_t)render::kProcGraphMaxNodes) {
      render::ProcNode n{};
      n.op = render::ProcNodeOp::Fbm2D;
      setOpDefaults(n);
      st.graph.nodes.push_back(n);
      st.graph.output = (int)st.graph.nodes.size() - 1;
      st.dirty = true;
    }
  }

  ImGui::SameLine();
  if (ImGui::Button("Remove last")) {
    if (!st.graph.nodes.empty()) {
      st.graph.nodes.pop_back();
      if (st.graph.nodes.empty()) st.graph.output = -1;
      else st.graph.output = std::min(st.graph.output, (int)st.graph.nodes.size() - 1);
      st.dirty = true;
    }
  }

  // Output selection
  if (!st.graph.nodes.empty()) {
    int out = st.graph.output;
    if (out < 0 || out >= (int)st.graph.nodes.size()) out = (int)st.graph.nodes.size() - 1;

    char preview[128];
    snprintf(preview, sizeof(preview), "#%d %s", out, render::procNodeOpName(st.graph.nodes[(std::size_t)out].op));

    if (ImGui::BeginCombo("Output", preview)) {
      for (int i = 0; i < (int)st.graph.nodes.size(); ++i) {
        char label[128];
        snprintf(label, sizeof(label), "#%d %s", i, render::procNodeOpName(st.graph.nodes[(std::size_t)i].op));
        const bool sel = (out == i);
        if (ImGui::Selectable(label, sel)) {
          st.graph.output = i;
          st.dirty = true;
        }
        if (sel) ImGui::SetItemDefaultFocus();
      }
      ImGui::EndCombo();
    }
  }

  for (int i = 0; i < (int)st.graph.nodes.size(); ++i) {
    auto& n = st.graph.nodes[(std::size_t)i];
    ImGui::PushID(i);

    char nodeLabel[128];
    snprintf(nodeLabel, sizeof(nodeLabel), "#%d %s", i, render::procNodeOpName(n.op));

    if (ImGui::TreeNode(nodeLabel)) {
      bool nodeChanged = false;
      nodeChanged |= drawOpCombo(n);
      nodeChanged |= drawNodeInputs(st.graph, i);
      nodeChanged |= drawNodeParams(n);

      unsigned long long seedULL = (unsigned long long)n.seed;
      if (ImGui::InputScalar("Node seed tweak", ImGuiDataType_U64, &seedULL)) {
        n.seed = (core::u64)seedULL;
        nodeChanged = true;
      }

      if (nodeChanged) st.dirty = true;

      ImGui::TreePop();
    }

    ImGui::PopID();
  }

  ImGui::EndDisabled();

  // Shader source
  ImGui::Separator();
  ImGui::Checkbox("Show generated GLSL", &st.showShaderSource);
  if (st.showShaderSource && st.baker.isReady()) {
    ImGui::BeginChild("glsl_src", ImVec2(0, 220), true);
    ImGui::TextUnformatted(st.baker.lastFragmentSource().c_str());
    ImGui::EndChild();
  }

  ImGui::End();
}

} // namespace stellar::game
