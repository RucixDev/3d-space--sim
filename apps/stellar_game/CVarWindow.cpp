#include "CVarWindow.h"

#include "stellar/core/CVar.h"
#include "stellar/ui/FuzzySearch.h"

#include <imgui.h>

#include <algorithm>
#include <cstdio>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::game {
namespace {

static bool anyNonSpace(const char* s) {
  if (!s) return false;
  for (const unsigned char* p = (const unsigned char*)s; *p; ++p) {
    if (std::isspace(*p) == 0) return true;
  }
  return false;
}

static std::string trimAscii(std::string_view sv) {
  std::size_t b = 0;
  while (b < sv.size() && std::isspace((unsigned char)sv[b])) ++b;
  std::size_t e = sv.size();
  while (e > b && std::isspace((unsigned char)sv[e - 1])) --e;
  return std::string(sv.substr(b, e - b));
}

static std::string flagsString(std::uint32_t flags) {
  std::string s;
  if (flags & stellar::core::CVar_Archive)  s += 'A';
  if (flags & stellar::core::CVar_ReadOnly) s += 'R';
  if (flags & stellar::core::CVar_Cheat)    s += 'C';
  if (s.empty()) s = "-";
  return s;
}

static std::string valueToString(const stellar::core::CVarValue& v) {
  // Keep this lightweight; we only need UI-friendly text.
  if (const bool* b = std::get_if<bool>(&v)) {
    return *b ? "true" : "false";
  }
  if (const std::int64_t* i = std::get_if<std::int64_t>(&v)) {
    return std::to_string(*i);
  }
  if (const double* d = std::get_if<double>(&v)) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.6g", *d);
    return std::string(buf);
  }
  if (const std::string* s = std::get_if<std::string>(&v)) {
    return *s;
  }
  return {};
}

static bool parseSetLine(std::string_view line, std::string& outName, std::string& outValue) {
  const std::string trimmed = trimAscii(line);
  if (trimmed.empty()) return false;

  // Prefer `name = value` parsing.
  const std::size_t eq = trimmed.find('=');
  if (eq != std::string::npos) {
    outName = trimAscii(std::string_view(trimmed).substr(0, eq));
    outValue = trimAscii(std::string_view(trimmed).substr(eq + 1));
    return !outName.empty();
  }

  // Fallback: split on first whitespace.
  std::size_t sp = 0;
  while (sp < trimmed.size() && !std::isspace((unsigned char)trimmed[sp])) ++sp;
  outName = trimmed.substr(0, sp);
  while (sp < trimmed.size() && std::isspace((unsigned char)trimmed[sp])) ++sp;
  outValue = trimmed.substr(sp);
  return !outName.empty() && !outValue.empty();
}

static void helpTooltip(const stellar::core::CVar& v) {
  if (!ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) return;
  ImGui::BeginTooltip();
  ImGui::TextUnformatted(v.name.c_str());
  if (!v.help.empty()) {
    ImGui::Separator();
    ImGui::TextWrapped("%s", v.help.c_str());
  }
  ImGui::Separator();
  ImGui::TextDisabled("Type: %s", stellar::core::CVarRegistry::typeName(v.type));
  ImGui::TextDisabled("Default: %s", valueToString(v.defaultValue).c_str());
  ImGui::TextDisabled("Flags: %s", flagsString(v.flags).c_str());
  ImGui::EndTooltip();
}

} // namespace

void openCVarWindow(CVarWindowState& st, bool focusFilter) {
  st.open = true;
  st.focusFilter = focusFilter;
}

void drawCVarWindow(CVarWindowState& st, const ToastFn& toast) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(860.0f, 620.0f), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("CVars", &st.open)) {
    ImGui::End();
    return;
  }

  ImGui::TextUnformatted("Runtime-tweakable console variables (CVars).");
  ImGui::TextDisabled("Tip: You can also use the in-game console commands: cvar.get / cvar.set / cvar.save / cvar.load");

  // Filter row
  ImGui::PushItemWidth(280.0f);
  if (st.focusFilter) {
    ImGui::SetKeyboardFocusHere();
    st.focusFilter = false;
  }
  ImGui::InputTextWithHint("##cvar_filter", "Filter (fuzzy)...", st.filter, sizeof(st.filter));
  ImGui::PopItemWidth();
  ImGui::SameLine();
  ImGui::Checkbox("Relevance sort", &st.sortByRelevance);
  ImGui::SameLine();
  ImGui::Checkbox("Changed only", &st.showChangedOnly);
  ImGui::SameLine();
  ImGui::Checkbox("Archive only", &st.showArchivedOnly);
  ImGui::SameLine();
  ImGui::Checkbox("Show cheats", &st.showCheats);

  // Quick set
  {
    ImGui::PushItemWidth(-140.0f);
    bool apply = ImGui::InputTextWithHint(
      "##cvar_set", "Set: name = value   (press Enter)", st.setLine, sizeof(st.setLine),
      ImGuiInputTextFlags_EnterReturnsTrue);
    ImGui::PopItemWidth();
    ImGui::SameLine();
    if (ImGui::Button("Apply")) apply = true;

    if (apply) {
      std::string name, value;
      if (parseSetLine(st.setLine, name, value)) {
        std::string err;
        if (!stellar::core::cvars().setFromString(name, value, &err)) {
          if (toast) toast(std::string("CVar error: ") + err, 2.6);
        } else {
          if (toast) toast(std::string("Set ") + name + " = " + value, 1.8);
          st.setLine[0] = '\0';
        }
      }
    }
  }

  // Load/save
  {
    ImGui::Separator();
    ImGui::PushItemWidth(220.0f);
    ImGui::InputTextWithHint("##cvar_cfg", "cvars.cfg", st.configPath, sizeof(st.configPath));
    ImGui::PopItemWidth();
    ImGui::SameLine();
    if (ImGui::Button("Save")) {
      std::string err;
      const bool ok = stellar::core::cvars().saveFile(st.configPath, &err);
      if (toast) toast(ok ? "Saved CVars." : (std::string("Save failed: ") + err), 2.2);
    }
    ImGui::SameLine();
    if (ImGui::Button("Load")) {
      std::string err;
      const bool ok = stellar::core::cvars().loadFile(st.configPath, &err);
      if (toast) toast(ok ? "Loaded CVars." : (std::string("Load failed: ") + err), 2.2);
    }
  }

  ImGui::Separator();

  const bool hasQuery = anyNonSpace(st.filter);
  const auto vars = stellar::core::cvars().list();

  struct Row {
    const stellar::core::CVar* v{nullptr};
    int score{0};
  };

  std::vector<Row> rows;
  rows.reserve(vars.size());

  for (const stellar::core::CVar* v : vars) {
    if (!v) continue;
    if (!st.showCheats && (v->flags & stellar::core::CVar_Cheat) != 0u) continue;
    if (st.showArchivedOnly && (v->flags & stellar::core::CVar_Archive) == 0u) continue;
    if (st.showChangedOnly && v->value == v->defaultValue) continue;

    int score = 0;
    if (hasQuery) {
      // Match against name + help to make discovery easier.
      std::string hay;
      hay.reserve(v->name.size() + v->help.size() + 1);
      hay += v->name;
      hay.push_back(' ');
      hay += v->help;
      score = stellar::ui::fuzzyMatchScore(st.filter, hay);
      if (score < 0) continue;
    }
    rows.push_back(Row{v, score});
  }

  if (hasQuery && st.sortByRelevance) {
    std::stable_sort(rows.begin(), rows.end(), [&](const Row& a, const Row& b) {
      if (a.score != b.score) return a.score > b.score;
      return a.v->name < b.v->name;
    });
  }

  ImGui::TextDisabled("Showing %d / %d", (int)rows.size(), (int)vars.size());

  const ImGuiTableFlags tableFlags =
    ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp |
    ImGuiTableFlags_Resizable | ImGuiTableFlags_ScrollY;

  const float tableHeight = std::max(120.0f, ImGui::GetContentRegionAvail().y - 12.0f);
  if (ImGui::BeginTable("##cvars_table", 6, tableFlags, ImVec2(0.0f, tableHeight))) {
    ImGui::TableSetupColumn("Name", ImGuiTableColumnFlags_WidthStretch, 0.30f);
    ImGui::TableSetupColumn("Type", ImGuiTableColumnFlags_WidthFixed, 72.0f);
    ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthStretch, 0.25f);
    ImGui::TableSetupColumn("Default", ImGuiTableColumnFlags_WidthStretch, 0.20f);
    ImGui::TableSetupColumn("Flags", ImGuiTableColumnFlags_WidthFixed, 56.0f);
    ImGui::TableSetupColumn("Actions", ImGuiTableColumnFlags_WidthFixed, 112.0f);
    ImGui::TableHeadersRow();

    for (const Row& r : rows) {
      const stellar::core::CVar& v = *r.v;
      const bool readOnly = (v.flags & stellar::core::CVar_ReadOnly) != 0u;

      ImGui::TableNextRow();

      // Name
      ImGui::TableSetColumnIndex(0);
      if (v.value != v.defaultValue) {
        ImGui::TextUnformatted(v.name.c_str());
      } else {
        ImGui::TextDisabled("%s", v.name.c_str());
      }
      helpTooltip(v);

      // Type
      ImGui::TableSetColumnIndex(1);
      ImGui::TextDisabled("%s", stellar::core::CVarRegistry::typeName(v.type));

      // Value (editable)
      ImGui::TableSetColumnIndex(2);
      if (readOnly) ImGui::BeginDisabled();

      std::string err;
      const std::string idBase = std::string("##cvar_") + v.name;
      switch (v.type) {
        case stellar::core::CVarType::Bool: {
          bool cur = std::get<bool>(v.value);
          if (ImGui::Checkbox((idBase + "_b").c_str(), &cur)) {
            if (!stellar::core::cvars().setBool(v.name, cur, &err)) {
              if (toast) toast(std::string("CVar error: ") + err, 2.4);
            }
          }
        } break;

        case stellar::core::CVarType::Int: {
          std::int64_t cur = std::get<std::int64_t>(v.value);
          ImGui::SetNextItemWidth(-FLT_MIN);
          if (ImGui::InputScalar((idBase + "_i").c_str(), ImGuiDataType_S64, &cur, nullptr, nullptr, "%lld")) {
            if (!stellar::core::cvars().setInt(v.name, cur, &err)) {
              if (toast) toast(std::string("CVar error: ") + err, 2.4);
            }
          }
        } break;

        case stellar::core::CVarType::Float: {
          double cur = std::get<double>(v.value);
          ImGui::SetNextItemWidth(-FLT_MIN);
          if (ImGui::InputDouble((idBase + "_f").c_str(), &cur, 0.0, 0.0, "%.6g")) {
            if (!stellar::core::cvars().setFloat(v.name, cur, &err)) {
              if (toast) toast(std::string("CVar error: ") + err, 2.4);
            }
          }
        } break;

        case stellar::core::CVarType::String: {
          const std::string cur = std::get<std::string>(v.value);
          std::string shown = cur;
          if (shown.size() > 64) {
            shown.resize(61);
            shown += "...";
          }
          ImGui::TextUnformatted(shown.c_str());
          ImGui::SameLine();
          if (ImGui::SmallButton(("Edit" + idBase).c_str())) {
            st.editName = v.name;
            std::memset(st.editBuf, 0, sizeof(st.editBuf));
            std::snprintf(st.editBuf, sizeof(st.editBuf), "%s", cur.c_str());
            st.requestEditPopup = true;
          }
        } break;
      }

      if (readOnly) ImGui::EndDisabled();

      // Default
      ImGui::TableSetColumnIndex(3);
      const std::string defStr = valueToString(v.defaultValue);
      if (defStr.size() > 64) {
        std::string trunc = defStr.substr(0, 61) + "...";
        ImGui::TextDisabled("%s", trunc.c_str());
      } else {
        ImGui::TextDisabled("%s", defStr.c_str());
      }

      // Flags
      ImGui::TableSetColumnIndex(4);
      ImGui::TextDisabled("%s", flagsString(v.flags).c_str());

      // Actions
      ImGui::TableSetColumnIndex(5);
      if (ImGui::SmallButton(("Reset" + idBase).c_str())) {
        std::string e;
        if (!stellar::core::cvars().reset(v.name, &e)) {
          if (toast) toast(std::string("Reset failed: ") + e, 2.4);
        } else {
          if (toast) toast(std::string("Reset ") + v.name, 1.5);
        }
      }
      ImGui::SameLine();
      if (ImGui::SmallButton(("Copy" + idBase).c_str())) {
        const std::string line = v.name + " = " + valueToString(v.value);
        ImGui::SetClipboardText(line.c_str());
        if (toast) toast("Copied to clipboard.", 1.2);
      }
    }

    ImGui::EndTable();
  }

  // String edit popup (modal)
  if (st.requestEditPopup) {
    ImGui::OpenPopup("Edit CVar String");
    st.requestEditPopup = false;
  }
  bool openPopup = true;
  if (ImGui::BeginPopupModal("Edit CVar String", &openPopup, ImGuiWindowFlags_AlwaysAutoResize)) {
    ImGui::TextUnformatted(st.editName.c_str());
    ImGui::Separator();
    ImGui::InputTextMultiline("##cvar_str", st.editBuf, sizeof(st.editBuf), ImVec2(560.0f, 180.0f));

    if (ImGui::Button("Apply")) {
      std::string err;
      if (!stellar::core::cvars().setString(st.editName, std::string(st.editBuf), &err)) {
        if (toast) toast(std::string("CVar error: ") + err, 2.4);
      } else {
        if (toast) toast(std::string("Set ") + st.editName, 1.5);
        ImGui::CloseCurrentPopup();
      }
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel")) {
      ImGui::CloseCurrentPopup();
    }

    ImGui::EndPopup();
  }

  ImGui::End();
}

} // namespace stellar::game
