#include "ConsoleWindow.h"

#include "stellar/core/CVar.h"

#include "stellar/ui/FuzzySearch.h"

#include <imgui.h>

#include <algorithm>
#include <chrono>
#include <ctime>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace stellar::game {

static std::string timestampNow() {
  using clock = std::chrono::system_clock;
  const auto now = clock::now();
  const auto t = clock::to_time_t(now);

  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif
  const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

  std::ostringstream oss;
  oss << std::put_time(&tm, "%H:%M:%S") << '.' << std::setw(3) << std::setfill('0') << ms.count();
  return oss.str();
}

static int levelIndex(core::LogLevel lvl) {
  switch (lvl) {
    case core::LogLevel::Trace: return 0;
    case core::LogLevel::Debug: return 1;
    case core::LogLevel::Info:  return 2;
    case core::LogLevel::Warn:  return 3;
    case core::LogLevel::Error: return 4;
    case core::LogLevel::Off:   return 2;
  }
  return 2;
}

static bool icontains(std::string_view hay, std::string_view needle) {
  if (needle.empty()) return true;
  // Simple ASCII case-insensitive substring.
  auto lower = [](unsigned char c) { return (unsigned char)std::tolower(c); };
  for (std::size_t i = 0; i < hay.size(); ++i) {
    std::size_t j = 0;
    while (i + j < hay.size() && j < needle.size() && lower((unsigned char)hay[i + j]) == lower((unsigned char)needle[j])) {
      ++j;
    }
    if (j == needle.size()) return true;
  }
  return false;
}

static void pushLine(ConsoleWindowState& st, ConsoleLine line) {
  std::lock_guard<std::mutex> lock(st.mutex);
  st.lines.push_back(std::move(line));
  while (st.lines.size() > st.maxLines) {
    st.lines.pop_front();
  }
  st.scrollToBottom = true;
}

void consolePrint(ConsoleWindowState& st, core::LogLevel lvl, std::string_view msg) {
  ConsoleLine line;
  line.level = lvl;
  line.timestamp = timestampNow();
  if (st.simTimeDaysPtr) {
    line.simTimeDays = *st.simTimeDaysPtr;
    line.hasSimTime = true;
  }
  line.text = std::string(msg);
  line.isCommand = false;
  pushLine(st, std::move(line));
}

void consoleClear(ConsoleWindowState& st) {
  std::lock_guard<std::mutex> lock(st.mutex);
  st.lines.clear();
  st.scrollToBottom = true;
}

void consoleAddCommand(ConsoleWindowState& st,
                       std::string name,
                       std::string help,
                       std::function<void(ConsoleWindowState&, const std::vector<std::string_view>& args)> fn) {
  ConsoleCommand cmd;
  cmd.name = std::move(name);
  cmd.help = std::move(help);
  cmd.fn = std::move(fn);
  st.commands.push_back(std::move(cmd));

  std::sort(st.commands.begin(), st.commands.end(), [](const ConsoleCommand& a, const ConsoleCommand& b) {
    return a.name < b.name;
  });
}

static core::LogLevel parseLogLevel(std::string_view s, bool* ok) {
  auto lower = [](std::string_view in) {
    std::string out;
    out.reserve(in.size());
    for (char ch : in) out.push_back((char)std::tolower((unsigned char)ch));
    return out;
  };
  const std::string v = lower(s);
  if (v == "trace") { if (ok) *ok = true; return core::LogLevel::Trace; }
  if (v == "debug") { if (ok) *ok = true; return core::LogLevel::Debug; }
  if (v == "info")  { if (ok) *ok = true; return core::LogLevel::Info; }
  if (v == "warn" || v == "warning") { if (ok) *ok = true; return core::LogLevel::Warn; }
  if (v == "error" || v == "err") { if (ok) *ok = true; return core::LogLevel::Error; }
  if (v == "off" || v == "none") { if (ok) *ok = true; return core::LogLevel::Off; }
  if (ok) *ok = false;
  return core::LogLevel::Info;
}

static void tokenize(std::string_view line, std::vector<std::string>& owned, std::vector<std::string_view>& out) {
  owned.clear();
  out.clear();

  std::size_t i = 0;
  auto skipWs = [&](void) {
    while (i < line.size() && std::isspace((unsigned char)line[i])) ++i;
  };

  skipWs();
  while (i < line.size()) {
    if (line[i] == '"') {
      // Quoted token
      ++i;
      std::string tok;
      while (i < line.size()) {
        char ch = line[i++];
        if (ch == '\\' && i < line.size()) {
          // Minimal escape support for \" and \\.
          char nxt = line[i++];
          tok.push_back(nxt);
          continue;
        }
        if (ch == '"') break;
        tok.push_back(ch);
      }
      owned.push_back(std::move(tok));
      out.push_back(owned.back());
    } else {
      // Unquoted token
      std::size_t start = i;
      while (i < line.size() && !std::isspace((unsigned char)line[i])) ++i;
      owned.emplace_back(line.substr(start, i - start));
      out.push_back(owned.back());
    }
    skipWs();
  }
}

static const ConsoleCommand* findCommand(const ConsoleWindowState& st, std::string_view name) {
  auto eqi = [](char a, char b) {
    return std::tolower((unsigned char)a) == std::tolower((unsigned char)b);
  };
  for (const auto& cmd : st.commands) {
    if (cmd.name.size() != name.size()) continue;
    bool ok = true;
    for (std::size_t i = 0; i < name.size(); ++i) {
      if (!eqi(cmd.name[i], name[i])) { ok = false; break; }
    }
    if (ok) return &cmd;
  }
  return nullptr;
}

static void execCommand(ConsoleWindowState& st, std::string_view line) {
  // Trim
  std::size_t b = 0;
  while (b < line.size() && std::isspace((unsigned char)line[b])) ++b;
  std::size_t e = line.size();
  while (e > b && std::isspace((unsigned char)line[e - 1])) --e;
  if (e <= b) return;

  const std::string trimmed(line.substr(b, e - b));

  // Record the command in the output.
  {
    ConsoleLine l;
    l.level = core::LogLevel::Info;
    l.timestamp = timestampNow();
    if (st.simTimeDaysPtr) { l.simTimeDays = *st.simTimeDaysPtr; l.hasSimTime = true; }
    l.text = "> " + trimmed;
    l.isCommand = true;
    pushLine(st, std::move(l));
  }

  // Update history (dedupe last occurrence).
  if (!trimmed.empty()) {
    for (auto it = st.history.begin(); it != st.history.end(); ++it) {
      if (*it == trimmed) { st.history.erase(it); break; }
    }
    st.history.push_back(trimmed);
    if (st.history.size() > 256) st.history.erase(st.history.begin());
  }
  st.historyPos = -1;

  // Parse
  std::vector<std::string> owned;
  std::vector<std::string_view> toks;
  tokenize(trimmed, owned, toks);
  if (toks.empty()) return;

  const std::string_view cmdName = toks[0];
  std::vector<std::string_view> args;
  if (toks.size() > 1) args.assign(toks.begin() + 1, toks.end());

  const ConsoleCommand* cmd = findCommand(st, cmdName);
  if (!cmd || !cmd->fn) {
    consolePrint(st, core::LogLevel::Warn, std::string("Unknown command: ") + std::string(cmdName));
    consolePrint(st, core::LogLevel::Info, "Type 'help' to list commands.");
    return;
  }

  cmd->fn(st, args);
}

// ---- Core log sink ----

static void coreLogSink(core::LogLevel level, std::string_view ts, std::string_view msg, void* user) {
  auto* st = reinterpret_cast<ConsoleWindowState*>(user);
  if (!st) return;

  ConsoleLine line;
  line.level = level;
  line.timestamp = std::string(ts);
  if (st->simTimeDaysPtr) { line.simTimeDays = *st->simTimeDaysPtr; line.hasSimTime = true; }
  line.text = std::string(msg);
  line.isCommand = false;
  pushLine(*st, std::move(line));
}

void consoleInstallCoreLogSink(ConsoleWindowState& st) {
  if (st.logSinkInstalled) return;
  st.logSinkInstalled = true;
  core::addLogSink(core::LogSink{&coreLogSink, &st});
}

void consoleRemoveCoreLogSink(ConsoleWindowState& st) {
  if (!st.logSinkInstalled) return;
  st.logSinkInstalled = false;
  core::removeLogSink(core::LogSink{&coreLogSink, &st});
}

// ---- Builtins ----

void consoleAddBuiltins(ConsoleWindowState& st) {
  core::installDefaultCVars();
  // help
  consoleAddCommand(st, "help", "List commands. Usage: help [command]",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      if (!args.empty()) {
        const ConsoleCommand* cmd = findCommand(c, args[0]);
        if (!cmd) {
          consolePrint(c, core::LogLevel::Warn, "Unknown command.");
          return;
        }
        consolePrint(c, core::LogLevel::Info, cmd->name + " - " + cmd->help);
        return;
      }
      consolePrint(c, core::LogLevel::Info, "Commands:");
      for (const auto& cmd : c.commands) {
        consolePrint(c, core::LogLevel::Info, "  " + cmd.name + " - " + cmd.help);
      }
      consolePrint(c, core::LogLevel::Info, "Tips: TAB to complete, ↑↓ for history.");
    });

  // clear
  consoleAddCommand(st, "clear", "Clear console output.",
    [](ConsoleWindowState& c, const std::vector<std::string_view>&) {
      consoleClear(c);
    });

  // echo
  consoleAddCommand(st, "echo", "Echo arguments back to the console.",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      std::string out;
      for (std::size_t i = 0; i < args.size(); ++i) {
        if (i) out.push_back(' ');
        out += args[i];
      }
      consolePrint(c, core::LogLevel::Info, out);
    });

  // loglevel
  consoleAddCommand(st, "loglevel", "Get/set log level. Usage: loglevel [trace|debug|info|warn|error|off]",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      if (args.empty()) {
        consolePrint(c, core::LogLevel::Info, std::string("Current log level: ") + std::string(core::toString(core::getLogLevel())));
        return;
      }
      bool ok = false;
      const core::LogLevel lvl = parseLogLevel(args[0], &ok);
      if (!ok) {
        consolePrint(c, core::LogLevel::Warn, "Unknown log level.");
        return;
      }
      core::setLogLevel(lvl);
      consolePrint(c, core::LogLevel::Info, "Log level updated.");
    });

  // ---- CVars ----
  // (Backed by stellar::core::CVarRegistry)
  consoleAddCommand(st, "cvars", "List CVars. Usage: cvars [filter]",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      core::installDefaultCVars();
      const std::string filter = args.empty() ? std::string() : std::string(args[0]);
      const auto vars = core::cvars().list(filter);

      if (vars.empty()) {
        consolePrint(c, core::LogLevel::Info, "No CVars.");
        return;
      }

      consolePrint(c, core::LogLevel::Info, "CVars:");
      for (const auto* v : vars) {
        std::string flags;
        if ((v->flags & core::CVar_Archive) != 0u) flags += 'A';
        if ((v->flags & core::CVar_ReadOnly) != 0u) flags += 'R';
        if ((v->flags & core::CVar_Cheat) != 0u) flags += 'C';

        std::string line = "  " + v->name + " = " + core::CVarRegistry::valueToString(*v);
        line += " (" + std::string(core::CVarRegistry::typeName(v->type)) + ")";
        if (!flags.empty()) line += " [" + flags + "]";
        if (!v->help.empty()) line += " - " + v->help;

        consolePrint(c, core::LogLevel::Info, line);
      }
    });

  consoleAddCommand(st, "cvar.get", "Get a CVar. Usage: cvar.get <name>",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      core::installDefaultCVars();
      if (args.empty()) {
        consolePrint(c, core::LogLevel::Warn, "Usage: cvar.get <name>");
        return;
      }

      const core::CVar* v = core::cvars().find(args[0]);
      if (!v) {
        consolePrint(c, core::LogLevel::Warn, "Unknown cvar.");
        const auto near = core::cvars().list(args[0]);
        if (!near.empty()) {
          consolePrint(c, core::LogLevel::Info, "Did you mean:");
          for (std::size_t i = 0; i < near.size() && i < 8; ++i) {
            consolePrint(c, core::LogLevel::Info, std::string("  ") + near[i]->name);
          }
        }
        return;
      }

      core::CVar tmp = *v;
      tmp.value = v->defaultValue;

      std::string out = v->name + " = " + core::CVarRegistry::valueToString(*v);
      out += " (" + std::string(core::CVarRegistry::typeName(v->type)) + ")";
      out += "  [default=" + core::CVarRegistry::valueToString(tmp) + "]";
      consolePrint(c, core::LogLevel::Info, out);
      if (!v->help.empty()) consolePrint(c, core::LogLevel::Info, std::string("  ") + v->help);
    });

  consoleAddCommand(st, "cvar.set", "Set a CVar. Usage: cvar.set <name> <value...>",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      core::installDefaultCVars();
      if (args.size() < 2) {
        consolePrint(c, core::LogLevel::Warn, "Usage: cvar.set <name> <value...>");
        return;
      }

      std::string value;
      for (std::size_t i = 1; i < args.size(); ++i) {
        if (i > 1) value.push_back(' ');
        value += args[i];
      }

      std::string err;
      if (!core::cvars().setFromString(args[0], value, &err)) {
        consolePrint(c, core::LogLevel::Warn, err);
        return;
      }

      const core::CVar* v = core::cvars().find(args[0]);
      if (v) {
        consolePrint(c, core::LogLevel::Info, v->name + " = " + core::CVarRegistry::valueToString(*v));
      } else {
        consolePrint(c, core::LogLevel::Info, "OK.");
      }
    });

  consoleAddCommand(st, "cvar.toggle", "Toggle a bool CVar. Usage: cvar.toggle <name>",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      core::installDefaultCVars();
      if (args.empty()) {
        consolePrint(c, core::LogLevel::Warn, "Usage: cvar.toggle <name>");
        return;
      }
      const core::CVar* v = core::cvars().find(args[0]);
      if (!v) {
        consolePrint(c, core::LogLevel::Warn, "Unknown cvar.");
        return;
      }
      if (v->type != core::CVarType::Bool) {
        consolePrint(c, core::LogLevel::Warn, "Not a bool cvar.");
        return;
      }
      const bool cur = std::get<bool>(v->value);
      std::string err;
      if (!core::cvars().setBool(v->name, !cur, &err)) {
        consolePrint(c, core::LogLevel::Warn, err);
        return;
      }
      const core::CVar* v2 = core::cvars().find(args[0]);
      if (v2) consolePrint(c, core::LogLevel::Info, v2->name + " = " + core::CVarRegistry::valueToString(*v2));
    });

  consoleAddCommand(st, "cvar.reset", "Reset a CVar to its default. Usage: cvar.reset <name>",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      core::installDefaultCVars();
      if (args.empty()) {
        consolePrint(c, core::LogLevel::Warn, "Usage: cvar.reset <name>");
        return;
      }
      std::string err;
      if (!core::cvars().reset(args[0], &err)) {
        consolePrint(c, core::LogLevel::Warn, err);
        return;
      }
      const core::CVar* v = core::cvars().find(args[0]);
      if (v) consolePrint(c, core::LogLevel::Info, v->name + " = " + core::CVarRegistry::valueToString(*v));
    });

  consoleAddCommand(st, "cvar.save", "Save archived CVars. Usage: cvar.save [path]",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      core::installDefaultCVars();
      const std::string path = args.empty() ? std::string("cvars.cfg") : std::string(args[0]);
      std::string err;
      if (!core::cvars().saveFile(path, &err)) {
        consolePrint(c, core::LogLevel::Warn, err);
        return;
      }
      consolePrint(c, core::LogLevel::Info, std::string("Saved CVars to ") + path);
    });

  consoleAddCommand(st, "cvar.load", "Load CVars from a file. Usage: cvar.load [path]",
    [](ConsoleWindowState& c, const std::vector<std::string_view>& args) {
      core::installDefaultCVars();
      const std::string path = args.empty() ? std::string("cvars.cfg") : std::string(args[0]);
      std::string err;
      if (!core::cvars().loadFile(path, &err)) {
        consolePrint(c, core::LogLevel::Warn, err.empty() ? "Failed to load." : err);
        return;
      }
      consolePrint(c, core::LogLevel::Info, std::string("Loaded CVars from ") + path);
    });
}

// ---- ImGui input callbacks ----

static int Stricmp(const char* s1, const char* s2) {
  int d;
  while ((d = std::tolower((unsigned char)*s2) - std::tolower((unsigned char)*s1)) == 0 && *s1) {
    s1++;
    s2++;
  }
  return d;
}

static int Strnicmp(const char* s1, const char* s2, int n) {
  int d = 0;
  while (n > 0 && (d = std::tolower((unsigned char)*s2) - std::tolower((unsigned char)*s1)) == 0 && *s1) {
    s1++;
    s2++;
    n--;
  }
  return d;
}

static int TextEditCallback(ImGuiInputTextCallbackData* data) {
  auto* st = reinterpret_cast<ConsoleWindowState*>(data->UserData);
  if (!st) return 0;

  switch (data->EventFlag) {
    case ImGuiInputTextFlags_CallbackCompletion: {
      // Find beginning of current word.
      const char* word_end = data->Buf + data->CursorPos;
      const char* word_start = word_end;
      while (word_start > data->Buf) {
        const char c = word_start[-1];
        if (c == ' ' || c == '\t' || c == ',') break;
        word_start--;
      }

      // Determine which token we're completing (0 = command, 1 = first argument, ...).
      int token_index = 0;
      bool in_tok = false;
      for (const char* p = data->Buf; p < word_start; ++p) {
        const char c = *p;
        const bool sep = (c == ' ' || c == '\t' || c == ',');
        if (sep) {
          if (in_tok) { in_tok = false; token_index++; }
        } else {
          in_tok = true;
        }
      }

      const int len = (int)(word_end - word_start);
      std::vector<const char*> candidates;

      auto complete = [&](bool add_space) {
        if (candidates.empty()) {
          consolePrint(*st, core::LogLevel::Warn, "No match.");
        } else if (candidates.size() == 1) {
          // Single match: complete and optionally add a space.
          data->DeleteChars((int)(word_start - data->Buf), len);
          data->InsertChars(data->CursorPos, candidates[0]);
          if (add_space) data->InsertChars(data->CursorPos, " ");
        } else {
          // Multiple matches: complete as much as possible.
          int match_len = len;
          for (;;) {
            int c = 0;
            bool all_matches = true;
            for (std::size_t i = 0; i < candidates.size() && all_matches; i++) {
              if (i == 0)
                c = std::tolower((unsigned char)candidates[i][match_len]);
              else if (c == 0 || c != std::tolower((unsigned char)candidates[i][match_len]))
                all_matches = false;
            }
            if (!all_matches)
              break;
            match_len++;
          }

          if (match_len > len) {
            data->DeleteChars((int)(word_start - data->Buf), len);
            data->InsertChars(data->CursorPos, candidates[0], candidates[0] + match_len);
          }

          consolePrint(*st, core::LogLevel::Info, "Possible matches:");
          std::vector<const char*> sorted = candidates;
          std::sort(sorted.begin(), sorted.end(), [](const char* a, const char* b) { return Stricmp(a, b) < 0; });
          for (const char* cand : sorted) {
            consolePrint(*st, core::LogLevel::Info, std::string("  ") + cand);
          }
        }
      };

      if (token_index == 0) {
        // Complete command name.
        candidates.reserve(st->commands.size());
        for (const auto& cmd : st->commands) {
          if (Strnicmp(cmd.name.c_str(), word_start, len) == 0) {
            candidates.push_back(cmd.name.c_str());
          }
        }
        complete(true);
      } else if (token_index == 1) {
        // Complete CVar names for cvar.* commands.
        const char* cmd_end = data->Buf;
        while (cmd_end < word_start) {
          const char c = *cmd_end;
          if (c == ' ' || c == '\t' || c == ',') break;
          ++cmd_end;
        }
        const std::string cmd_name(data->Buf, cmd_end);
        const bool wants_cvar =
          Stricmp(cmd_name.c_str(), "cvar.get") == 0 ||
          Stricmp(cmd_name.c_str(), "cvar.set") == 0 ||
          Stricmp(cmd_name.c_str(), "cvar.toggle") == 0 ||
          Stricmp(cmd_name.c_str(), "cvar.reset") == 0;

        if (wants_cvar) {
          core::installDefaultCVars();
          const auto vars = core::cvars().list();
          candidates.reserve(vars.size());
          for (const auto* v : vars) {
            if (Strnicmp(v->name.c_str(), word_start, len) == 0) {
              candidates.push_back(v->name.c_str());
            }
          }
          complete(true);
        }
      }

      break;
    }

    case ImGuiInputTextFlags_CallbackHistory: {
      const int prev_history_pos = st->historyPos;
      if (data->EventKey == ImGuiKey_UpArrow) {
        if (st->historyPos == -1)
          st->historyPos = (int)st->history.size() - 1;
        else if (st->historyPos > 0)
          st->historyPos--;
      } else if (data->EventKey == ImGuiKey_DownArrow) {
        if (st->historyPos != -1)
          if (++st->historyPos >= (int)st->history.size())
            st->historyPos = -1;
      }

      if (prev_history_pos != st->historyPos) {
        const char* history_str = (st->historyPos >= 0) ? st->history[(std::size_t)st->historyPos].c_str() : "";
        data->DeleteChars(0, data->BufTextLen);
        data->InsertChars(0, history_str);
      }
      break;
    }
  }

  return 0;
}

// ---- UI drawing ----

static std::string formatLine(const ConsoleLine& l, bool showWallClock, bool showSimTime, bool showLevel) {
  std::string out;
  out.reserve(l.text.size() + 48);

  if (showWallClock && !l.timestamp.empty()) {
    out += '[';
    out += l.timestamp;
    out += "] ";
  }

  if (showSimTime && l.hasSimTime) {
    // Display in days with 3 decimals (good for long sims).
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << l.simTimeDays;
    out += "(D";
    out += oss.str();
    out += ") ";
  }

  if (showLevel) {
    out += '[';
    out += std::string(core::toString(l.level));
    out += "] ";
  }

  out += l.text;
  return out;
}

static std::string buildVisibleDump(const ConsoleWindowState& st) {
  std::string dump;
  dump.reserve(st.lines.size() * 80);

  const std::string_view f = st.filter;
  for (const auto& l : st.lines) {
    const int idx = levelIndex(l.level);
    if (idx >= 0 && idx < 5 && !st.levelEnabled[idx]) continue;
    if (!f.empty()) {
      if (!icontains(l.text, f) && !icontains(l.timestamp, f)) continue;
    }
    dump += formatLine(l, st.showWallClock, st.showSimTime, st.showLevel);
    dump.push_back('\n');
  }
  return dump;
}

void drawConsoleWindow(ConsoleWindowState& st) {
  if (!st.open) return;

  ImGui::SetNextWindowSize(ImVec2(720.0f, 420.0f), ImGuiCond_FirstUseEver);
  if (!ImGui::Begin("Console", &st.open)) {
    ImGui::End();
    return;
  }

  // Shortcuts: Ctrl+L to clear.
  if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows)) {
    if (ImGui::GetIO().KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_L)) {
      consoleClear(st);
    }
  }

  // Toolbar
  if (ImGui::Button("Clear")) {
    consoleClear(st);
  }
  ImGui::SameLine();
  if (ImGui::Button("Copy")) {
    const std::string dump = buildVisibleDump(st);
    ImGui::SetClipboardText(dump.c_str());
    consolePrint(st, core::LogLevel::Info, "Copied visible log to clipboard.");
  }
  ImGui::SameLine();
  if (ImGui::Button("Save")) {
    const std::string dump = buildVisibleDump(st);
    const char* path = "console_log.txt";
    std::ofstream ofs(path, std::ios::binary);
    if (ofs) {
      ofs.write(dump.data(), (std::streamsize)dump.size());
      consolePrint(st, core::LogLevel::Info, std::string("Saved to ") + path);
    } else {
      consolePrint(st, core::LogLevel::Warn, "Failed to save console_log.txt");
    }
  }

  ImGui::SameLine();
  ImGui::Checkbox("Auto-scroll", &st.autoScroll);
  ImGui::SameLine();
  ImGui::Checkbox("Wrap", &st.wrapLines);
  ImGui::SameLine();
  ImGui::Checkbox("Wall clock", &st.showWallClock);
  ImGui::SameLine();
  ImGui::Checkbox("Sim time", &st.showSimTime);
  ImGui::SameLine();
  ImGui::Checkbox("Level", &st.showLevel);

  ImGui::Separator();

  // Filters
  ImGui::TextDisabled("Filter:");
  ImGui::SameLine();
  ImGui::SetNextItemWidth(260.0f);
  ImGui::InputTextWithHint("##console_filter", "substring...", st.filter, sizeof(st.filter));

  ImGui::SameLine();
  ImGui::TextDisabled("Levels:");
  ImGui::SameLine();
  ImGui::Checkbox("T", &st.levelEnabled[0]);
  ImGui::SameLine();
  ImGui::Checkbox("D", &st.levelEnabled[1]);
  ImGui::SameLine();
  ImGui::Checkbox("I", &st.levelEnabled[2]);
  ImGui::SameLine();
  ImGui::Checkbox("W", &st.levelEnabled[3]);
  ImGui::SameLine();
  ImGui::Checkbox("E", &st.levelEnabled[4]);

  ImGui::SameLine();
  ImGui::TextDisabled("(TAB complete | ↑↓ history | Ctrl+L clear)");

  ImGui::Separator();

  // Output region
  const float footerH = ImGui::GetStyle().ItemSpacing.y + ImGui::GetFrameHeightWithSpacing();
  if (ImGui::BeginChild("##console_scroller", ImVec2(0, -footerH), true,
                        ImGuiWindowFlags_HorizontalScrollbar)) {
    if (st.wrapLines) {
      ImGui::PushTextWrapPos(0.0f);
    }

    std::lock_guard<std::mutex> lock(st.mutex);

    // We keep it simple: iterate and render. maxLines is small enough.
    const std::string_view f = st.filter;
    for (const auto& l : st.lines) {
      const int idx = levelIndex(l.level);
      if (idx >= 0 && idx < 5 && !st.levelEnabled[idx]) continue;
      if (!f.empty()) {
        if (!icontains(l.text, f) && !icontains(l.timestamp, f)) continue;
      }

      const std::string formatted = formatLine(l, st.showWallClock, st.showSimTime, st.showLevel);
      ImGui::TextUnformatted(formatted.c_str());
    }

    if (st.wrapLines) {
      ImGui::PopTextWrapPos();
    }

    if (st.scrollToBottom || (st.autoScroll && ImGui::GetScrollY() >= ImGui::GetScrollMaxY())) {
      ImGui::SetScrollHereY(1.0f);
      st.scrollToBottom = false;
    }
  }
  ImGui::EndChild();

  // Input
  ImGui::Separator();

  ImGuiInputTextFlags flags = ImGuiInputTextFlags_EnterReturnsTrue
                              | ImGuiInputTextFlags_CallbackCompletion
                              | ImGuiInputTextFlags_CallbackHistory
                              | ImGuiInputTextFlags_AutoSelectAll;

  if (st.focusInput) {
    ImGui::SetKeyboardFocusHere();
    st.focusInput = false;
  }

  if (ImGui::InputText("Input", st.input, sizeof(st.input), flags, &TextEditCallback, &st)) {
    execCommand(st, st.input);
    st.input[0] = '\0';
    st.focusInput = true;
  }

  ImGui::End();
}

} // namespace stellar::game
