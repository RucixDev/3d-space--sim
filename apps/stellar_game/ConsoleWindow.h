#pragma once

#include "stellar/core/Log.h"

#include <deque>
#include <functional>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>

namespace stellar::game {

// Lightweight in-game console window:
//  - Scrollable log output (can hook into core::log via LogSink).
//  - Filter + per-level visibility.
//  - Command input with history and tab completion.

struct ConsoleLine {
  stellar::core::LogLevel level{stellar::core::LogLevel::Info};
  std::string timestamp; // wall-clock timestamp ("HH:MM:SS.mmm"), may be empty
  double simTimeDays{0.0};
  bool hasSimTime{false};
  std::string text;
  bool isCommand{false}; // user-entered command line
};

struct ConsoleWindowState;

struct ConsoleCommand {
  std::string name;
  std::string help;
  std::function<void(ConsoleWindowState&, const std::vector<std::string_view>& args)> fn;
};

struct ConsoleWindowState {
  bool open{false};
  bool focusInput{false};

  // Display options
  bool autoScroll{true};
  bool showWallClock{true};
  bool showSimTime{false};
  bool showLevel{true};
  bool wrapLines{false};

  // Trace, Debug, Info, Warn, Error
  bool levelEnabled[5]{true, true, true, true, true};

  char filter[128]{};
  char input[256]{};

  // Command history
  std::vector<std::string> history;
  int historyPos{-1}; // -1 = new line
  std::size_t maxHistory{256};

  // History persistence (optional).
  // Path is treated as UTF-8 on all platforms.
  char historyPath[256]{"console_history.txt"};
  bool historyAutoSaveOnExit{true};
  bool historyDirty{false};
  bool historyLoaded{false};

  // Registered commands
  std::vector<ConsoleCommand> commands;

  // Output lines
  std::deque<ConsoleLine> lines;
  std::size_t maxLines{4000};

  // Optional pointer for stamping sim-time onto messages.
  const double* simTimeDaysPtr{nullptr};

  // Core log sink integration
  bool logSinkInstalled{false};

  // Internal UI flags
  bool scrollToBottom{false};

  // Protect concurrent writes (e.g., if log() is called off-thread).
  mutable std::mutex mutex;
};

// Append a line to the console.
void consolePrint(ConsoleWindowState& st, stellar::core::LogLevel lvl, std::string_view msg);

// Clear output/history.
void consoleClear(ConsoleWindowState& st);

// ---- History persistence ----
// Save/load command history. The format is a simple line-based text file (one command per line).
// Lines starting with '#' are treated as comments.
// Returns false and sets outError on failure.
bool consoleSaveHistoryToFile(const ConsoleWindowState& st, const std::string& path, std::string* outError = nullptr);
bool consoleLoadHistoryFromFile(ConsoleWindowState& st, const std::string& path, std::string* outError = nullptr);

// Clear only the command history (does not clear output lines).
void consoleClearHistory(ConsoleWindowState& st);

// Register a command.
void consoleAddCommand(ConsoleWindowState& st,
                       std::string name,
                       std::string help,
                       std::function<void(ConsoleWindowState&, const std::vector<std::string_view>& args)> fn);

// Add built-in commands: help, clear, echo, loglevel.
void consoleAddBuiltins(ConsoleWindowState& st);

// Execute a command line programmatically (same path as typing into the console).
// Useful for integrating other UI (e.g., a command palette) with the console.
void consoleExecLine(ConsoleWindowState& st, std::string_view line);

// Install/remove a core::Log sink that forwards logs into the console.
void consoleInstallCoreLogSink(ConsoleWindowState& st);
void consoleRemoveCoreLogSink(ConsoleWindowState& st);

// Draw the console window. (No-op if st.open == false.)
void drawConsoleWindow(ConsoleWindowState& st);

} // namespace stellar::game
