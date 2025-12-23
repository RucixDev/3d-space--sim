#include "stellar/core/Log.h"

#include <atomic>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>

namespace stellar::core {

static std::atomic<LogLevel> g_level{LogLevel::Info};
static std::mutex g_logMutex;

void setLogLevel(LogLevel level) { g_level.store(level, std::memory_order_relaxed); }
LogLevel getLogLevel() { return g_level.load(std::memory_order_relaxed); }

std::string_view toString(LogLevel level) {
  switch (level) {
    case LogLevel::Trace: return "TRACE";
    case LogLevel::Debug: return "DEBUG";
    case LogLevel::Info:  return "INFO ";
    case LogLevel::Warn:  return "WARN ";
    case LogLevel::Error: return "ERROR";
    case LogLevel::Off:   return "OFF  ";
  }
  return "?????";
}

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

void log(LogLevel level, std::string_view message) {
  const LogLevel cur = getLogLevel();
  if (cur == LogLevel::Off) return;
  if (static_cast<int>(level) < static_cast<int>(cur)) return;

  std::lock_guard<std::mutex> lock(g_logMutex);
  std::cerr << "[" << timestampNow() << "][" << toString(level) << "] " << message << "\n";
}

} // namespace stellar::core
