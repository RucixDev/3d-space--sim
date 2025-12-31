#include "stellar/core/Log.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <vector>

namespace stellar::core {

static std::atomic<LogLevel> g_level{LogLevel::Info};
static std::mutex g_logMutex;
static std::vector<LogSink> g_sinks;

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

void addLogSink(LogSink sink) {
  if (!sink.fn) return;
  std::lock_guard<std::mutex> lock(g_logMutex);
  g_sinks.push_back(sink);
}

void removeLogSink(LogSink sink) {
  std::lock_guard<std::mutex> lock(g_logMutex);
  g_sinks.erase(std::remove_if(g_sinks.begin(), g_sinks.end(), [&](const LogSink& s) {
    return s.fn == sink.fn && s.user == sink.user;
  }), g_sinks.end());
}

void log(LogLevel level, std::string_view message) {
  const LogLevel cur = getLogLevel();
  if (cur == LogLevel::Off) return;
  if (static_cast<int>(level) < static_cast<int>(cur)) return;

  const std::string ts = timestampNow();

  // Copy sinks under the mutex, then invoke outside the lock to avoid
  // re-entrancy/deadlocks if a sink itself logs.
  std::vector<LogSink> sinksCopy;
  {
    std::lock_guard<std::mutex> lock(g_logMutex);
    std::cerr << "[" << ts << "][" << toString(level) << "] " << message << "\n";
    sinksCopy = g_sinks;
  }

  for (const LogSink& s : sinksCopy) {
    if (!s.fn) continue;
    s.fn(level, ts, message, s.user);
  }
}

} // namespace stellar::core
