#pragma once

#include <string_view>

namespace stellar::core {

enum class LogLevel {
  Trace = 0,
  Debug = 1,
  Info  = 2,
  Warn  = 3,
  Error = 4,
  Off   = 5
};

void setLogLevel(LogLevel level);
LogLevel getLogLevel();

std::string_view toString(LogLevel level);

// Optional callback sink for log messages.
//
// Sinks are invoked after the message has been written to stderr.
// The timestamp and message views are only valid for the duration of the callback.
struct LogSink {
  using Fn = void (*)(LogLevel level, std::string_view timestamp, std::string_view message, void* user);
  Fn fn{nullptr};
  void* user{nullptr};
};

// Register/unregister a sink.
//
// Notes:
//  - This is a minimal facility intended for tooling/UI integration.
//  - Sinks obey the current log level filter.
//  - addLogSink() is idempotent only if you avoid registering duplicates.
void addLogSink(LogSink sink);
void removeLogSink(LogSink sink);

// Thread-safe enough for a starter: writes to stderr with a timestamp.
void log(LogLevel level, std::string_view message);

} // namespace stellar::core

#define STELLAR_LOG_TRACE(msg) ::stellar::core::log(::stellar::core::LogLevel::Trace, (msg))
#define STELLAR_LOG_DEBUG(msg) ::stellar::core::log(::stellar::core::LogLevel::Debug, (msg))
#define STELLAR_LOG_INFO(msg)  ::stellar::core::log(::stellar::core::LogLevel::Info,  (msg))
#define STELLAR_LOG_WARN(msg)  ::stellar::core::log(::stellar::core::LogLevel::Warn,  (msg))
#define STELLAR_LOG_ERROR(msg) ::stellar::core::log(::stellar::core::LogLevel::Error, (msg))
