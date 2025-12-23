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

// Thread-safe enough for a starter: writes to stderr with a timestamp.
void log(LogLevel level, std::string_view message);

} // namespace stellar::core

#define STELLAR_LOG_TRACE(msg) ::stellar::core::log(::stellar::core::LogLevel::Trace, (msg))
#define STELLAR_LOG_DEBUG(msg) ::stellar::core::log(::stellar::core::LogLevel::Debug, (msg))
#define STELLAR_LOG_INFO(msg)  ::stellar::core::log(::stellar::core::LogLevel::Info,  (msg))
#define STELLAR_LOG_WARN(msg)  ::stellar::core::log(::stellar::core::LogLevel::Warn,  (msg))
#define STELLAR_LOG_ERROR(msg) ::stellar::core::log(::stellar::core::LogLevel::Error, (msg))
