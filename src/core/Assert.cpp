#include "stellar/core/Assert.h"
#include "stellar/core/Log.h"

#include <sstream>

namespace stellar::core {

[[noreturn]] void panic(std::string_view message, const char* file, int line) {
  std::ostringstream oss;
  oss << "PANIC: " << message << " (" << file << ":" << line << ")";
  log(LogLevel::Error, oss.str());
  std::abort();
}

} // namespace stellar::core
