#include "test_registry.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstring>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace {

static bool icontains(std::string_view hay, std::string_view needle) {
  if (needle.empty()) return true;
  if (hay.empty()) return false;

  auto lower = [](char c) -> char {
    return (char)std::tolower((unsigned char)c);
  };

  for (std::size_t i = 0; i < hay.size(); ++i) {
    std::size_t j = 0;
    while (i + j < hay.size() && j < needle.size() && lower(hay[i + j]) == lower(needle[j])) {
      ++j;
    }
    if (j == needle.size()) return true;
  }
  return false;
}

static void printHelp(const char* exe) {
  std::cout << "stellar_tests\n"
            << "Usage:\n"
            << "  " << (exe ? exe : "stellar_tests") << " [options] [pattern ...]\n\n"
            << "Options:\n"
            << "  --list                 List all available tests and exit\n"
            << "  --filter <pattern>     Run only tests whose name contains <pattern> (case-insensitive)\n"
            << "  --failfast             Stop after the first failing test\n"
            << "  -h, --help             Show this help\n\n"
            << "Patterns:\n"
            << "  Any extra positional arguments are treated as filter patterns (OR'd together).\n"
            << "Examples:\n"
            << "  stellar_tests --list\n"
            << "  stellar_tests --filter trade\n"
            << "  stellar_tests orbit docking\n";
}

} // namespace

int main(int argc, char** argv) {
  bool list = false;
  bool failfast = false;
  std::vector<std::string> patterns;

  for (int i = 1; i < argc; ++i) {
    const char* a = argv[i];
    if (!a) continue;

    if (std::strcmp(a, "--list") == 0) {
      list = true;
      continue;
    }
    if (std::strcmp(a, "--failfast") == 0) {
      failfast = true;
      continue;
    }
    if (std::strcmp(a, "--filter") == 0) {
      if (i + 1 >= argc) {
        std::cerr << "[stellar_tests] --filter requires a pattern\n";
        return 1;
      }
      patterns.emplace_back(argv[++i]);
      continue;
    }
    if (std::strcmp(a, "--help") == 0 || std::strcmp(a, "-h") == 0) {
      printHelp(argv[0]);
      return 0;
    }

    // Positional patterns.
    patterns.emplace_back(a);
  }

  std::size_t count = 0;
  const stellar::test::Case* cases = stellar::test::cases(&count);
  if (!cases || count == 0) {
    std::cerr << "[stellar_tests] No tests registered.\n";
    return 1;
  }

  if (list) {
    for (std::size_t i = 0; i < count; ++i) {
      std::cout << cases[i].name << "\n";
    }
    return 0;
  }

  auto matches = [&](const char* name) -> bool {
    if (!name) return false;
    if (patterns.empty()) return true;
    for (const auto& p : patterns) {
      if (icontains(name, p)) return true;
    }
    return false;
  };

  using clock = std::chrono::high_resolution_clock;

  int totalFails = 0;
  int ran = 0;
  int skipped = 0;

  for (std::size_t i = 0; i < count; ++i) {
    const stellar::test::Case& tc = cases[i];
    if (!matches(tc.name)) {
      ++skipped;
      continue;
    }

    ++ran;
    const auto t0 = clock::now();

    int fails = 0;
    if (!tc.fn) {
      fails = 1;
      std::cerr << "[stellar_tests] [FAIL] " << (tc.name ? tc.name : "<null>")
                << " (null function pointer)\n";
    } else {
      try {
        fails = tc.fn();
      } catch (const std::exception& e) {
        fails = 1;
        std::cerr << "[stellar_tests] [EXC ] " << (tc.name ? tc.name : "<null>")
                  << ": " << e.what() << "\n";
      } catch (...) {
        fails = 1;
        std::cerr << "[stellar_tests] [EXC ] " << (tc.name ? tc.name : "<null>")
                  << ": unknown exception\n";
      }
    }

    const auto t1 = clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    if (fails == 0) {
      std::cout << "[stellar_tests] [PASS] " << tc.name << " (" << ms << " ms)\n";
    } else {
      std::cerr << "[stellar_tests] [FAIL] " << tc.name << " fails=" << fails
                << " (" << ms << " ms)\n";
      totalFails += fails;
      if (failfast) break;
    }
  }

  if (ran == 0) {
    std::cerr << "[stellar_tests] No tests matched the given pattern(s).\n";
    printHelp(argv[0]);
    return 1;
  }

  if (totalFails == 0) {
    std::cout << "[stellar_tests] ALL PASS (" << ran << " tests";
    if (skipped > 0) std::cout << ", " << skipped << " skipped";
    std::cout << ")\n";
    return 0;
  }

  std::cerr << "[stellar_tests] FAILS=" << totalFails << " (" << ran << " tests";
  if (skipped > 0) std::cerr << ", " << skipped << " skipped";
  std::cerr << ")\n";
  return 1;
}
