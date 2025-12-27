#include "stellar/core/Args.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

static std::vector<char*> makeArgv(std::initializer_list<const char*> items) {
  std::vector<char*> argv;
  argv.reserve(items.size());
  for (const char* s : items) {
    argv.push_back(const_cast<char*>(s));
  }
  return argv;
}

int test_args() {
  int fails = 0;

  using stellar::core::Args;

  // Negative numeric values should be accepted as values for --key and not be misread as switches.
  {
    auto argv = makeArgv({"app",
                          "--pos", "-1", "-2.5", "3",
                          "--radius", "-50.5"});
    Args args;
    args.setArity("pos", 3);
    args.parse((int)argv.size(), argv.data());

    const auto pos = args.values("pos");
    if (pos.size() != 3 || pos[0] != "-1" || pos[1] != "-2.5" || pos[2] != "3") {
      std::cerr << "[test_args] expected --pos -1 -2.5 3, got size=" << pos.size() << "\n";
      ++fails;
    }

    double radius = 0.0;
    if (!args.getDouble("radius", radius) || std::abs(radius - (-50.5)) > 1e-9) {
      std::cerr << "[test_args] expected --radius -50.5 to parse\n";
      ++fails;
    }
  }

  // The end-of-options marker should force everything after it to be positional.
  {
    auto argv = makeArgv({"app", "--flag", "--", "--notAFlag", "-x", "pos"});
    Args args;
    args.parse((int)argv.size(), argv.data());

    if (!args.hasFlag("flag")) {
      std::cerr << "[test_args] expected --flag to be recognized\n";
      ++fails;
    }
    if (args.hasFlag("notAFlag") || args.hasFlag("x")) {
      std::cerr << "[test_args] expected tokens after -- to NOT be parsed as flags\n";
      ++fails;
    }
    const auto& pos = args.positional();
    if (pos.size() != 3 || pos[0] != "--notAFlag" || pos[1] != "-x" || pos[2] != "pos") {
      std::cerr << "[test_args] expected 3 positional args after --, got size=" << pos.size() << "\n";
      ++fails;
    }
  }

  // A single '-' token is commonly used as a stdin/stdout placeholder and should not be
  // treated as an option introducer for --key value pairs.
  {
    auto argv = makeArgv({"app", "--out", "-"});
    Args args;
    args.parse((int)argv.size(), argv.data());

    if (args.hasFlag("out")) {
      std::cerr << "[test_args] expected --out - to be parsed as a key/value, not a flag\n";
      ++fails;
    }

    const auto v = args.last("out");
    if (!v || *v != "-") {
      std::cerr << "[test_args] expected --out to have value '-'\n";
      ++fails;
    }
  }

  // Short options can take values if configured via setArity (e.g. -o file.txt).
  {
    auto argv = makeArgv({"app", "-o", "file.txt"});
    Args args;
    args.setArity("o", 1);
    args.parse((int)argv.size(), argv.data());

    const auto v = args.last("o");
    if (!v || *v != "file.txt") {
      std::cerr << "[test_args] expected -o file.txt to parse as key 'o'\n";
      ++fails;
    }
  }

  // Support a minimal -k=value form for single-letter options (useful for values that begin with '-').
  {
    auto argv = makeArgv({"app", "-o=-"});
    Args args;
    args.parse((int)argv.size(), argv.data());

    const auto v = args.last("o");
    if (!v || *v != "-") {
      std::cerr << "[test_args] expected -o=- to parse value '-'\n";
      ++fails;
    }
  }

  // Negative numbers that appear as positional args should stay positional (not become short flags).
  {
    auto argv = makeArgv({"app", "-1", "-0.25", "-.5"});
    Args args;
    args.parse((int)argv.size(), argv.data());
    const auto& pos = args.positional();
    if (pos.size() != 3 || pos[0] != "-1" || pos[1] != "-0.25" || pos[2] != "-.5") {
      std::cerr << "[test_args] expected negative positional args to remain positional\n";
      ++fails;
    }
  }

  // Short flag parsing should still work.
  {
    auto argv = makeArgv({"app", "-hv"});
    Args args;
    args.parse((int)argv.size(), argv.data());
    if (!args.hasFlag("h") || !args.hasFlag("v")) {
      std::cerr << "[test_args] expected -hv to set flags h and v\n";
      ++fails;
    }
  }

  if (fails == 0) std::cout << "[test_args] pass\n";
  return fails;
}
