#include "stellar/core/CVar.h"

#include "test_harness.h"

#include <cstdio>
#include <string>

int test_cvars() {
  int failures = 0;

  using stellar::core::CVarRegistry;
  using stellar::core::CVarType;

  // ---- Define + typed get/set ----
  {
    CVarRegistry r;
    CHECK(r.defineBool("a.bool", true, stellar::core::CVar_Archive, "test") != nullptr);
    CHECK(r.defineInt("a.int", 42, stellar::core::CVar_None) != nullptr);
    CHECK(r.defineFloat("a.float", 1.5, stellar::core::CVar_None) != nullptr);
    CHECK(r.defineString("a.str", "hello") != nullptr);

    CHECK(r.getBool("a.bool", false) == true);
    CHECK(r.getInt("a.int", 0) == 42);
    CHECK(r.getFloat("a.float", 0.0) > 1.4);
    CHECK(r.getString("a.str", "") == "hello");

    std::string err;
    CHECK(r.setFromString("a.bool", "false", &err));
    CHECK(r.getBool("a.bool", true) == false);

    CHECK(r.setFromString("a.int", "-7", &err));
    CHECK(r.getInt("a.int", 0) == -7);

    CHECK(r.setFromString("a.float", "2.25", &err));
    CHECK(r.getFloat("a.float", 0.0) > 2.2);

    CHECK(r.setFromString("a.str", "\"hi there\"", &err));
    CHECK(r.getString("a.str", "") == "hi there");
  }

  // ---- Pending assignment (load before define) ----
  {
    const std::string path = "stellar_test_cvars_pending.cfg";
    {
      std::FILE* f = std::fopen(path.c_str(), "wb");
      CHECK(f != nullptr);
      if (f) {
        std::fputs("# test\n", f);
        std::fputs("pending.int = 123\n", f);
        std::fputs("pending.str = \"hello world\"\n", f);
        std::fclose(f);
      }
    }

    CVarRegistry r;
    std::string err;
    CHECK(r.loadFile(path, &err));
    CHECK(r.hasPending("pending.int"));
    CHECK(r.hasPending("pending.str"));

    // Defining should apply pending assignments automatically.
    CHECK(r.defineInt("pending.int", 0) != nullptr);
    CHECK(r.defineString("pending.str", "x") != nullptr);

    CHECK(r.getInt("pending.int", 0) == 123);
    CHECK(r.getString("pending.str", "") == "hello world");

    std::remove(path.c_str());
  }

  // ---- Save + reload round-trip ----
  {
    const std::string path = "stellar_test_cvars_roundtrip.cfg";

    CVarRegistry a;
    CHECK(a.defineBool("rt.b", true, stellar::core::CVar_Archive) != nullptr);
    CHECK(a.defineInt("rt.i", 1, stellar::core::CVar_Archive) != nullptr);
    CHECK(a.defineFloat("rt.f", 1.0, stellar::core::CVar_Archive) != nullptr);
    CHECK(a.defineString("rt.s", "hello world", stellar::core::CVar_Archive) != nullptr);

    std::string err;
    CHECK(a.setFromString("rt.b", "0", &err));
    CHECK(a.setFromString("rt.i", "99", &err));
    CHECK(a.setFromString("rt.f", "3.5", &err));
    CHECK(a.setFromString("rt.s", "\"with spaces\"", &err));

    CHECK(a.saveFile(path, &err));

    // Load into a fresh registry, but load BEFORE define to exercise pending.
    CVarRegistry b;
    CHECK(b.loadFile(path, &err));

    CHECK(b.defineBool("rt.b", true, stellar::core::CVar_Archive) != nullptr);
    CHECK(b.defineInt("rt.i", 0, stellar::core::CVar_Archive) != nullptr);
    CHECK(b.defineFloat("rt.f", 0.0, stellar::core::CVar_Archive) != nullptr);
    CHECK(b.defineString("rt.s", "", stellar::core::CVar_Archive) != nullptr);

    CHECK(b.getBool("rt.b", true) == false);
    CHECK(b.getInt("rt.i", 0) == 99);
    CHECK(b.getFloat("rt.f", 0.0) > 3.4);
    CHECK(b.getString("rt.s", "") == "with spaces");

    std::remove(path.c_str());
  }

  return failures;
}
