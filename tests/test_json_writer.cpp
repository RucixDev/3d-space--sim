// Minimal unit tests for stellar::core::JsonWriter.
//
// The project uses a simple home-grown test harness (tests/test_harness.h)
// where each translation unit exports an `int test_<name>()` entry point.
//
// This file intentionally avoids external test frameworks.

#include "test_harness.h"

#include <stellar/core/JsonWriter.h>

#include <string>

static int failures = 0;

static void expectEq(const std::string& a, const std::string& b, const char* msg) {
  if (a != b) {
    std::cerr << "FAILED: " << msg << "\n  expected: " << b << "\n  got     : " << a
              << "\n";
    ++failures;
  }
}

int test_json_writer() {
  failures = 0;

  // 1) Object with a few primitive fields.
  {
    stellar::core::JsonWriter w;
    w.beginObject();
    w.key("answer");
    w.value(42);
    w.key("pi");
    w.value(3.5);
    w.key("ok");
    w.value(true);
    w.endObject();
    expectEq(w.takeString(),
             std::string("{\"answer\":42,\"pi\":3.5,\"ok\":true}"),
             "object primitives");
  }

  // 2) Nested array + object.
  {
    stellar::core::JsonWriter w;
    w.beginObject();
    w.key("arr");
    w.beginArray();
    w.value(1);
    w.value(2);
    w.beginObject();
    w.key("x");
    w.value("y");
    w.endObject();
    w.endArray();
    w.endObject();
    expectEq(w.takeString(), std::string("{\"arr\":[1,2,{\"x\":\"y\"}]}"),
             "nested structures");
  }

  // 3) String escaping.
  {
    stellar::core::JsonWriter w;
    w.beginObject();
    w.key("s");
    w.value("a\n\t\\\"b");
    w.endObject();
    expectEq(w.takeString(),
             std::string("{\"s\":\"a\\n\\t\\\\\\\"b\"}"),
             "string escaping");
  }

  // 4) Null.
  {
    stellar::core::JsonWriter w;
    w.beginArray();
    w.nullValue();
    w.endArray();
    expectEq(w.takeString(), std::string("[null]"), "null literal");
  }

  if (failures == 0) {
    std::cout << "All JsonWriter tests passed\n";
  }
  return failures;
}
