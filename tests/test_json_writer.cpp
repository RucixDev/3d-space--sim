// Minimal unit tests for stellar::core::JsonWriter.
//
// The project uses a simple home-grown test harness (tests/test_harness.h)
// where each translation unit exports an `int test_<name>()` entry point.
//
// This file intentionally avoids external test frameworks.

#include "test_harness.h"

#include <stellar/core/JsonWriter.h>

#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <optional>
#include <vector>
#include <unordered_set>
#include <set>
#include <span>
#include <array>
#include <map>
#include <unordered_map>
#include <tuple>
#include <variant>

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

  // 1b) Non-finite numbers (NaN/Inf) should serialize as JSON null (JSON has no such numeric literals).
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("nan", nan);
    w.field("inf", inf);
    w.field("ninf", -inf);
    w.endObject();

    expectEq(w.takeString(),
             std::string("{\"nan\":null,\"inf\":null,\"ninf\":null}"),
             "non-finite -> null");
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

  // 5) nullptr / null c-strings should serialize as JSON null (not empty string).
  {
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("a", nullptr);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"a\":null}"), "nullptr literal");
  }

  {
    const char* s = nullptr;
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("s", s);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"s\":null}"), "nullptr c-string");
  }

  // 6) std::optional<T> should serialize to either null or the underlying value.
  {
    std::optional<int> v;
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("v", v);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"v\":null}"), "optional empty -> null");
  }

  {
    std::optional<int> v = 7;
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("v", v);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"v\":7}"), "optional int");
  }

  {
    std::optional<const char*> s;
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("s", s);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"s\":null}"), "optional c-string empty");
  }

  {
    std::optional<const char*> s = "yo";
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("s", s);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"s\":\"yo\"}"), "optional c-string value");
  }

  // 7) std::vector<T> should serialize as a JSON array.
  {
    std::vector<int> v{1, 2, 3};
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("v", v);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"v\":[1,2,3]}"), "vector int");
  }

  // 8) std::optional<std::vector<T>> should serialize as null or an array.
  {
    std::optional<std::vector<int>> v;
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("v", v);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"v\":null}"), "optional<vector> empty");
  }

  {
    std::optional<std::vector<int>> v = std::vector<int>{4, 5};
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("v", v);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"v\":[4,5]}"), "optional<vector> value");
  }


  // 9) std::map<std::string, T> should serialize as a JSON object.
  {
    std::map<std::string, int> m{{"b", 2}, {"a", 1}};
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("m", m);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"m\":{\"a\":1,\"b\":2}}"), "map<string,int>");
  }

  // 10) std::unordered_map<std::string, T> should serialize deterministically (sorted by key).
  {
    std::unordered_map<std::string, int> m;
    m["b"] = 2;
    m["a"] = 1;
    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("m", m);
    w.endObject();
    expectEq(w.takeString(), std::string("{\"m\":{\"a\":1,\"b\":2}}"), "unordered_map deterministic");
  }

  

  // 10b) std::array<T,N> should serialize as a JSON array.
  {
    std::array<int, 3> a{{1, 2, 3}};

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("a", a);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"a\":[1,2,3]}"), "array<int,3>");
  }

  // 10c) std::span<T> should serialize as a JSON array.
  {
    std::array<int, 3> backing{{1, 2, 3}};
    std::span<int> s(backing.data(), backing.size());

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("s", s);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"s\":[1,2,3]}"), "span<int>");
  }

  // 10d) std::set<T> should serialize as a JSON array (naturally sorted).
  {
    std::set<int> s{3, 1, 2};

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("s", s);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"s\":[1,2,3]}"), "set<int>");
  }

  // 10e) std::unordered_set<T> should serialize deterministically when elements are comparable.
  {
    std::unordered_set<int> s;
    s.insert(3);
    s.insert(1);
    s.insert(2);

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("s", s);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"s\":[1,2,3]}"), "unordered_set<int> deterministic");
  }


  // 11) 8-bit integer types should serialize as JSON numbers (not characters).
  // This is especially important for uint8_t/int8_t (aliases of unsigned/signed char).
  {
    std::uint8_t u8 = 65;
    std::int8_t i8 = -5;
    unsigned char uc = 120;
    signed char sc = -7;
    char ch = 42;

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("u8", u8);
    w.field("i8", i8);
    w.field("uc", uc);
    w.field("sc", sc);
    w.field("ch", ch);
    w.endObject();

    expectEq(w.takeString(),
             std::string("{\"u8\":65,\"i8\":-5,\"uc\":120,\"sc\":-7,\"ch\":42}"),
             "8-bit integers numeric");
  }

  // 12) std::pair and std::tuple should serialize as JSON arrays.
  {
    std::pair<int, const char*> p{1, "a"};

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("p", p);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"p\":[1,\"a\"]}"), "pair<int,const char*>");
  }

  {
    const auto t = std::make_tuple(1, std::string("x"), false);

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("t", t);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"t\":[1,\"x\",false]}"), "tuple<int,string,bool>");
  }

  {
    const std::tuple<> empty;

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("t", empty);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"t\":[]}"), "tuple<> empty");
  }

  // 13) std::variant should serialize as the held alternative.
  {
    std::variant<std::monostate, int, std::string> v;

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("v", v);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"v\":null}"), "variant monostate -> null");
  }

  {
    std::variant<std::monostate, int, std::string> v = 3;

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("v", v);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"v\":3}"), "variant int");
  }

  {
    std::variant<std::monostate, int, std::string> v = std::string("hi");

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("v", v);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"v\":\"hi\"}"), "variant string");
  }


  // 14) fieldIf(...) should omit absent values instead of writing null.
  {
    std::optional<int> v;

    stellar::core::JsonWriter w;
    w.beginObject();
    w.fieldIf("v", v);
    w.endObject();

    expectEq(w.takeString(), std::string("{}"), "fieldIf optional empty omitted");
  }

  {
    std::optional<int> v = 7;

    stellar::core::JsonWriter w;
    w.beginObject();
    w.fieldIf("v", v);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"v\":7}"), "fieldIf optional value");
  }

  {
    const char* s = nullptr;

    stellar::core::JsonWriter w;
    w.beginObject();
    w.fieldIf("s", s);
    w.endObject();

    expectEq(w.takeString(), std::string("{}"), "fieldIf c-string null omitted");
  }

  {
    const char* s = "hi";

    stellar::core::JsonWriter w;
    w.beginObject();
    w.fieldIf("s", s);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"s\":\"hi\"}"), "fieldIf c-string");
  }

  // 15) RAII scopes should close objects/arrays automatically.
  {
    stellar::core::JsonWriter w;
    {
      auto obj = w.scopedObject();
      (void)obj;
      w.field("a", 1);
    }
    expectEq(w.takeString(), std::string("{\"a\":1}"), "scopedObject root");
  }

  {
    stellar::core::JsonWriter w;
    {
      auto arr = w.scopedArray();
      (void)arr;
      w.value(1);
      w.value(2);
    }
    expectEq(w.takeString(), std::string("[1,2]"), "scopedArray root");
  }

  {
    stellar::core::JsonWriter w;
    w.beginObject();
    {
      auto obj = w.scopedObject("inner");
      (void)obj;
      w.field("x", 1);
    }
    w.endObject();
    expectEq(w.takeString(), std::string("{\"inner\":{\"x\":1}}"), "scopedObject key");
  }

  {
    stellar::core::JsonWriter w;
    w.beginObject();
    {
      auto arr = w.scopedArray("a");
      (void)arr;
      w.value(1);
      w.value(2);
    }
    w.endObject();
    expectEq(w.takeString(), std::string("{\"a\":[1,2]}"), "scopedArray key");
  }

  // 16) Pointer fieldIf overload should dereference non-null pointers.
  {
    int x = 5;
    int* p = &x;

    stellar::core::JsonWriter w;
    w.beginObject();
    w.fieldIf("x", p);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"x\":5}"), "fieldIf pointer");
  }

  {
    int* p = nullptr;

    stellar::core::JsonWriter w;
    w.beginObject();
    w.fieldIf("x", p);
    w.endObject();

    expectEq(w.takeString(), std::string("{}"), "fieldIf pointer null omitted");
  }







  // 17) std::filesystem::path should serialize as a JSON string (UTF-8), using
  // a stable generic representation with '/' separators.
  {
    const std::filesystem::path p = std::filesystem::path("foo") / "bar";

    stellar::core::JsonWriter w;
    w.beginObject();
    w.field("p", p);
    w.endObject();

    expectEq(w.takeString(), std::string("{\"p\":\"foo/bar\"}"), "filesystem path");
  }
  if (failures == 0) {
    std::cout << "All JsonWriter tests passed\n";
  }
  return failures;
}
