#include "test_harness.h"

#include "stellar/ui/Format.h"

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

#include <cstdint>
#include <limits>
#include <string>
#include <filesystem>

int test_ui_json_writer() {
  int failures = 0;

  std::string out;

  // 1) Literal nullptr should serialize as JSON null.
  {
    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("a", nullptr);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"a\": null\n}"));
  }

  // 2) A typed c-string pointer that is null should also serialize as JSON null.
  // (This is the common case when serializing optional const char* fields.)
  {
    const char* s = nullptr;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("s", s);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"s\": null\n}"));
  }

  // 3) Non-null c-string should serialize as a JSON string.
  {
    const char* s = "hi";

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("s", s);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"s\": \"hi\"\n}"));
  }



  // 4) std::optional<T> should serialize to either null or the underlying value.
  {
    std::optional<int> v;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("v", v);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"v\": null\n}"));
  }

  {
    std::optional<int> v = 7;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("v", v);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"v\": 7\n}"));
  }

  {
    std::optional<const char*> s;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("s", s);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"s\": null\n}"));
  }

  {
    std::optional<const char*> s = "yo";

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("s", s);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"s\": \"yo\"\n}"));
  }



  // 5) std::vector<T> should serialize as a JSON array.
  {
    std::vector<int> a{1, 2, 3};

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("a", a);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "a": [
    1,
    2,
    3
  ]
})"));
  }

  

  // 5b) std::array<T,N> should serialize as a JSON array.
  {
    std::array<int, 3> a{{1, 2, 3}};

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("a", a);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "a": [
    1,
    2,
    3
  ]
})"));
  }

  // 5c) std::span<T> should serialize as a JSON array.
  {
    std::array<int, 3> backing{{1, 2, 3}};
    std::span<int> s(backing.data(), backing.size());

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("s", s);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "s": [
    1,
    2,
    3
  ]
})"));
  }

  // 5d) std::set<T> should serialize as a JSON array (naturally sorted).
  {
    std::set<int> s{3, 1, 2};

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("s", s);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "s": [
    1,
    2,
    3
  ]
})"));
  }

  // 5e) std::unordered_set<T> should serialize deterministically when elements are comparable.
  {
    std::unordered_set<int> s;
    s.insert(3);
    s.insert(1);
    s.insert(2);

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("s", s);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "s": [
    1,
    2,
    3
  ]
})"));
  }


  // 6) std::optional<std::vector<T>> should serialize as null or an array.
  {
    std::optional<std::vector<int>> v;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("v", v);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"v\": null\n}"));
  }

  {
    std::optional<std::vector<int>> v = std::vector<int>{4, 5};

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("v", v);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "v": [
    4,
    5
  ]
})"));
  }


  // 7) std::map<std::string, T> should serialize as a JSON object.
  {
    std::map<std::string, int> m{{"b", 2}, {"a", 1}};

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("m", m);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "m": {
    "a": 1,
    "b": 2
  }
})"));
  }

  // 8) std::unordered_map<std::string, T> should serialize deterministically (sorted by key).
  {
    std::unordered_map<std::string, int> m;
    m["b"] = 2;
    m["a"] = 1;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("m", m);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "m": {
    "a": 1,
    "b": 2
  }
})"));
  }

  // 8b) std::pair and std::tuple should serialize as JSON arrays.
  {
    const std::pair<int, const char*> p{1, "a"};

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("p", p);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "p": [
    1,
    "a"
  ]
})"));
  }

  {
    const auto t = std::make_tuple(1, std::string("x"), false);

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("t", t);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "t": [
    1,
    "x",
    false
  ]
})"));
  }

  {
    const std::tuple<> empty;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("t", empty);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"t\": []\n}"));
  }

  // 8c) std::variant should serialize as the held alternative.
  {
    std::variant<std::monostate, int, std::string> v;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("v", v);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"v\": null\n}"));
  }

  {
    std::variant<std::monostate, int, std::string> v = 3;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("v", v);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"v\": 3\n}"));
  }

  {
    std::variant<std::monostate, int, std::string> v = std::string("hi");

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("v", v);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"v\": \"hi\"\n}"));
  }


  // 9) 8-bit integer types should serialize as JSON numbers (not characters).
  // This is especially important for uint8_t/int8_t (aliases of unsigned/signed char).
  {
    std::uint8_t u8 = 65;
    std::int8_t i8 = -5;
    unsigned char uc = 120;
    signed char sc = -7;
    char ch = 42;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("u8", u8);
    w.member("i8", i8);
    w.member("uc", uc);
    w.member("sc", sc);
    w.member("ch", ch);
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "u8": 65,
  "i8": -5,
  "uc": 120,
  "sc": -7,
  "ch": 42
})"));
  }





  // 10) Floating point values should serialize as JSON numbers.
  {
    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("pi", 3.5);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"pi\": 3.5\n}"));
  }

  // 11) Non-finite floating point values should serialize as JSON null (JSON has no NaN/Inf).
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("nan", nan);
    w.member("inf", inf);
    w.member("ninf", -inf);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"nan\": null,\n  \"inf\": null,\n  \"ninf\": null\n}"));
  }


  // 12) memberIf(...) should omit absent values instead of writing null.
  {
    std::optional<int> v;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.memberIf("v", v);
    w.endObject();

    CHECK_EQ(out, std::string("{}"));
  }

  {
    std::optional<int> v = 7;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.memberIf("v", v);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"v\": 7\n}"));
  }

  {
    const char* s = nullptr;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.memberIf("s", s);
    w.endObject();

    CHECK_EQ(out, std::string("{}"));
  }

  {
    const char* s = "hi";

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.memberIf("s", s);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"s\": \"hi\"\n}"));
  }

  {
    int x = 5;
    int* p = &x;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.memberIf("x", p);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"x\": 5\n}"));
  }

  {
    int* p = nullptr;

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.memberIf("x", p);
    w.endObject();

    CHECK_EQ(out, std::string("{}"));
  }

  // 13) RAII scopes should close objects/arrays automatically.
  {
    stellar::ui::JsonWriter w(out);
    {
      auto obj = w.scopedObject();
      (void)obj;
      w.member("a", 1);
    }

    CHECK_EQ(out, std::string("{\n  \"a\": 1\n}"));
  }

  {
    stellar::ui::JsonWriter w(out);
    {
      auto arr = w.scopedArray();
      (void)arr;
      w.value(1);
      w.value(2);
    }

    CHECK_EQ(out, std::string("[\n  1,\n  2\n]"));
  }

  {
    stellar::ui::JsonWriter w(out);
    w.beginObject();
    {
      auto obj = w.scopedObject("inner");
      (void)obj;
      w.member("x", 1);
    }
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "inner": {
    "x": 1
  }
})"));
  }

  {
    stellar::ui::JsonWriter w(out);
    w.beginObject();
    {
      auto arr = w.scopedArray("a");
      (void)arr;
      w.value(1);
      w.value(2);
    }
    w.endObject();

    CHECK_EQ(out, std::string(R"({
  "a": [
    1,
    2
  ]
})"));
  }



  // N) std::filesystem::path should serialize as a JSON string (UTF-8), using
  // a stable generic representation with '/' separators.
  {
    const std::filesystem::path p = std::filesystem::path("foo") / "bar";

    stellar::ui::JsonWriter w(out);
    w.beginObject();
    w.member("p", p);
    w.endObject();

    CHECK_EQ(out, std::string("{\n  \"p\": \"foo/bar\"\n}"));
  }

  return failures;
}
