#include "test_harness.h"

#include "stellar/core/Base64.h"

#include <string>

using namespace stellar;

int test_base64() {
  int failures = 0;

  auto roundTrip = [&](const std::string& in) {
    const std::string enc = core::base64Encode(in);
    std::string dec;
    const bool ok = core::base64Decode(enc, &dec);
    CHECK(ok);
    CHECK(dec == in);
  };

  roundTrip("");
  roundTrip("hello");
  roundTrip("hello world");
  roundTrip("line1\nline2\r\nline3");

  // UTF-8 payload should round-trip as bytes.
  roundTrip("Nebula: \xF0\x9F\x8C\x8C  Star: \xE2\x98\x85");

  // A few known encodings (spot checks).
  CHECK(core::base64Encode("f") == "Zg==");
  CHECK(core::base64Encode("fo") == "Zm8=");
  CHECK(core::base64Encode("foo") == "Zm9v");
  CHECK(core::base64Encode("foobar") == "Zm9vYmFy");

  // Decoder should reject obviously invalid input.
  {
    std::string out;
    CHECK(!core::base64Decode("@@@", &out));
  }

  return failures;
}
