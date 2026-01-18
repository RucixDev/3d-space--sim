#include "stellar/core/Base64.h"

#include <array>
#include <cctype>
#include <cstdint>

namespace stellar::core {

namespace {

static constexpr char kAlphabet[] =
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static std::array<int8_t, 256> makeDecodeTable() {
  std::array<int8_t, 256> t{};
  t.fill(-1);
  for (int i = 0; i < 64; ++i) {
    t[(unsigned char)kAlphabet[i]] = static_cast<int8_t>(i);
  }
  return t;
}

static const std::array<int8_t, 256> kDecodeTable = makeDecodeTable();

} // namespace

std::string base64Encode(std::string_view bytes) {
  if (bytes.empty()) return {};

  const std::size_t n = bytes.size();
  const unsigned char* in = reinterpret_cast<const unsigned char*>(bytes.data());

  std::string out;
  out.reserve(((n + 2) / 3) * 4);

  std::size_t i = 0;
  while (i + 3 <= n) {
    const uint32_t b0 = in[i++];
    const uint32_t b1 = in[i++];
    const uint32_t b2 = in[i++];

    out.push_back(kAlphabet[(b0 >> 2) & 0x3Fu]);
    out.push_back(kAlphabet[((b0 & 0x3u) << 4) | ((b1 >> 4) & 0xFu)]);
    out.push_back(kAlphabet[((b1 & 0xFu) << 2) | ((b2 >> 6) & 0x3u)]);
    out.push_back(kAlphabet[b2 & 0x3Fu]);
  }

  const std::size_t rem = n - i;
  if (rem == 1) {
    const uint32_t b0 = in[i++];
    out.push_back(kAlphabet[(b0 >> 2) & 0x3Fu]);
    out.push_back(kAlphabet[((b0 & 0x3u) << 4) & 0x3Fu]);
    out.push_back('=');
    out.push_back('=');
  } else if (rem == 2) {
    const uint32_t b0 = in[i++];
    const uint32_t b1 = in[i++];
    out.push_back(kAlphabet[(b0 >> 2) & 0x3Fu]);
    out.push_back(kAlphabet[((b0 & 0x3u) << 4) | ((b1 >> 4) & 0xFu)]);
    out.push_back(kAlphabet[((b1 & 0xFu) << 2) & 0x3Fu]);
    out.push_back('=');
  }

  return out;
}

bool base64Decode(std::string_view b64, std::string* outBytes) {
  if (!outBytes) return false;
  outBytes->clear();

  if (b64.empty()) return true;

  // For efficiency, reserve the maximum possible size.
  // (It can be slightly smaller after padding is applied.)
  outBytes->reserve((b64.size() / 4) * 3);

  int sextets[4] = {0, 0, 0, 0};
  int count = 0;
  int pad = 0;
  bool done = false;

  auto flushGroup = [&]() -> bool {
    // count is assumed to be 4.
    const uint32_t triple =
        ((uint32_t)sextets[0] << 18) |
        ((uint32_t)sextets[1] << 12) |
        ((uint32_t)sextets[2] << 6) |
        ((uint32_t)sextets[3]);

    outBytes->push_back((char)((triple >> 16) & 0xFFu));
    if (pad < 2) outBytes->push_back((char)((triple >> 8) & 0xFFu));
    if (pad < 1) outBytes->push_back((char)(triple & 0xFFu));

    if (pad != 0) done = true;
    count = 0;
    pad = 0;
    return true;
  };

  for (char ch : b64) {
    const unsigned char c = static_cast<unsigned char>(ch);

    if (std::isspace(c)) {
      continue;
    }

    if (done) {
      // Only allow trailing whitespace after a padded block.
      return false;
    }

    if (c == '=') {
      sextets[count++] = 0;
      ++pad;
      if (pad > 2) return false;
    } else {
      const int8_t v = kDecodeTable[c];
      if (v < 0) return false;
      sextets[count++] = (int)v;
    }

    if (count == 4) {
      if (!flushGroup()) return false;
    }
  }

  // Handle missing padding (non-RFC input, but common in the wild).
  if (count != 0) {
    if (count == 1) return false;

    if (count == 2) {
      // 2 sextets -> 12 bits -> 1 byte
      const uint32_t triple = ((uint32_t)sextets[0] << 18) | ((uint32_t)sextets[1] << 12);
      outBytes->push_back((char)((triple >> 16) & 0xFFu));
    } else if (count == 3) {
      // 3 sextets -> 18 bits -> 2 bytes
      const uint32_t triple =
          ((uint32_t)sextets[0] << 18) |
          ((uint32_t)sextets[1] << 12) |
          ((uint32_t)sextets[2] << 6);
      outBytes->push_back((char)((triple >> 16) & 0xFFu));
      outBytes->push_back((char)((triple >> 8) & 0xFFu));
    }
  }

  return true;
}

} // namespace stellar::core
