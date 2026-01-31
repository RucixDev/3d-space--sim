#include "stellar/sim/TimeTrialGhostCodec.h"

#include "stellar/core/Base64.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace stellar::sim {

namespace {

static constexpr char kMagic[4] = {'T', 'T', 'G', 'H'}; // Time Trial Ghost
static constexpr core::u8 kWireVersion = 1;

static constexpr std::size_t kFloatsPerSample = 8; // t + pos(3) + orient(4)
static constexpr std::size_t kBytesPerSample = kFloatsPerSample * 4;

static constexpr core::u32 kMaxSamples = 200000; // safety limit (~6.4MB raw)

inline void writeU8(std::string& out, core::u8 v) {
  out.push_back(static_cast<char>(v));
}

inline void writeU32LE(std::string& out, core::u32 v) {
  out.push_back(static_cast<char>(v & 0xFFu));
  out.push_back(static_cast<char>((v >> 8) & 0xFFu));
  out.push_back(static_cast<char>((v >> 16) & 0xFFu));
  out.push_back(static_cast<char>((v >> 24) & 0xFFu));
}

inline bool readU8(const std::string& in, std::size_t& idx, core::u8* out) {
  if (!out) return false;
  if (idx + 1 > in.size()) return false;
  *out = static_cast<core::u8>(static_cast<unsigned char>(in[idx]));
  idx += 1;
  return true;
}

inline bool readU32LE(const std::string& in, std::size_t& idx, core::u32* out) {
  if (!out) return false;
  if (idx + 4 > in.size()) return false;
  const auto b0 = static_cast<core::u32>(static_cast<unsigned char>(in[idx + 0]));
  const auto b1 = static_cast<core::u32>(static_cast<unsigned char>(in[idx + 1]));
  const auto b2 = static_cast<core::u32>(static_cast<unsigned char>(in[idx + 2]));
  const auto b3 = static_cast<core::u32>(static_cast<unsigned char>(in[idx + 3]));
  *out = (b0) | (b1 << 8) | (b2 << 16) | (b3 << 24);
  idx += 4;
  return true;
}

inline void writeF32LE(std::string& out, float f) {
  static_assert(sizeof(float) == 4);
  core::u32 bits = 0;
  std::memcpy(&bits, &f, 4);
  writeU32LE(out, bits);
}

inline bool readF32LE(const std::string& in, std::size_t& idx, float* out) {
  if (!out) return false;
  core::u32 bits = 0;
  if (!readU32LE(in, idx, &bits)) return false;
  float f = 0.0f;
  std::memcpy(&f, &bits, 4);
  *out = f;
  return true;
}

inline bool finiteF(float f) {
  return std::isfinite(static_cast<double>(f));
}

} // namespace

std::string encodeTimeTrialGhostSamplesB64(const std::vector<TimeTrialGhostSample>& samples) {
  if (samples.empty()) return {};

  const core::u32 n = static_cast<core::u32>(std::min<std::size_t>(samples.size(), kMaxSamples));

  std::string bytes;
  bytes.reserve(4 + 1 + 4 + (std::size_t)n * kBytesPerSample);

  bytes.append(kMagic, kMagic + 4);
  writeU8(bytes, kWireVersion);
  writeU32LE(bytes, n);

  for (core::u32 i = 0; i < n; ++i) {
    const auto& s = samples[(std::size_t)i];

    const float t = static_cast<float>(std::clamp(s.tSec, 0.0, 1e9));
    const float px = static_cast<float>(s.posKm.x);
    const float py = static_cast<float>(s.posKm.y);
    const float pz = static_cast<float>(s.posKm.z);

    math::Quatd q = s.orient;
    const double lenSq = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
    if (!(lenSq > 1e-24) || !std::isfinite(lenSq)) {
      q = math::Quatd::identity();
    } else {
      q = q.normalized();
    }

    writeF32LE(bytes, t);
    writeF32LE(bytes, finiteF(px) ? px : 0.0f);
    writeF32LE(bytes, finiteF(py) ? py : 0.0f);
    writeF32LE(bytes, finiteF(pz) ? pz : 0.0f);

    writeF32LE(bytes, static_cast<float>(q.w));
    writeF32LE(bytes, static_cast<float>(q.x));
    writeF32LE(bytes, static_cast<float>(q.y));
    writeF32LE(bytes, static_cast<float>(q.z));
  }

  return core::base64Encode(bytes);
}

bool decodeTimeTrialGhostSamplesB64(std::string_view b64,
                                   std::vector<TimeTrialGhostSample>* outSamples,
                                   std::string* err) {
  if (!outSamples) {
    if (err) *err = "null outSamples";
    return false;
  }

  outSamples->clear();

  if (b64.empty() || b64 == "~") {
    return true;
  }

  std::string bytes;
  if (!core::base64Decode(b64, &bytes)) {
    if (err) *err = "invalid base64";
    return false;
  }

  if (bytes.size() < 9) {
    if (err) *err = "ghost blob too small";
    return false;
  }

  if (!(bytes[0] == kMagic[0] && bytes[1] == kMagic[1] && bytes[2] == kMagic[2] && bytes[3] == kMagic[3])) {
    if (err) *err = "bad ghost magic";
    return false;
  }

  std::size_t idx = 4;
  core::u8 ver = 0;
  if (!readU8(bytes, idx, &ver)) {
    if (err) *err = "truncated ghost header";
    return false;
  }

  if (ver != kWireVersion) {
    if (err) *err = "unsupported ghost version";
    return false;
  }

  core::u32 n = 0;
  if (!readU32LE(bytes, idx, &n)) {
    if (err) *err = "truncated ghost count";
    return false;
  }

  if (n > kMaxSamples) {
    if (err) *err = "ghost too large";
    return false;
  }

  const std::size_t required = 4 + 1 + 4 + (std::size_t)n * kBytesPerSample;
  if (bytes.size() < required) {
    if (err) *err = "truncated ghost data";
    return false;
  }

  outSamples->reserve(n);

  for (core::u32 i = 0; i < n; ++i) {
    float t = 0.0f;
    float px = 0.0f, py = 0.0f, pz = 0.0f;
    float qw = 1.0f, qx = 0.0f, qy = 0.0f, qz = 0.0f;

    if (!readF32LE(bytes, idx, &t) ||
        !readF32LE(bytes, idx, &px) ||
        !readF32LE(bytes, idx, &py) ||
        !readF32LE(bytes, idx, &pz) ||
        !readF32LE(bytes, idx, &qw) ||
        !readF32LE(bytes, idx, &qx) ||
        !readF32LE(bytes, idx, &qy) ||
        !readF32LE(bytes, idx, &qz)) {
      outSamples->clear();
      if (err) *err = "truncated ghost sample";
      return false;
    }

    TimeTrialGhostSample s{};
    s.tSec = std::isfinite((double)t) ? std::max(0.0, (double)t) : 0.0;
    s.posKm = {
        std::isfinite((double)px) ? (double)px : 0.0,
        std::isfinite((double)py) ? (double)py : 0.0,
        std::isfinite((double)pz) ? (double)pz : 0.0,
    };

    math::Quatd q{(double)qw, (double)qx, (double)qy, (double)qz};
    const double lenSq = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
    if (!(lenSq > 1e-24) || !std::isfinite(lenSq)) {
      s.orient = math::Quatd::identity();
    } else {
      s.orient = q.normalized();
    }

    outSamples->push_back(s);
  }

  return true;
}

} // namespace stellar::sim
