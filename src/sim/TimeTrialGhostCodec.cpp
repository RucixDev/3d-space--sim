#include "stellar/sim/TimeTrialGhostCodec.h"

#include "stellar/core/Base64.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace stellar::sim {

namespace {

static constexpr char kMagic[4] = {'T', 'T', 'G', 'H'}; // Time Trial Ghost
static constexpr core::u8 kWireVersionV1 = 1;
static constexpr core::u8 kWireVersionV2 = 2;
static constexpr core::u8 kWireVersionCurrent = kWireVersionV2;

static constexpr std::size_t kFloatsPerSampleV1 = 8; // t + pos(3) + orient(4)
static constexpr std::size_t kBytesPerSampleV1 = kFloatsPerSampleV1 * 4;

// v2 (delta+quantized) record sizes (bytes):
//   Keyframe: type(1) + full sample (8 f32 = 32) = 33
//   Delta:    type(1) + t f32(4) + dpos 3*i16(6) + quat idx u8(1) + quat 3*i16(6) = 18
static constexpr std::size_t kBytesPerV2Keyframe = 1 + 32;
static constexpr std::size_t kBytesPerV2Delta = 1 + 4 + 6 + 1 + 6;

static constexpr core::u32 kMaxSamples = 200000; // safety limit (~6.4MB raw)

// v2 position delta quantization:
//  - unit: km
//  - quantum: 1e-4 km = 10cm
//  - max delta per sample: 32767 * quantum ~= 3.2767 km
static constexpr double kPosQuantumKm = 1e-4;
static constexpr double kInvPosQuantumKm = 1.0 / kPosQuantumKm;
static constexpr core::i32 kPosDeltaMaxQ = 32767;

// v2 quaternion compression: "smallest three".
// For a unit quaternion, if you drop the largest-magnitude component, the remaining three
// are in [-1/sqrt(2), +1/sqrt(2)]. We scale by sqrt(2) to map to [-1, +1] for quantization.
static constexpr double kInvSqrt2 = 0.70710678118654752440084436210485;
static constexpr double kSqrt2 = 1.4142135623730950488016887242097;
static constexpr double kQuatQuantScale = 32767.0;

enum class V2RecordType : core::u8 {
  Delta = 0,
  Keyframe = 1,
};

inline void writeU8(std::string& out, core::u8 v) {
  out.push_back(static_cast<char>(v));
}

inline void writeU16LE(std::string& out, core::u16 v) {
  out.push_back(static_cast<char>(v & 0xFFu));
  out.push_back(static_cast<char>((v >> 8) & 0xFFu));
}

inline void writeU32LE(std::string& out, core::u32 v) {
  out.push_back(static_cast<char>(v & 0xFFu));
  out.push_back(static_cast<char>((v >> 8) & 0xFFu));
  out.push_back(static_cast<char>((v >> 16) & 0xFFu));
  out.push_back(static_cast<char>((v >> 24) & 0xFFu));
}

inline void writeI16LE(std::string& out, core::i16 v) {
  writeU16LE(out, static_cast<core::u16>(v));
}

inline bool readU8(const std::string& in, std::size_t& idx, core::u8* out) {
  if (!out) return false;
  if (idx + 1 > in.size()) return false;
  *out = static_cast<core::u8>(static_cast<unsigned char>(in[idx]));
  idx += 1;
  return true;
}

inline bool readU16LE(const std::string& in, std::size_t& idx, core::u16* out) {
  if (!out) return false;
  if (idx + 2 > in.size()) return false;
  const auto b0 = static_cast<core::u16>(static_cast<unsigned char>(in[idx + 0]));
  const auto b1 = static_cast<core::u16>(static_cast<unsigned char>(in[idx + 1]));
  *out = static_cast<core::u16>((b0) | (b1 << 8));
  idx += 2;
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

inline bool readI16LE(const std::string& in, std::size_t& idx, core::i16* out) {
  if (!out) return false;
  core::u16 u = 0;
  if (!readU16LE(in, idx, &u)) return false;
  *out = static_cast<core::i16>(u);
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

inline bool finiteD(double d) {
  return std::isfinite(d);
}

inline math::Vec3d sanitizePos(const math::Vec3d& p) {
  return {
      finiteD(p.x) ? p.x : 0.0,
      finiteD(p.y) ? p.y : 0.0,
      finiteD(p.z) ? p.z : 0.0,
  };
}

inline math::Quatd sanitizeQuat(const math::Quatd& qIn) {
  math::Quatd q = qIn;
  const double lenSq = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
  if (!(lenSq > 1e-24) || !std::isfinite(lenSq)) {
    return math::Quatd::identity();
  }
  return q.normalized();
}

inline float safeF32(double d) {
  if (!finiteD(d)) return 0.0f;
  const double clamped = std::clamp(d, -1e30, 1e30);
  return static_cast<float>(clamped);
}

struct PackedQuatSmallestThree {
  core::u8 largestIndex{0}; // 0=w,1=x,2=y,3=z
  core::i16 a{0};
  core::i16 b{0};
  core::i16 c{0};
};

inline core::i16 quantizeQuatComponentSmallestThree(double v) {
  // v is expected in [-1/sqrt(2), +1/sqrt(2)]. Scale by sqrt(2) -> [-1, +1].
  if (!finiteD(v)) v = 0.0;
  double x = std::clamp(v * kSqrt2, -1.0, 1.0);
  long long q = std::llround(x * kQuatQuantScale);
  if (q < -32767) q = -32767;
  if (q > 32767) q = 32767;
  return static_cast<core::i16>(q);
}

inline double dequantizeQuatComponentSmallestThree(core::i16 q) {
  const double x = std::clamp(static_cast<double>(q) / kQuatQuantScale, -1.0, 1.0);
  return x * kInvSqrt2;
}

inline PackedQuatSmallestThree packQuatSmallestThree(const math::Quatd& qIn) {
  math::Quatd q = sanitizeQuat(qIn);
  std::array<double, 4> c = {q.w, q.x, q.y, q.z};

  int idx = 0;
  double maxAbs = std::abs(c[0]);
  for (int i = 1; i < 4; ++i) {
    const double a = std::abs(c[i]);
    if (a > maxAbs) {
      maxAbs = a;
      idx = i;
    }
  }

  // Canonicalize sign: ensure the largest component is positive.
  if (c[idx] < 0.0) {
    for (double& v : c) v = -v;
  }

  double a = 0.0, b = 0.0, cc = 0.0;
  switch (idx) {
    case 0: a = c[1]; b = c[2]; cc = c[3]; break; // drop w
    case 1: a = c[0]; b = c[2]; cc = c[3]; break; // drop x
    case 2: a = c[0]; b = c[1]; cc = c[3]; break; // drop y
    default: a = c[0]; b = c[1]; cc = c[2]; break; // drop z
  }

  PackedQuatSmallestThree out;
  out.largestIndex = static_cast<core::u8>(idx);
  out.a = quantizeQuatComponentSmallestThree(a);
  out.b = quantizeQuatComponentSmallestThree(b);
  out.c = quantizeQuatComponentSmallestThree(cc);
  return out;
}

inline math::Quatd unpackQuatSmallestThree(core::u8 largestIndex, core::i16 aQ, core::i16 bQ, core::i16 cQ) {
  const int idx = static_cast<int>(largestIndex);
  if (idx < 0 || idx > 3) {
    return math::Quatd::identity();
  }

  const double a = dequantizeQuatComponentSmallestThree(aQ);
  const double b = dequantizeQuatComponentSmallestThree(bQ);
  const double c = dequantizeQuatComponentSmallestThree(cQ);

  const double sumSq = a * a + b * b + c * c;
  const double largest = std::sqrt(std::max(0.0, 1.0 - sumSq));

  std::array<double, 4> out = {0.0, 0.0, 0.0, 0.0};
  switch (idx) {
    case 0: out[0] = largest; out[1] = a; out[2] = b; out[3] = c; break;
    case 1: out[1] = largest; out[0] = a; out[2] = b; out[3] = c; break;
    case 2: out[2] = largest; out[0] = a; out[1] = b; out[3] = c; break;
    default: out[3] = largest; out[0] = a; out[1] = b; out[2] = c; break;
  }

  return math::Quatd{out[0], out[1], out[2], out[3]}.normalized();
}

inline void writeFullSampleF32(std::string& out,
                               float t,
                               const math::Vec3d& posKm,
                               const math::Quatd& orient) {
  const float px = safeF32(posKm.x);
  const float py = safeF32(posKm.y);
  const float pz = safeF32(posKm.z);

  math::Quatd q = sanitizeQuat(orient);

  writeF32LE(out, t);
  writeF32LE(out, finiteF(px) ? px : 0.0f);
  writeF32LE(out, finiteF(py) ? py : 0.0f);
  writeF32LE(out, finiteF(pz) ? pz : 0.0f);
  writeF32LE(out, safeF32(q.w));
  writeF32LE(out, safeF32(q.x));
  writeF32LE(out, safeF32(q.y));
  writeF32LE(out, safeF32(q.z));
}

inline bool decodeV1(const std::string& bytes,
                     std::size_t idx,
                     core::u32 n,
                     std::vector<TimeTrialGhostSample>* outSamples,
                     std::string* err) {
  if (!outSamples) {
    if (err) *err = "null outSamples";
    return false;
  }

  const std::size_t required = 4 + 1 + 4 + (std::size_t)n * kBytesPerSampleV1;
  if (bytes.size() < required) {
    if (err) *err = "truncated ghost data";
    return false;
  }

  outSamples->reserve(n);
  double prevT = 0.0;
  bool hasPrevT = false;

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

    const double tSec = std::isfinite((double)t) ? std::max(0.0, (double)t) : 0.0;
    if (hasPrevT && tSec + 1e-9 < prevT) {
      outSamples->clear();
      if (err) *err = "non-monotonic ghost time";
      return false;
    }
    prevT = tSec;
    hasPrevT = true;

    TimeTrialGhostSample s{};
    s.tSec = tSec;
    s.posKm = {
        std::isfinite((double)px) ? (double)px : 0.0,
        std::isfinite((double)py) ? (double)py : 0.0,
        std::isfinite((double)pz) ? (double)pz : 0.0,
    };

    math::Quatd q{(double)qw, (double)qx, (double)qy, (double)qz};
    s.orient = sanitizeQuat(q);

    outSamples->push_back(s);
  }

  return true;
}

inline bool decodeV2(const std::string& bytes,
                     std::size_t idx,
                     core::u32 n,
                     std::vector<TimeTrialGhostSample>* outSamples,
                     std::string* err) {
  if (!outSamples) {
    if (err) *err = "null outSamples";
    return false;
  }

  outSamples->reserve(n);

  math::Vec3d prevPosKm{0.0, 0.0, 0.0};
  bool hasPrevPos = false;
  double prevT = 0.0;
  bool hasPrevT = false;

  for (core::u32 i = 0; i < n; ++i) {
    core::u8 typeU8 = 0;
    if (!readU8(bytes, idx, &typeU8)) {
      outSamples->clear();
      if (err) *err = "truncated ghost record";
      return false;
    }

    const V2RecordType type = static_cast<V2RecordType>(typeU8);

    if (type == V2RecordType::Keyframe) {
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
        if (err) *err = "truncated ghost keyframe";
        return false;
      }

      const double tSec = std::isfinite((double)t) ? std::max(0.0, (double)t) : 0.0;
      if (hasPrevT && tSec + 1e-9 < prevT) {
        outSamples->clear();
        if (err) *err = "non-monotonic ghost time";
        return false;
      }
      prevT = tSec;
      hasPrevT = true;

      TimeTrialGhostSample s{};
      s.tSec = tSec;
      s.posKm = {
          std::isfinite((double)px) ? (double)px : 0.0,
          std::isfinite((double)py) ? (double)py : 0.0,
          std::isfinite((double)pz) ? (double)pz : 0.0,
      };

      s.orient = sanitizeQuat(math::Quatd{(double)qw, (double)qx, (double)qy, (double)qz});

      prevPosKm = s.posKm;
      hasPrevPos = true;

      outSamples->push_back(s);
    } else if (type == V2RecordType::Delta) {
      if (!hasPrevPos) {
        outSamples->clear();
        if (err) *err = "delta record before first keyframe";
        return false;
      }

      float t = 0.0f;
      core::i16 dxQ = 0, dyQ = 0, dzQ = 0;
      core::u8 qIdx = 0;
      core::i16 qa = 0, qb = 0, qc = 0;

      if (!readF32LE(bytes, idx, &t) ||
          !readI16LE(bytes, idx, &dxQ) ||
          !readI16LE(bytes, idx, &dyQ) ||
          !readI16LE(bytes, idx, &dzQ) ||
          !readU8(bytes, idx, &qIdx) ||
          !readI16LE(bytes, idx, &qa) ||
          !readI16LE(bytes, idx, &qb) ||
          !readI16LE(bytes, idx, &qc)) {
        outSamples->clear();
        if (err) *err = "truncated ghost delta";
        return false;
      }

      const double tSec = std::isfinite((double)t) ? std::max(0.0, (double)t) : 0.0;
      if (hasPrevT && tSec + 1e-9 < prevT) {
        outSamples->clear();
        if (err) *err = "non-monotonic ghost time";
        return false;
      }
      prevT = tSec;
      hasPrevT = true;

      TimeTrialGhostSample s{};
      s.tSec = tSec;
      s.posKm = prevPosKm + math::Vec3d{(double)dxQ, (double)dyQ, (double)dzQ} * kPosQuantumKm;
      s.orient = unpackQuatSmallestThree(qIdx, qa, qb, qc);

      prevPosKm = s.posKm;
      outSamples->push_back(s);
    } else {
      outSamples->clear();
      if (err) *err = "unknown ghost record type";
      return false;
    }
  }

  return true;
}

} // namespace

std::string encodeTimeTrialGhostSamplesB64(const std::vector<TimeTrialGhostSample>& samples) {
  if (samples.empty()) return {};

  const core::u32 n = static_cast<core::u32>(std::min<std::size_t>(samples.size(), kMaxSamples));

  std::string bytes;
  bytes.reserve(4 + 1 + 4 + (std::size_t)n * (kBytesPerV2Delta + 2));

  bytes.append(kMagic, kMagic + 4);
  writeU8(bytes, kWireVersionCurrent);
  writeU32LE(bytes, n);

  math::Vec3d predPosKm{0.0, 0.0, 0.0};
  bool hasPredPos = false;
  double prevT = 0.0;

  for (core::u32 i = 0; i < n; ++i) {
    const auto& s = samples[(std::size_t)i];

    double tSec = s.tSec;
    if (!finiteD(tSec)) tSec = 0.0;
    tSec = std::clamp(tSec, 0.0, 1e9);
    if (i > 0 && tSec < prevT) tSec = prevT; // enforce monotonicity

    const float t = static_cast<float>(tSec);
    const math::Vec3d posKm = sanitizePos(s.posKm);
    const math::Quatd q = sanitizeQuat(s.orient);

    // First sample must be a keyframe so the stream is self-contained.
    if (i == 0 || !hasPredPos) {
      writeU8(bytes, static_cast<core::u8>(V2RecordType::Keyframe));
      writeFullSampleF32(bytes, t, posKm, q);
      predPosKm = posKm;
      hasPredPos = true;
      prevT = tSec;
      continue;
    }

    // Delta record: quantize position delta from the encoder's *predicted* decoded position.
    // This keeps error bounded (<= ~0.5 quantum per component) without requiring periodic keyframes.
    const math::Vec3d deltaKm = posKm - predPosKm;
    const long long dxQll = std::llround(deltaKm.x * kInvPosQuantumKm);
    const long long dyQll = std::llround(deltaKm.y * kInvPosQuantumKm);
    const long long dzQll = std::llround(deltaKm.z * kInvPosQuantumKm);

    const bool deltaInRange =
        (dxQll >= -kPosDeltaMaxQ && dxQll <= kPosDeltaMaxQ) &&
        (dyQll >= -kPosDeltaMaxQ && dyQll <= kPosDeltaMaxQ) &&
        (dzQll >= -kPosDeltaMaxQ && dzQll <= kPosDeltaMaxQ);

    if (!deltaInRange) {
      // Large jump: fall back to full keyframe.
      writeU8(bytes, static_cast<core::u8>(V2RecordType::Keyframe));
      writeFullSampleF32(bytes, t, posKm, q);
      predPosKm = posKm;
      prevT = tSec;
      continue;
    }

    const core::i16 dxQ = static_cast<core::i16>(dxQll);
    const core::i16 dyQ = static_cast<core::i16>(dyQll);
    const core::i16 dzQ = static_cast<core::i16>(dzQll);

    // Update predicted decoded position.
    predPosKm = predPosKm + math::Vec3d{(double)dxQ, (double)dyQ, (double)dzQ} * kPosQuantumKm;

    const PackedQuatSmallestThree pq = packQuatSmallestThree(q);

    writeU8(bytes, static_cast<core::u8>(V2RecordType::Delta));
    writeF32LE(bytes, t);
    writeI16LE(bytes, dxQ);
    writeI16LE(bytes, dyQ);
    writeI16LE(bytes, dzQ);
    writeU8(bytes, pq.largestIndex);
    writeI16LE(bytes, pq.a);
    writeI16LE(bytes, pq.b);
    writeI16LE(bytes, pq.c);

    prevT = tSec;
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

  core::u32 n = 0;
  if (!readU32LE(bytes, idx, &n)) {
    if (err) *err = "truncated ghost count";
    return false;
  }

  if (n > kMaxSamples) {
    if (err) *err = "ghost too large";
    return false;
  }

  if (ver == kWireVersionV1) {
    return decodeV1(bytes, idx, n, outSamples, err);
  }
  if (ver == kWireVersionV2) {
    return decodeV2(bytes, idx, n, outSamples, err);
  }

  if (err) *err = "unsupported ghost version";
  return false;
}

} // namespace stellar::sim
