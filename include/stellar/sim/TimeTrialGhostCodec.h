#pragma once

#include "stellar/core/Types.h"
#include "stellar/math/Quat.h"
#include "stellar/math/Vec3.h"

#include <string>
#include <string_view>
#include <vector>

namespace stellar::sim {

// Compact ghost recording codec used by the Time Trials gameplay loop.
//
// Design goals:
//  - SaveGame-friendly: a single whitespace-free Base64 token per ghost.
//  - Deterministic and robust: explicit little-endian encoding + validation.
//  - Lightweight: stores only what is required to replay a ghost ship (pose over time).

// Current / maximum supported codec version.
// The decoder is backwards compatible with earlier versions.
static constexpr core::u32 kTimeTrialGhostCodecVersion = 2;

struct TimeTrialGhostSample {
  double tSec{0.0};
  math::Vec3d posKm{};
  math::Quatd orient{math::Quatd::identity()};
};

// Encodes samples into a Base64 token.
// Returns an empty string if `samples` is empty.
std::string encodeTimeTrialGhostSamplesB64(const std::vector<TimeTrialGhostSample>& samples);

// Decodes a Base64 token into `outSamples`.
// Returns false on malformed input; `outSamples` is cleared on failure.
// Treats empty string or "~" as "no ghost" and returns true.
bool decodeTimeTrialGhostSamplesB64(std::string_view b64,
                                   std::vector<TimeTrialGhostSample>* outSamples,
                                   std::string* err = nullptr);

} // namespace stellar::sim
