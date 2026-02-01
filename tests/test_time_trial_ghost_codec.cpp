#include "stellar/sim/SaveGame.h"
#include "stellar/sim/TimeTrialGhostCodec.h"

#include <cmath>
#include <filesystem>
#include <iostream>

static bool approxEq(double a, double b, double eps) { return std::abs(a - b) <= eps; }

int test_time_trial_ghost_codec() {
  int fails = 0;

  // --- Codec round-trip (float32 quantization expected) ---
  std::vector<stellar::sim::TimeTrialGhostSample> in;
  in.reserve(64);
  for (int i = 0; i < 64; ++i) {
    stellar::sim::TimeTrialGhostSample s;
    s.tSec = i * 0.25;
    s.posKm = {1000.0 + i * 1.5, -50.0 + i * 0.1, 25.0 - i * 0.33};
    s.orient = stellar::math::Quatd::fromAxisAngle({0, 1, 0}, i * 0.01).normalized();
    in.push_back(s);
  }

  const std::string b64 = stellar::sim::encodeTimeTrialGhostSamplesB64(in);
  if (b64.empty()) {
    std::cerr << "[test_time_trial_ghost_codec] expected non-empty base64\n";
    ++fails;
  } else {
    std::vector<stellar::sim::TimeTrialGhostSample> out;
    std::string err;
    if (!stellar::sim::decodeTimeTrialGhostSamplesB64(b64, &out, &err)) {
      std::cerr << "[test_time_trial_ghost_codec] decode failed: " << err << "\n";
      ++fails;
    } else if (out.size() != in.size()) {
      std::cerr << "[test_time_trial_ghost_codec] size mismatch: " << out.size() << " vs " << in.size() << "\n";
      ++fails;
    } else {
      const double posEps = 1e-4;
      const double tEps = 1e-6;
      for (std::size_t i = 0; i < in.size(); ++i) {
        if (!approxEq(out[i].tSec, in[i].tSec, tEps)) {
          std::cerr << "[test_time_trial_ghost_codec] t mismatch at i=" << i << "\n";
          ++fails;
          break;
        }
        if ((out[i].posKm - in[i].posKm).length() > posEps) {
          std::cerr << "[test_time_trial_ghost_codec] pos mismatch at i=" << i << "\n";
          ++fails;
          break;
        }
        // Orientation may flip sign (q and -q are equivalent). Compare by absolute dot.
        const auto a = in[i].orient.normalized();
        const auto b = out[i].orient.normalized();
        const double dot = a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z;
        if (std::abs(dot) < 0.9999) {
          std::cerr << "[test_time_trial_ghost_codec] orient mismatch at i=" << i << " dot=" << dot << "\n";
          ++fails;
          break;
        }
      }
    }
  }

  // --- SaveGame round-trip for time-trial best records (including ghost blob) ---
  {
    const std::filesystem::path p = "tmp_time_trial_best.sav";
    std::filesystem::remove(p);

    stellar::sim::SaveGame s;
    s.version = 36;
    stellar::sim::TimeTrialBestRecord r;
    r.courseKey = 123456789ULL;
    r.bestTimeSec = 42.5;
    r.ghostCodec = stellar::sim::kTimeTrialGhostCodecVersion;
    r.ghostB64 = b64;
    s.timeTrialBest.push_back(r);

    if (!stellar::sim::saveToFile(s, p.string())) {
      std::cerr << "[test_time_trial_ghost_codec] saveToFile failed\n";
      ++fails;
    } else {
      stellar::sim::SaveGame out;
      if (!stellar::sim::loadFromFile(p.string(), out)) {
        std::cerr << "[test_time_trial_ghost_codec] loadFromFile failed\n";
        ++fails;
      } else if (out.timeTrialBest.size() != 1) {
        std::cerr << "[test_time_trial_ghost_codec] expected 1 best record, got " << out.timeTrialBest.size() << "\n";
        ++fails;
      } else {
        const auto& rr = out.timeTrialBest[0];
        if (rr.courseKey != r.courseKey || !approxEq(rr.bestTimeSec, r.bestTimeSec, 1e-9) || rr.ghostCodec != r.ghostCodec) {
          std::cerr << "[test_time_trial_ghost_codec] best record mismatch\n";
          ++fails;
        }
        if (rr.ghostB64 != r.ghostB64) {
          std::cerr << "[test_time_trial_ghost_codec] ghost blob mismatch\n";
          ++fails;
        }
      }
    }

    std::filesystem::remove(p);
  }

  return fails;
}
